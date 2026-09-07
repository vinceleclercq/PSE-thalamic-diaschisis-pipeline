%% PSE_PROJECT_ROI_TO_MODALITIES_03.m
% PSE acute multimodal ROI projection + quantitative extraction
%
% Purpose
% -------
% Project the four validated FreeSurfer ROIs from native T1 space into the
% native voxel grids of the already-coregistered quantitative modalities:
%
%   PET_FDG
%   ASL_CBF
%   CTP_CBF
%   CTP_CBV
%   CTP_MTT
%   CTP_TMAX
%   CTP_TTP
%
% IMPORTANT METHODOLOGICAL POINT:
%   Quantitative PET/ASL/CTP images are NEVER resampled or modified.
%   Only the binary ROI masks are resampled to the target modality grid,
%   using nearest-neighbour interpolation in world coordinates.
%
% Inputs expected
% ---------------
% FreeSurfer native masks:
%   <PSE_ROOT>/derivatives/roi_masks_native/
%       sub-P00X/acute/
%           Left_Thalamus_mask.nii.gz
%           Right_Thalamus_mask.nii.gz
%           Left_Hippocampus_mask.nii.gz
%           Right_Hippocampus_mask.nii.gz
%
% Coregistered modality copies:
%   PET/ASL:
%     <PSE_ROOT>/derivatives/coreg/sub-P00X/acute/<MODALITY>/
%
%   CTP:
%     <PSE_ROOT>/derivatives/coreg_ctp_centered/
%         sub-P00X/acute/<MODALITY>/
%
% Outputs
% -------
% 1) Projected masks:
%   <PSE_ROOT>/derivatives/roi_masks_modalities/
%
% 2) Quantitative tables:
%   <PSE_ROOT>/derivatives/roi_values/
%
%   PSE_ACUTE_ROI_MODALITY_VALUES_LONG_v03.csv
%   PSE_ACUTE_ROI_MODALITY_SUMMARY_LR_v03.csv
%   PSE_ACUTE_ROI_MODALITY_STATUS_v03.csv
%   PSE_ACUTE_ROI_MODALITY_RESULTS_v03.mat
%
% Definitions
% -----------
% Ratio_LR = Left / Right
% AI_LR    = (Left - Right) / (Left + Right)
%
% These are descriptive LEFT/RIGHT metrics only.
% Ipsi/contra metrics should be calculated later after adding stroke side.
%
% Missing modalities are skipped and recorded in the status table.
%
% CTP validity rule:
%   negative CTP values are excluded from quantitative statistics because
%   they are non-physiologic fill/sentinel values in these exported maps.
%   Zeros are retained. Coverage after this exclusion is explicitly reported.
%
% Requires SPM (tested for SPM25-compatible calling syntax).

clear; clc;

%% ========================================================================
% Configuration

ROOT = string(getenv("PSE_ROOT"));
if strlength(ROOT)==0, ROOT = fullfile(string(getenv("HOME")),"PSE"); end

SPM_DIR = string(getenv("SPM_DIR"));
if strlength(SPM_DIR)>0
    SPM_DIR_CANDIDATES = SPM_DIR;
else
    SPM_DIR_CANDIDATES = [fullfile(string(getenv("HOME")),"spm25"), ...
                          fullfile(string(getenv("HOME")),"spm")];
end

ROI_NATIVE_ROOT = fullfile(ROOT, "derivatives", "roi_masks_native");
COREG_ROOT      = fullfile(ROOT, "derivatives", "coreg");
CTP_COREG_ROOT  = fullfile(ROOT, "derivatives", "coreg_ctp_centered");

OUT_MASK_ROOT   = fullfile(ROOT, "derivatives", "roi_masks_modalities");
OUT_VALUE_ROOT  = fullfile(ROOT, "derivatives", "roi_values");
WORK_ROOT       = fullfile(OUT_MASK_ROOT, "_work_unzipped");

if ~isfolder(OUT_MASK_ROOT), mkdir(OUT_MASK_ROOT); end
if ~isfolder(OUT_VALUE_ROOT), mkdir(OUT_VALUE_ROOT); end
if ~isfolder(WORK_ROOT), mkdir(WORK_ROOT); end

subjects = discoverSubjectDirs(ROI_NATIVE_ROOT);
if isempty(subjects)
    error("No subject ROI directories were found under: %s", ROI_NATIVE_ROOT);
end

roiNames = [ ...
    "Left_Thalamus", ...
    "Right_Thalamus", ...
    "Left_Hippocampus", ...
    "Right_Hippocampus"];

roiFiles = [ ...
    "Left_Thalamus_mask.nii.gz", ...
    "Right_Thalamus_mask.nii.gz", ...
    "Left_Hippocampus_mask.nii.gz", ...
    "Right_Hippocampus_mask.nii.gz"];

roiSides = ["Left","Right","Left","Right"];
roiStructures = ["Thalamus","Thalamus","Hippocampus","Hippocampus"];

% PET/ASL use the standard coregistration output.
% CTP maps use the centered CTP registration output.
modalities = struct( ...
    'name', { ...
        'PET_FDG', ...
        'ASL_CBF', ...
        'CTP_CBF', ...
        'CTP_CBV', ...
        'CTP_MTT', ...
        'CTP_TMAX', ...
        'CTP_TTP'}, ...
    'rootType', { ...
        'coreg', ...
        'coreg', ...
        'ctp', ...
        'ctp', ...
        'ctp', ...
        'ctp', ...
        'ctp'} );

%% ========================================================================
% Initialize SPM

if exist('spm_vol', 'file') ~= 2
    for iSPM = 1:numel(SPM_DIR_CANDIDATES)
        if isfolder(SPM_DIR_CANDIDATES(iSPM))
            addpath(char(SPM_DIR_CANDIDATES(iSPM)));
            break;
        end
    end
end

if exist('spm_vol', 'file') ~= 2
    error('SPM not found. Checked: %s and %s', ...
        SPM_DIR_CANDIDATES(1), SPM_DIR_CANDIDATES(2));
end

% IMPORTANT for SPM25 compatibility:
% SPM APIs expect character vectors here, not MATLAB string scalars.
spm('defaults', 'PET');

fprintf("============================================================\n");
fprintf("PSE ACUTE ROI -> PET / ASL / CTP\n");
fprintf("SPM: %s\n", spm('Ver'));
fprintf("============================================================\n");

%% ========================================================================
% Output row storage

longRows = {};
statusRows = {};

%% ========================================================================
% Main loop

for isub = 1:numel(subjects)

    SUB = subjects(isub);

    fprintf("\n============================================================\n");
    fprintf("[%d/%d] %s\n", isub, numel(subjects), SUB);
    fprintf("============================================================\n");

    nativeMaskDir = fullfile(ROI_NATIVE_ROOT, SUB, "acute");

    if ~isfolder(nativeMaskDir)
        warning("Native ROI directory missing for %s", SUB);
        continue;
    end

    % Prepare the four native ROI masks for SPM.
    nativeMaskSPM = strings(numel(roiFiles),1);

    for ir = 1:numel(roiFiles)
        p = fullfile(nativeMaskDir, roiFiles(ir));

        if ~isfile(p)
            error("Missing native ROI mask: %s", p);
        end

        nativeMaskSPM(ir) = ensureUnzippedForSPM( ...
            p, fullfile(WORK_ROOT, SUB, "native_masks"));
    end

    for im = 1:numel(modalities)

        MOD = string(modalities(im).name);
        rootType = string(modalities(im).rootType);

        if rootType == "coreg"
            modalityRoot = COREG_ROOT;
        else
            modalityRoot = CTP_COREG_ROOT;
        end

        [targetFile, findMessage] = findModalityNifti( ...
            modalityRoot, SUB, "acute", MOD);

        if strlength(targetFile) == 0
            fprintf("  %-10s : MISSING (%s)\n", MOD, findMessage);

            statusRows(end+1,:) = { ...
                char(SUB), char(MOD), "MISSING", "", char(findMessage)}; %#ok<AGROW>
            continue;
        end

        fprintf("\n  %-10s : %s\n", MOD, targetFile);

        targetWorkDir = fullfile(WORK_ROOT, SUB, MOD);
        targetSPM = ensureUnzippedForSPM(targetFile, targetWorkDir);

        VtAll = spm_vol(char(targetSPM));

        if numel(VtAll) ~= 1
            warning("%s %s has %d volumes. Using the first volume only.", ...
                SUB, MOD, numel(VtAll));
            Vt = VtAll(1);
        else
            Vt = VtAll;
        end

        % Read quantitative data through SPM so NIfTI scaling is respected.
        Y = spm_read_vols(Vt);

        if ndims(Y) ~= 3
            error("Expected a 3D quantitative image for %s %s.", SUB, MOD);
        end

        targetVoxelVol = abs(det(Vt.mat(1:3,1:3)));

        modOutDir = fullfile(OUT_MASK_ROOT, SUB, "acute", MOD);
        if ~isfolder(modOutDir), mkdir(modOutDir); end

        for ir = 1:numel(roiNames)

            ROI = roiNames(ir);

            VmAll = spm_vol(char(nativeMaskSPM(ir)));
            if numel(VmAll) ~= 1
                error("Native ROI mask is not 3D: %s", nativeMaskSPM(ir));
            end
            Vm = VmAll;

            % -------------------------------------------------------------
            % Project ROI mask to target modality grid.
            % Nearest-neighbour sampling only.
            projectedMask = resampleMaskToTargetGrid(Vm, Vt);

            outMaskFile = fullfile( ...
                modOutDir, ...
                ROI + "_mask_in_" + MOD + ".nii");

            writeBinaryMaskLikeTarget(Vt, projectedMask, outMaskFile, ...
                sprintf("%s mask projected to %s grid", ROI, MOD));

            % -------------------------------------------------------------
            % Extract quantitative values without modifying target image.
            maskLogical = projectedMask > 0;

            nROI = nnz(maskLogical);

            valuesAll = Y(maskLogical);

            finiteMask = isfinite(valuesAll);
            nFinite = nnz(finiteMask);
            finiteFraction = safeDivide(nFinite, nROI);

            % CTP parametric maps should not contain physiologically valid
            % negative CBF/CBV/time values. In these data, repeated negative
            % fill/sentinel values (e.g. -102 / -1024-like values) can be
            % present outside the valid perfusion-map support and severely
            % bias ordinary means. Therefore:
            %   PET / ASL : all finite ROI values are retained.
            %   CTP       : finite values >= 0 are retained for statistics.
            %
            % We DO retain zeros; only negative CTP values are excluded.
            if startsWith(MOD, "CTP_")
                statsMask = finiteMask & (valuesAll >= 0);
            else
                statsMask = finiteMask;
            end

            negativeMask = finiteMask & (valuesAll < 0);

            nStats = nnz(statsMask);
            statsCoverageFraction = safeDivide(nStats, nROI);

            nNegativeExcluded = nnz(negativeMask);
            negativeExcludedFraction = safeDivide(nNegativeExcluded, nROI);

            values = valuesAll(statsMask);

            rawFiniteValues = valuesAll(finiteMask);

            if isempty(rawFiniteValues)
                rawMin = NaN;
                rawMax = NaN;
            else
                rawMin = min(rawFiniteValues);
                rawMax = max(rawFiniteValues);
            end

            if isempty(values)
                vMean = NaN;
                vMedian = NaN;
                vSD = NaN;
                vMin = NaN;
                vMax = NaN;
            else
                vMean = mean(values);
                vMedian = median(values);
                vSD = std(values);
                vMin = min(values);
                vMax = max(values);
            end

            projectedMaskVolume = nROI * targetVoxelVol;

            if nROI == 0
                qcFlag = "FAIL_EMPTY_MASK";
            elseif nROI < 5
                qcFlag = "WARN_VERY_SMALL_MASK";
            elseif startsWith(MOD, "CTP_") && statsCoverageFraction < 0.50
                qcFlag = "FAIL_CTP_VALID_COVERAGE_LT50";
            elseif startsWith(MOD, "CTP_") && statsCoverageFraction < 0.80
                qcFlag = "WARN_CTP_VALID_COVERAGE_LT80";
            elseif ~startsWith(MOD, "CTP_") && finiteFraction < 0.80
                qcFlag = "WARN_FINITE_COVERAGE_LT80";
            else
                qcFlag = "OK_NUMERIC";
            end

            fprintf("    %-20s N=%5d  finite=%6.1f%%  stats=%6.1f%%  neg_excl=%6.1f%%  mean=%g\n", ...
                ROI, nROI, 100*finiteFraction, 100*statsCoverageFraction, ...
                100*negativeExcludedFraction, vMean);

            longRows(end+1,:) = { ... %#ok<AGROW>
                char(SUB), ...
                char(MOD), ...
                char(ROI), ...
                char(roiStructures(ir)), ...
                char(roiSides(ir)), ...
                char(targetFile), ...
                char(outMaskFile), ...
                nROI, ...
                nFinite, ...
                finiteFraction, ...
                nStats, ...
                statsCoverageFraction, ...
                nNegativeExcluded, ...
                negativeExcludedFraction, ...
                targetVoxelVol, ...
                projectedMaskVolume, ...
                vMean, ...
                vMedian, ...
                vSD, ...
                vMin, ...
                vMax, ...
                rawMin, ...
                rawMax, ...
                char(qcFlag)};
        end

        statusRows(end+1,:) = { ...
            char(SUB), char(MOD), "PROCESSED", char(targetFile), "OK"}; %#ok<AGROW>
    end
end

%% ========================================================================
% Long-format table

longVarNames = { ...
    'Subject', ...
    'Modality', ...
    'ROI', ...
    'Structure', ...
    'Side', ...
    'TargetFile', ...
    'ProjectedMaskFile', ...
    'Nvox_ROI', ...
    'Nvox_Finite', ...
    'FiniteFraction', ...
    'Nvox_StatsEligible', ...
    'StatsCoverageFraction', ...
    'Nvox_NegativeExcluded', ...
    'NegativeExcludedFraction', ...
    'TargetVoxelVolume_mm3', ...
    'ProjectedMaskVolume_mm3', ...
    'Mean', ...
    'Median', ...
    'SD', ...
    'Min', ...
    'Max', ...
    'RawMin', ...
    'RawMax', ...
    'NumericQC'};

if isempty(longRows)
    error("No modality/ROI values were extracted.");
end

Tlong = cell2table(longRows, 'VariableNames', longVarNames);

%% ========================================================================
% Build paired Left/Right summary table

summaryRows = {};

uniqueSubjects = unique(string(Tlong.Subject), "stable");
uniqueModalities = unique(string(Tlong.Modality), "stable");
structures = ["Thalamus","Hippocampus"];

for is = 1:numel(uniqueSubjects)
    SUB = uniqueSubjects(is);

    for im = 1:numel(uniqueModalities)
        MOD = uniqueModalities(im);

        for it = 1:numel(structures)
            STR = structures(it);

            idxL = string(Tlong.Subject)==SUB & ...
                   string(Tlong.Modality)==MOD & ...
                   string(Tlong.Structure)==STR & ...
                   string(Tlong.Side)=="Left";

            idxR = string(Tlong.Subject)==SUB & ...
                   string(Tlong.Modality)==MOD & ...
                   string(Tlong.Structure)==STR & ...
                   string(Tlong.Side)=="Right";

            if nnz(idxL) ~= 1 || nnz(idxR) ~= 1
                continue;
            end

            Lmean = Tlong.Mean(idxL);
            Rmean = Tlong.Mean(idxR);
            Lmedian = Tlong.Median(idxL);
            Rmedian = Tlong.Median(idxR);

            Lcoverage = Tlong.StatsCoverageFraction(idxL);
            Rcoverage = Tlong.StatsCoverageFraction(idxR);

            if Lcoverage < 0.50 || Rcoverage < 0.50
                meanRatioLR = NaN;
                meanAILR = NaN;
                medianRatioLR = NaN;
                medianAILR = NaN;
                pairQC = "FAIL_COVERAGE_LT50";
            else
                meanRatioLR = safeDivide(Lmean, Rmean);
                meanAILR = safeAI(Lmean, Rmean);

                medianRatioLR = safeDivide(Lmedian, Rmedian);
                medianAILR = safeAI(Lmedian, Rmedian);

                if Lcoverage < 0.80 || Rcoverage < 0.80
                    pairQC = "WARN_COVERAGE_LT80";
                else
                    pairQC = "OK";
                end
            end

            summaryRows(end+1,:) = { ... %#ok<AGROW>
                char(SUB), ...
                char(MOD), ...
                char(STR), ...
                Lmean, ...
                Rmean, ...
                meanRatioLR, ...
                meanAILR, ...
                Lmedian, ...
                Rmedian, ...
                medianRatioLR, ...
                medianAILR, ...
                Tlong.Nvox_ROI(idxL), ...
                Tlong.Nvox_ROI(idxR), ...
                Lcoverage, ...
                Rcoverage, ...
                Tlong.NegativeExcludedFraction(idxL), ...
                Tlong.NegativeExcludedFraction(idxR), ...
                char(pairQC)};
        end
    end
end

summaryVarNames = { ...
    'Subject', ...
    'Modality', ...
    'Structure', ...
    'LeftMean', ...
    'RightMean', ...
    'Mean_Ratio_LR', ...
    'Mean_AI_LR', ...
    'LeftMedian', ...
    'RightMedian', ...
    'Median_Ratio_LR', ...
    'Median_AI_LR', ...
    'Left_Nvox', ...
    'Right_Nvox', ...
    'Left_StatsCoverageFraction', ...
    'Right_StatsCoverageFraction', ...
    'Left_NegativeExcludedFraction', ...
    'Right_NegativeExcludedFraction', ...
    'PairQC'};

Tsummary = cell2table(summaryRows, 'VariableNames', summaryVarNames);

%% ========================================================================
% Modality status table

statusVarNames = {'Subject','Modality','Status','TargetFile','Message'};
Tstatus = cell2table(statusRows, 'VariableNames', statusVarNames);

%% ========================================================================
% Save results

longCSV = fullfile(OUT_VALUE_ROOT, ...
    "PSE_ACUTE_ROI_MODALITY_VALUES_LONG_v03.csv");

summaryCSV = fullfile(OUT_VALUE_ROOT, ...
    "PSE_ACUTE_ROI_MODALITY_SUMMARY_LR_v03.csv");

statusCSV = fullfile(OUT_VALUE_ROOT, ...
    "PSE_ACUTE_ROI_MODALITY_STATUS_v03.csv");

matFile = fullfile(OUT_VALUE_ROOT, ...
    "PSE_ACUTE_ROI_MODALITY_RESULTS_v03.mat");

writetable(Tlong, longCSV);
writetable(Tsummary, summaryCSV);
writetable(Tstatus, statusCSV);

save(matFile, "Tlong", "Tsummary", "Tstatus");

fprintf("\n============================================================\n");
fprintf("MULTIMODAL ROI EXTRACTION COMPLETE\n");
fprintf("============================================================\n");
fprintf("Long table:    %s\n", longCSV);
fprintf("LR summary:    %s\n", summaryCSV);
fprintf("Status table:  %s\n", statusCSV);
fprintf("MAT results:   %s\n", matFile);
fprintf("Projected ROI masks:\n  %s\n", OUT_MASK_ROOT);
fprintf("============================================================\n\n");

disp(Tstatus);
disp(Tsummary);

%% ========================================================================
% Local functions
% ========================================================================

function outPath = ensureUnzippedForSPM(inPath, outDir)

    inPath = string(inPath);
    outDir = string(outDir);

    if ~isfolder(outDir)
        mkdir(outDir);
    end

    if endsWith(lower(inPath), ".nii.gz")

        [~, nameGz, extGz] = fileparts(inPath); %#ok<ASGLU>
        % nameGz is e.g. "mask.nii"
        outPath = fullfile(outDir, nameGz);

        if ~isfile(outPath)
            files = gunzip(inPath, outDir);

            if isempty(files)
                error("gunzip failed for %s", inPath);
            end

            outPath = string(files{1});
        end

    else
        outPath = inPath;
    end
end

function [targetFile, message] = findModalityNifti(rootDir, subject, timepoint, modality)

    targetFile = "";
    message = "";

    expectedDir = fullfile(rootDir, subject, timepoint, modality);

    if isfolder(expectedDir)
        candidates = recursiveNiftiFiles(expectedDir);
    else
        subjectDir = fullfile(rootDir, subject, timepoint);

        if ~isfolder(subjectDir)
            message = "subject/timepoint coreg directory not found";
            return;
        end

        allCandidates = recursiveNiftiFiles(subjectDir);

        % Fallback: exact modality token in the full path.
        candidates = allCandidates( ...
            contains(lower(allCandidates), lower(modality)) );
    end

    if isempty(candidates)
        message = "no NIfTI found";
        return;
    end

    % Exclude QC / derived masks / Q.Clear from primary PET search.
    bad = contains(lower(candidates), "qc") | ...
          contains(lower(candidates), "mask") | ...
          contains(lower(candidates), "resliced");

    if modality == "PET_FDG"
        bad = bad | contains(lower(candidates), "qclear");
    end

    candidates = candidates(~bad);

    if isempty(candidates)
        message = "only excluded/QC NIfTI files found";
        return;
    end

    % If more than one file remains, use the largest file on disk.
    % A converted 3D modality volume should normally be the largest NIfTI.
    if numel(candidates) > 1
        sizes = nan(numel(candidates),1);

        for i = 1:numel(candidates)
            d = dir(candidates(i));
            if ~isempty(d)
                sizes(i) = d(1).bytes;
            end
        end

        [~, ix] = max(sizes);
        chosen = candidates(ix);

        message = sprintf("multiple candidates (%d); largest selected", ...
            numel(candidates));
    else
        chosen = candidates(1);
        message = "single candidate";
    end

    targetFile = chosen;
end

function files = recursiveNiftiFiles(folder)

    d1 = dir(fullfile(folder, "**", "*.nii"));
    d2 = dir(fullfile(folder, "**", "*.nii.gz"));

    d = [d1; d2];

    files = strings(numel(d),1);

    for i = 1:numel(d)
        files(i) = string(fullfile(d(i).folder, d(i).name));
    end
end

function maskTarget = resampleMaskToTargetGrid(Vm, Vt)
% Nearest-neighbour sampling of a binary source mask onto the target grid.
%
% Both volumes are related through their NIfTI/SPM world-coordinate affine
% matrices. The target quantitative image itself is never changed.

    dims = Vt.dim(1:3);

    nx = dims(1);
    ny = dims(2);
    nz = dims(3);

    maskTarget = false(nx, ny, nz);

    [X, Y] = ndgrid(1:nx, 1:ny);

    invMaskMat = inv(Vm.mat);

    for z = 1:nz

        nxy = numel(X);

        targetVox = [ ...
            X(:)'; ...
            Y(:)'; ...
            repmat(z, 1, nxy); ...
            ones(1, nxy)];

        world = Vt.mat * targetVox;
        sourceVox = invMaskMat * world;

        sampled = spm_sample_vol( ...
            Vm, ...
            sourceVox(1,:), ...
            sourceVox(2,:), ...
            sourceVox(3,:), ...
            0);  % nearest neighbour

        sampled = reshape(sampled, nx, ny);

        maskTarget(:,:,z) = sampled >= 0.5;
    end
end

function writeBinaryMaskLikeTarget(Vt, mask, outFile, description)

    Vo = Vt;
    Vo.fname = char(outFile);
    Vo.dt = [spm_type('uint8') 0];
    Vo.pinfo = [1; 0; 0];
    Vo.descrip = char(description);

    if isfield(Vo, "n")
        Vo.n = [1 1];
    end

    spm_write_vol(Vo, double(mask));
end

function y = safeDivide(a,b)
    if isempty(a) || isempty(b) || ~isfinite(a) || ~isfinite(b) || b == 0
        y = NaN;
    else
        y = a / b;
    end
end

function y = safeAI(left,right)
    denom = left + right;

    if isempty(left) || isempty(right) || ...
       ~isfinite(left) || ~isfinite(right) || ...
       denom == 0
        y = NaN;
    else
        y = (left - right) / denom;
    end
end

function subjects = discoverSubjectDirs(parentDir)
% Discover pseudonymized subject directories (sub-*) without assuming a
% fixed cohort size or numbering scheme.

    d = dir(fullfile(parentDir, "sub-*"));
    d = d([d.isdir]);
    subjects = sort(string({d.name}));
end
