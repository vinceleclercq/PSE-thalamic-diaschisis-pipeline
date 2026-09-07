%% PSE_COREG_ESTIMATE_01.m
% PSE PET-MR project - rigid multimodal coregistration to T1 (acute).
%
% IMPORTANT:
%   This script performs ESTIMATE ONLY coregistration.
%   It does NOT resample PET, ASL or CTP quantitative images.
%
% Rationale:
%   Quantitative images are kept in their native voxel grids to avoid
%   interpolation of PET/ASL/CTP values. Only the NIfTI affine headers of
%   COPIES in derivatives/coreg are updated. Later, T1-derived ROI masks
%   will be resampled into each quantitative modality's native grid using
%   nearest-neighbour interpolation.
%
% Source NIfTIs under derivatives/nifti are never modified.
%
% Reference:
%   T1 native image for each subject.
%
% Registration:
%   PET_FDG -> T1              (independent rigid NMI registration)
%   ASL_CBF -> T1              (independent rigid NMI registration)
%   CTP_REF -> T1              (rigid NMI registration)
%      same transform applied to CTP_CBF / CBV / TMAX / MTT / TTP
%
% Outputs:
%   <PSE_ROOT>/derivatives/coreg/sub-Pxxx/acute/...
%   <PSE_ROOT>/derivatives/coreg_logs/
%       PSE_COREG_ESTIMATE_LOG.csv
%
% Requirements:
%   MATLAB + SPM

clear; clc;

%% ----------------------------- SETTINGS --------------------------------

rootDir = getenv('PSE_ROOT');
if isempty(rootDir), rootDir = fullfile(getenv('HOME'),'PSE'); end
sessionName = 'acute';

niftiRoot = fullfile(rootDir, 'derivatives', 'nifti');
coregRoot = fullfile(rootDir, 'derivatives', 'coreg');
logRoot   = fullfile(rootDir, 'derivatives', 'coreg_logs');

overwriteExisting = true;

ctpMaps = {'CTP_CBF','CTP_CBV','CTP_TMAX','CTP_MTT','CTP_TTP'};

%% ---------------------------- PREFLIGHT --------------------------------

if isempty(which('spm')) || isempty(which('spm_coreg')) || isempty(which('spm_get_space'))
    error('SPM is not available on the MATLAB path.');
end

if ~isfolder(niftiRoot)
    error('NIfTI root not found: %s', niftiRoot);
end

if ~isfolder(coregRoot), mkdir(coregRoot); end
if ~isfolder(logRoot), mkdir(logRoot); end

spm('defaults','fmri');
try
    spm_get_defaults('cmdline', true);
catch
end

flags = spm_get_defaults('coreg.estimate');
flags.cost_fun = 'nmi';
flags.graphics = 0;

fprintf('\n============================================================\n');
fprintf('PSE COREGISTRATION - ESTIMATE ONLY\n');
fprintf('Reference: native T1\n');
fprintf('No quantitative image will be resampled.\n');
fprintf('============================================================\n\n');

subjects = dir(fullfile(niftiRoot,'sub-P*'));
subjects = subjects([subjects.isdir]);

if isempty(subjects)
    error('No sub-P* folders found under %s', niftiRoot);
end

[~,ord] = sort({subjects.name});
subjects = subjects(ord);

rows = cell(0,1);
rowN = 0;

%% ------------------------------ LOOP ------------------------------------

for s = 1:numel(subjects)

    subject = subjects(s).name;
    inSession = fullfile(niftiRoot, subject, sessionName);
    outSession = fullfile(coregRoot, subject, sessionName);

    fprintf('\n------------------------------------------------------------\n');
    fprintf('%s\n', subject);
    fprintf('------------------------------------------------------------\n');

    if ~isfolder(inSession)
        fprintf('Missing input session folder.\n');
        continue;
    end

    if ~isfolder(outSession), mkdir(outSession); end

    t1 = findSingleNifti(fullfile(inSession,'T1'));

    if strlength(t1)==0
        error('%s: no unique T1 NIfTI found.', subject);
    end

    % ---------------- PET_FDG -> T1 ----------------
    [R, transform] = coregSingleModality( ...
        subject, sessionName, 'PET_FDG', t1, ...
        fullfile(inSession,'PET_FDG'), ...
        fullfile(outSession,'PET_FDG'), ...
        flags, overwriteExisting);

    rowN = rowN + 1;
    rows{rowN} = R; %#ok<SAGROW>

    if ~isempty(transform)
        save(fullfile(outSession,'PET_FDG','PET_FDG_to_T1_transform.mat'), ...
            '-struct','transform');
    end

    printRow(R);

    % ---------------- ASL_CBF -> T1 ----------------
    [R, transform] = coregSingleModality( ...
        subject, sessionName, 'ASL_CBF', t1, ...
        fullfile(inSession,'ASL_CBF'), ...
        fullfile(outSession,'ASL_CBF'), ...
        flags, overwriteExisting);

    rowN = rowN + 1;
    rows{rowN} = R; %#ok<SAGROW>

    if ~isempty(transform)
        save(fullfile(outSession,'ASL_CBF','ASL_CBF_to_T1_transform.mat'), ...
            '-struct','transform');
    end

    printRow(R);

    % ---------------- CTP_REF + maps -> T1 ----------------
    [R, transform] = coregCTPGroup( ...
        subject, sessionName, t1, inSession, outSession, ...
        ctpMaps, flags, overwriteExisting);

    rowN = rowN + 1;
    rows{rowN} = R; %#ok<SAGROW>

    if ~isempty(transform)
        ctpOut = fullfile(outSession,'CTP');
        if ~isfolder(ctpOut), mkdir(ctpOut); end
        save(fullfile(ctpOut,'CTP_to_T1_transform.mat'), ...
            '-struct','transform');
    end

    printRow(R);
end

%% ------------------------------ LOG -------------------------------------

Log = struct2table(vertcat(rows{:}));

csvPath = fullfile(logRoot,'PSE_COREG_ESTIMATE_LOG.csv');
matPath = fullfile(logRoot,'PSE_COREG_ESTIMATE_LOG.mat');

writetable(Log,csvPath);
save(matPath,'Log');

fprintf('\n============================================================\n');
fprintf('COREGISTRATION ESTIMATION COMPLETE\n');
fprintf('Coreg copies: %s\n', coregRoot);
fprintf('Log CSV     : %s\n', csvPath);
fprintf('============================================================\n\n');

disp(Log(:,{'Subject','Modality','Status', ...
    'Tx_mm','Ty_mm','Tz_mm','Rx_deg','Ry_deg','Rz_deg','Warning'}));

%% ============================= FUNCTIONS ================================

function [R, transform] = coregSingleModality(subject, sessionName, modality, ...
        t1, sourceFolder, outFolder, flags, overwriteExisting)

    R = initRow(subject,sessionName,modality);
    transform = [];

    src = findSingleNifti(sourceFolder);

    if strlength(src)==0
        R.Status = "MISSING";
        R.Warning = "No unique source NIfTI found.";
        return;
    end

    if ~isfolder(outFolder), mkdir(outFolder); end

    dst = fullfile(outFolder, modality + "_coreg.nii");

    try
        copyNiftiAndJson(src,dst,overwriteExisting);

        VG = spm_vol(char(t1));
        VF = spm_vol(char(dst));

        oldMat = spm_get_space(char(dst));

        x = spm_coreg(VG,VF,flags);
        M = spm_matrix(x);

        newMat = M \ oldMat;
        spm_get_space(char(dst),newMat);

        R = addParams(R,x);
        R.Status = "OK";
        R.Output = string(dst);

        warningText = transformMagnitudeWarning(x);
        if strlength(warningText)>0
            R.Status = "WARNING";
            R.Warning = warningText;
        end

        transform = struct( ...
            'subject', string(subject), ...
            'modality', string(modality), ...
            'reference', string(t1), ...
            'source', string(src), ...
            'output', string(dst), ...
            'x', x, ...
            'M_spm_matrix', M, ...
            'oldMat', oldMat, ...
            'newMat', newMat);

    catch ME
        R.Status = "ERROR";
        R.Warning = string(ME.message);
    end
end

function [R, transform] = coregCTPGroup(subject, sessionName, t1, ...
        inSession, outSession, ctpMaps, flags, overwriteExisting)

    R = initRow(subject,sessionName,'CTP_REF');
    transform = [];

    refSrc = findSingleNifti(fullfile(inSession,'CTP_REF'));

    if strlength(refSrc)==0
        R.Status = "MISSING";
        R.Warning = "No CTP_REF NIfTI found.";
        return;
    end

    ctpOut = fullfile(outSession,'CTP');
    if ~isfolder(ctpOut), mkdir(ctpOut); end

    refDst = fullfile(ctpOut,'CTP_REF_coreg.nii');

    try
        copyNiftiAndJson(refSrc,refDst,overwriteExisting);

        % Copy all available parametric maps BEFORE applying transformation.
        mapInfo = struct('modality',{},'src',{},'dst',{},'oldMat',{},'newMat',{});

        for k = 1:numel(ctpMaps)

            mod = ctpMaps{k};
            mapSrc = findSingleNifti(fullfile(inSession,mod));

            if strlength(mapSrc)==0
                continue;
            end

            mapDst = fullfile(ctpOut, mod + "_coreg.nii");
            copyNiftiAndJson(mapSrc,mapDst,overwriteExisting);

            mapInfo(end+1).modality = mod; %#ok<AGROW>
            mapInfo(end).src = mapSrc;
            mapInfo(end).dst = string(mapDst);
            mapInfo(end).oldMat = spm_get_space(char(mapDst));
            mapInfo(end).newMat = [];
        end

        VG = spm_vol(char(t1));
        VF = spm_vol(char(refDst));

        oldRefMat = spm_get_space(char(refDst));

        x = spm_coreg(VG,VF,flags);
        M = spm_matrix(x);

        % Apply the rigid transform to the CTP reference COPY.
        newRefMat = M \ oldRefMat;
        spm_get_space(char(refDst),newRefMat);

        % Apply exactly the same world-space transform to every CTP map COPY.
        for k = 1:numel(mapInfo)
            mapInfo(k).newMat = M \ mapInfo(k).oldMat;
            spm_get_space(char(mapInfo(k).dst),mapInfo(k).newMat);
        end

        R = addParams(R,x);
        R.Status = "OK";
        R.Output = string(refDst);

        warningText = transformMagnitudeWarning(x);
        if strlength(warningText)>0
            R.Status = "WARNING";
            R.Warning = warningText;
        end

        transform = struct( ...
            'subject', string(subject), ...
            'modality', "CTP", ...
            'reference', string(t1), ...
            'source', string(refSrc), ...
            'output', string(refDst), ...
            'x', x, ...
            'M_spm_matrix', M, ...
            'oldRefMat', oldRefMat, ...
            'newRefMat', newRefMat, ...
            'maps', mapInfo);

    catch ME
        R.Status = "ERROR";
        R.Warning = string(ME.message);
    end
end

function R = initRow(subject,sessionName,modality)
    R = struct( ...
        'Subject',string(subject), ...
        'Session',string(sessionName), ...
        'Modality',string(modality), ...
        'Status',"", ...
        'Tx_mm',NaN, ...
        'Ty_mm',NaN, ...
        'Tz_mm',NaN, ...
        'Rx_deg',NaN, ...
        'Ry_deg',NaN, ...
        'Rz_deg',NaN, ...
        'Output',"", ...
        'Warning',"");
end

function R = addParams(R,x)
    R.Tx_mm = x(1);
    R.Ty_mm = x(2);
    R.Tz_mm = x(3);
    R.Rx_deg = rad2deg(x(4));
    R.Ry_deg = rad2deg(x(5));
    R.Rz_deg = rad2deg(x(6));
end

function warningText = transformMagnitudeWarning(x)

    translation = norm(x(1:3));
    rotationDeg = norm(rad2deg(x(4:6)));

    warnings = strings(0,1);

    if translation > 80
        warnings(end+1,1) = sprintf( ...
            'Large translation estimate: %.1f mm; visual QC required.',translation); %#ok<AGROW>
    end

    if rotationDeg > 45
        warnings(end+1,1) = sprintf( ...
            'Large rotation estimate: %.1f deg; visual QC required.',rotationDeg); %#ok<AGROW>
    end

    if isempty(warnings)
        warningText = "";
    else
        warningText = strjoin(warnings,' | ');
    end
end

function nii = findSingleNifti(folder)

    nii = "";

    if ~isfolder(folder)
        return;
    end

    d = dir(fullfile(folder,'*.nii'));

    if numel(d) ~= 1
        return;
    end

    nii = string(fullfile(d(1).folder,d(1).name));
end

function copyNiftiAndJson(src,dst,overwriteExisting)

    src = string(src);
    dst = string(dst);

    if isfile(dst)
        if overwriteExisting
            delete(dst);
        else
            error('Output already exists: %s',dst);
        end
    end

    copyfile(src,dst);

    [srcFolder,srcBase,~] = fileparts(src);
    [dstFolder,dstBase,~] = fileparts(dst);

    srcJson = fullfile(srcFolder,srcBase + ".json");
    dstJson = fullfile(dstFolder,dstBase + ".json");

    if isfile(srcJson)
        if isfile(dstJson) && overwriteExisting
            delete(dstJson);
        end
        copyfile(srcJson,dstJson);
    end
end

function printRow(R)

    fprintf('%-12s : %-8s',R.Modality,R.Status);

    if isfinite(R.Tx_mm)
        fprintf(' | T=[%.1f %.1f %.1f] mm | R=[%.1f %.1f %.1f] deg', ...
            R.Tx_mm,R.Ty_mm,R.Tz_mm,R.Rx_deg,R.Ry_deg,R.Rz_deg);
    end

    if strlength(R.Warning)>0
        fprintf(' | %s',R.Warning);
    end

    fprintf('\n');
end
