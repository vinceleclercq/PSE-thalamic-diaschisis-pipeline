%% PSE_NIFTI_QC_01.m
% PSE PET-MR project - NIfTI quality-control audit, acute time point.
%
% Non-destructive. Reads generated NIfTI files and writes:
%   derivatives/qc/PSE_NIFTI_QC.csv
%   derivatives/qc/montages/*.png
%
% It checks:
%   - image dimensions
%   - voxel sizes
%   - datatype
%   - finite / non-zero voxels
%   - basic intensity statistics
%   - obvious empty/constant volumes
%
% It also creates simple orthogonal montage PNGs for visual inspection.
%
% Requirements:
%   MATLAB + Image Processing Toolbox

clear; clc;

rootDir = getenv('PSE_ROOT');
if isempty(rootDir), rootDir = fullfile(getenv('HOME'),'PSE'); end
niftiRoot = fullfile(rootDir, 'derivatives', 'nifti');
qcRoot = fullfile(rootDir, 'derivatives', 'qc');
montageDir = fullfile(qcRoot, 'montages');

if ~isfolder(niftiRoot)
    error('NIfTI folder not found: %s', niftiRoot);
end
if ~isfolder(qcRoot), mkdir(qcRoot); end
if ~isfolder(montageDir), mkdir(montageDir); end

subjects = dir(fullfile(niftiRoot, 'sub-P*'));
subjects = subjects([subjects.isdir]);

if isempty(subjects)
    error('No converted subject folders found under %s', niftiRoot);
end

% Modalities we want to inspect in the first acute QC pass.
modalities = { ...
    'T1','FLAIR','DWI','ADC','PET_FDG','PET_FDG_QCLEAR', ...
    'ASL_CBF','CTP_REF','CTP_CBF','CTP_CBV','CTP_TMAX','CTP_MTT','CTP_TTP'};

rows = cell(0,1);
rowN = 0;

fprintf('\n============================================================\n');
fprintf('PSE NIFTI QC - ACUTE\n');
fprintf('NIfTI root: %s\n', niftiRoot);
fprintf('============================================================\n\n');

for s = 1:numel(subjects)

    subject = subjects(s).name;
    sessionDir = fullfile(niftiRoot, subject, 'acute');

    fprintf('\n------------------------------------------------------------\n');
    fprintf('%s\n', subject);
    fprintf('------------------------------------------------------------\n');

    for m = 1:numel(modalities)

        modality = modalities{m};
        modDir = fullfile(sessionDir, modality);

        if ~isfolder(modDir)
            continue;
        end

        niiFiles = dir(fullfile(modDir, '*.nii'));

        if isempty(niiFiles)
            continue;
        end

        [~, ord] = sort({niiFiles.name});
        niiFiles = niiFiles(ord);

        for n = 1:numel(niiFiles)

            niiPath = fullfile(niiFiles(n).folder, niiFiles(n).name);
            rowN = rowN + 1;

            R = initialiseRow(subject, 'acute', modality, niiFiles(n).name, niiPath);

            try
                info = niftiinfo(niiPath);
                img = niftiread(info);

                R.Status = "OK";
                R.ImageSize = numericVectorToText(info.ImageSize);

                if isfield(info, 'PixelDimensions')
                    R.PixelDimensions = numericVectorToText(info.PixelDimensions);
                end

                if isfield(info, 'Datatype')
                    R.Datatype = string(info.Datatype);
                end

                x = double(img(:));
                finiteMask = isfinite(x);
                finiteVals = x(finiteMask);

                R.TotalVoxels = numel(x);
                R.FiniteVoxels = sum(finiteMask);
                R.NaNVoxels = sum(isnan(x));
                R.InfVoxels = sum(isinf(x));

                if isempty(finiteVals)
                    R.Status = "WARNING";
                    R.Warning = "No finite voxel values.";
                else
                    R.Min = min(finiteVals);
                    R.Max = max(finiteVals);
                    R.Mean = mean(finiteVals);
                    R.Median = median(finiteVals);
                    R.Std = std(finiteVals);
                    R.NonZeroVoxels = sum(finiteVals ~= 0);

                    if R.Max == R.Min
                        R.Status = "WARNING";
                        R.Warning = "Constant-valued image.";
                    elseif R.NonZeroVoxels == 0
                        R.Status = "WARNING";
                        R.Warning = "All finite voxels are zero.";
                    end
                end

                % Create montage for 3-D images, and use first volume for 4-D.
                img3 = first3DVolume(img);
                pngPath = fullfile(montageDir, ...
                    sprintf('%s_%s_%s.png', subject, modality, erase(niiFiles(n).name,'.nii')));

                makeOrthogonalMontage(img3, pngPath);
                R.MontagePNG = string(pngPath);

            catch ME
                R.Status = "ERROR";
                R.Warning = string(ME.message);
            end

            rows{rowN} = R; %#ok<SAGROW>

            fprintf('%-18s %-18s : %-8s | size=%s | range=[%g, %g]', ...
                modality, niiFiles(n).name, R.Status, R.ImageSize, R.Min, R.Max);

            if strlength(R.Warning) > 0
                fprintf(' | %s', R.Warning);
            end
            fprintf('\n');
        end
    end
end

if isempty(rows)
    error('No NIfTI files were found.');
end

T = struct2table(vertcat(rows{:}));

csvPath = fullfile(qcRoot, 'PSE_NIFTI_QC.csv');
matPath = fullfile(qcRoot, 'PSE_NIFTI_QC.mat');

writetable(T, csvPath);
save(matPath, 'T');

fprintf('\n============================================================\n');
fprintf('QC COMPLETE\n');
fprintf('CSV      : %s\n', csvPath);
fprintf('Montages : %s\n', montageDir);
fprintf('============================================================\n\n');

disp(T(:, {'Subject','Modality','File','Status','ImageSize', ...
    'PixelDimensions','Min','Max','Mean','Warning'}));

%% ============================ FUNCTIONS =================================

function R = initialiseRow(subject, session, modality, fileName, filePath)
    R = struct( ...
        'Subject', string(subject), ...
        'Session', string(session), ...
        'Modality', string(modality), ...
        'File', string(fileName), ...
        'Path', string(filePath), ...
        'Status', "", ...
        'ImageSize', "", ...
        'PixelDimensions', "", ...
        'Datatype', "", ...
        'TotalVoxels', 0, ...
        'FiniteVoxels', 0, ...
        'NaNVoxels', 0, ...
        'InfVoxels', 0, ...
        'NonZeroVoxels', 0, ...
        'Min', NaN, ...
        'Max', NaN, ...
        'Mean', NaN, ...
        'Median', NaN, ...
        'Std', NaN, ...
        'MontagePNG', "", ...
        'Warning', "");
end

function txt = numericVectorToText(v)
    txt = strjoin(string(double(v(:)')), ' x ');
end

function img3 = first3DVolume(img)
    dims = ndims(img);
    if dims <= 3
        img3 = img;
    else
        idx = repmat({':'},1,dims);
        for d = 4:dims
            idx{d} = 1;
        end
        img3 = squeeze(img(idx{:}));
    end
end

function makeOrthogonalMontage(img, pngPath)

    img = double(img);

    if ndims(img) ~= 3 || any(size(img) < 2)
        return;
    end

    finiteVals = img(isfinite(img));
    if isempty(finiteVals)
        return;
    end

    % Robust display window, resistant to isolated extreme values.
    lo = percentileLocal(finiteVals, 1);
    hi = percentileLocal(finiteVals, 99);

    if ~isfinite(lo) || ~isfinite(hi) || hi <= lo
        lo = min(finiteVals);
        hi = max(finiteVals);
    end

    if hi <= lo
        return;
    end

    img(~isfinite(img)) = lo;
    img = (img - lo) ./ (hi - lo);
    img = min(max(img,0),1);

    sx = size(img,1);
    sy = size(img,2);
    sz = size(img,3);

    sagittal = squeeze(img(round(sx/2),:,:))';
    coronal  = squeeze(img(:,round(sy/2),:))';
    axial    = squeeze(img(:,:,round(sz/2)))';

    sagittal = flipud(sagittal);
    coronal  = flipud(coronal);
    axial    = flipud(axial);

    targetH = max([size(sagittal,1), size(coronal,1), size(axial,1)]);

    sagittal = resizeToHeight(sagittal, targetH);
    coronal  = resizeToHeight(coronal, targetH);
    axial    = resizeToHeight(axial, targetH);

    spacer = zeros(targetH, 8);
    canvas = [sagittal, spacer, coronal, spacer, axial];

    imwrite(canvas, pngPath);
end

function out = resizeToHeight(img, targetH)
    if size(img,1) == targetH
        out = img;
    else
        scale = targetH / size(img,1);
        out = imresize(img, scale, 'nearest');
    end
end

function p = percentileLocal(x, pct)
    x = sort(x(:));
    n = numel(x);
    if n == 0
        p = NaN;
        return;
    end
    idx = 1 + (n-1) * pct / 100;
    lo = floor(idx);
    hi = ceil(idx);
    if lo == hi
        p = x(lo);
    else
        p = x(lo) + (idx-lo) * (x(hi)-x(lo));
    end
end
