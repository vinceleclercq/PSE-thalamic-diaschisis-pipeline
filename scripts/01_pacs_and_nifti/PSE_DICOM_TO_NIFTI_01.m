%% PSE_DICOM_TO_NIFTI_01.m
% PSE PET-MR project - DICOM to NIfTI conversion, acute time point only.
%
% This script:
%   1) reads DICOM files from <PSE_ROOT>/sub-Pxxx/acute
%   2) automatically selects only PERFUSION_ASL images for ASL_CBF
%   3) converts available series with SPM
%   4) NEVER modifies the source PACS export
%   5) writes only to <PSE_ROOT>/derivatives/nifti
%   6) writes a conversion log CSV
%
% Resting-state fMRI is deliberately NOT converted in this first pass.
%
% Requirements:
%   - MATLAB
%   - Image Processing Toolbox
%   - SPM on the MATLAB path

clear; clc;

%% ----------------------------- SETTINGS --------------------------------

rootDir = getenv('PSE_ROOT');
if isempty(rootDir), rootDir = fullfile(getenv('HOME'),'PSE'); end
sessionName = 'acute';

% Safety: if a modality already has generated NIfTI files, skip it.
overwriteExisting = false;

% rs-fMRI is not needed for the present ROI analysis and may create many
% 3D NIfTI volumes. Keep false for this first pass.
convertRestingState = false;

% Modalities converted directly from the folder with the same name.
directModalities = { ...
    'T1', ...
    'DWI', ...
    'ADC', ...
    'FLAIR', ...
    'PET_FDG', ...
    'PET_FDG_QCLEAR', ...
    'CTP_REF', ...
    'CTP_CBF', ...
    'CTP_CBV', ...
    'CTP_TMAX', ...
    'CTP_MTT', ...
    'CTP_TTP'};

% All possible locations in which PACS-exported ASL may be stored.
aslCandidateFolders = { ...
    'ASL_CBF', ...
    'ASL_UNSORTED', ...
    'ASL_SERIES_1', ...
    'ASL_SERIES_2', ...
    'ASL_SOURCE'};

%% ---------------------------- PREFLIGHT ---------------------------------

if ~isfolder(rootDir)
    error('Root folder not found: %s', rootDir);
end

if isempty(which('dicominfo'))
    error('dicominfo not found. Image Processing Toolbox is required.');
end

if isempty(which('spm')) || isempty(which('spm_dicom_headers')) || ...
        isempty(which('spm_dicom_convert'))
    error(['SPM is not available on the MATLAB path. Start MATLAB with SPM ' ...
           'configured, then run this script again.']);
end

spm('defaults','fmri');
try
    spm_get_defaults('cmdline', true);
catch
end

outRoot = fullfile(rootDir, 'derivatives', 'nifti');
logDir  = fullfile(rootDir, 'derivatives', 'conversion_logs');

if ~isfolder(outRoot), mkdir(outRoot); end
if ~isfolder(logDir),  mkdir(logDir);  end

fprintf('\n============================================================\n');
fprintf('PSE DICOM TO NIFTI - ACUTE DATA\n');
fprintf('Root: %s\n', rootDir);
fprintf('SPM : %s\n', spm('Version'));
fprintf('============================================================\n\n');

%% -------------------------- FIND SUBJECTS -------------------------------

subjects = dir(fullfile(rootDir, 'sub-P*'));
subjects = subjects([subjects.isdir]);

if isempty(subjects)
    error('No subject folders matching sub-P* were found under %s.', rootDir);
end

[~, order] = sort({subjects.name});
subjects = subjects(order);

%% ---------------------------- LOG TABLE ---------------------------------

logRows = cell(0,1);
logN = 0;

%% -------------------------- MAIN CONVERSION -----------------------------

for s = 1:numel(subjects)

    subject = subjects(s).name;
    sessionDir = fullfile(rootDir, subject, sessionName);

    fprintf('\n------------------------------------------------------------\n');
    fprintf('SUBJECT: %s    SESSION: %s\n', subject, sessionName);
    fprintf('------------------------------------------------------------\n');

    if ~isfolder(sessionDir)
        fprintf('Session folder missing: %s\n', sessionDir);
        continue;
    end

    % ---------- Standard modalities ----------
    for m = 1:numel(directModalities)

        modality = directModalities{m};
        sourceDir = fullfile(sessionDir, modality);
        outDir = fullfile(outRoot, subject, sessionName, modality);

        [status, nDicom, nNifti, generated, warningText] = ...
            convertFolderWithSPM(sourceDir, outDir, modality, overwriteExisting);

        logN = logN + 1;
        logRows{logN} = makeLogRow(subject, sessionName, modality, ...
            sourceDir, outDir, status, nDicom, nNifti, generated, warningText);

        printResult(modality, status, nDicom, nNifti, warningText);
    end

    % ---------- ASL: metadata-based automatic selection ----------
    modality = 'ASL_CBF';
    outDir = fullfile(outRoot, subject, sessionName, modality);

    [aslFiles, aslSourceFolders, aslWarning] = ...
        collectPerfusionASL(sessionDir, aslCandidateFolders);

    nDicom = numel(aslFiles);

    if isempty(aslFiles)
        status = "MISSING";
        nNifti = 0;
        generated = "";
        warningText = aslWarning;
        if strlength(warningText) == 0
            warningText = "No DICOM with ImageType containing PERFUSION_ASL was found.";
        end
    else
        [status, nNifti, generated, warningText2] = ...
            convertSelectedFilesWithSPM(aslFiles, outDir, modality, overwriteExisting);

        warningText = strtrim(joinNonEmpty([aslWarning, warningText2], " | "));
    end

    sourceSummary = strjoin(aslSourceFolders, ' || ');
    if isempty(sourceSummary)
        sourceSummary = fullfile(sessionDir, 'ASL_*');
    end

    logN = logN + 1;
    logRows{logN} = makeLogRow(subject, sessionName, modality, ...
        sourceSummary, outDir, status, nDicom, nNifti, generated, warningText);

    printResult(modality, status, nDicom, nNifti, warningText);

    % ---------- Optional rs-fMRI ----------
    if convertRestingState
        modality = 'rs_fMRI';
        sourceDir = fullfile(sessionDir, modality);
        outDir = fullfile(outRoot, subject, sessionName, modality);

        [status, nDicom, nNifti, generated, warningText] = ...
            convertFolderWithSPM(sourceDir, outDir, modality, overwriteExisting);

        logN = logN + 1;
        logRows{logN} = makeLogRow(subject, sessionName, modality, ...
            sourceDir, outDir, status, nDicom, nNifti, generated, warningText);

        printResult(modality, status, nDicom, nNifti, warningText);
    end
end

%% ----------------------------- SAVE LOG ---------------------------------

if isempty(logRows)
    error('No conversion attempts were logged.');
end

Log = struct2table(vertcat(logRows{:}));

csvPath = fullfile(logDir, 'PSE_NIFTI_CONVERSION_LOG.csv');
matPath = fullfile(logDir, 'PSE_NIFTI_CONVERSION_LOG.mat');

writetable(Log, csvPath);
save(matPath, 'Log');

fprintf('\n============================================================\n');
fprintf('CONVERSION COMPLETE\n');
fprintf('NIfTI root: %s\n', outRoot);
fprintf('Log CSV   : %s\n', csvPath);
fprintf('Log MAT   : %s\n', matPath);
fprintf('============================================================\n\n');

disp(Log(:, {'Subject','Session','Modality','Status', ...
    'SelectedDICOM','NIfTIProduced','Warning'}));

%% =========================== LOCAL FUNCTIONS ============================

function [status, nDicom, nNifti, generated, warningText] = ...
        convertFolderWithSPM(sourceDir, outDir, modality, overwriteExisting)

    nDicom = 0;
    nNifti = 0;
    generated = "";
    warningText = "";

    if ~isfolder(sourceDir)
        status = "MISSING";
        warningText = "Source folder absent.";
        return;
    end

    files = collectReadableDICOM(sourceDir);
    nDicom = numel(files);

    if isempty(files)
        status = "MISSING";
        warningText = "No readable DICOM found.";
        return;
    end

    [status, nNifti, generated, warningText] = ...
        convertSelectedFilesWithSPM(files, outDir, modality, overwriteExisting);
end

function [status, nNifti, generated, warningText] = ...
        convertSelectedFilesWithSPM(files, outDir, modality, overwriteExisting)

    nNifti = 0;
    generated = "";
    warningText = "";

    if isfolder(outDir)
        existing = dir(fullfile(outDir, '*.nii'));
    else
        existing = [];
    end

    if ~isempty(existing) && ~overwriteExisting
        status = "SKIPPED";
        nNifti = numel(existing);
        generated = strjoin(string({existing.name}), ' || ');
        warningText = "Output already exists; conversion skipped.";
        return;
    end

    if isfolder(outDir) && overwriteExisting
        delete(fullfile(outDir, '*.nii'));
        delete(fullfile(outDir, '*.json'));
    elseif ~isfolder(outDir)
        mkdir(outDir);
    end

    tmpDir = fullfile(outDir, '_spm_tmp');

    if isfolder(tmpDir)
        rmdir(tmpDir, 's');
    end
    mkdir(tmpDir);

    try
        fileChar = char(files(:));
        hdr = spm_dicom_headers(fileChar);

        out = spm_dicom_convert(hdr, 'all', 'flat', 'nii', tmpDir, true);

        produced = string(out.files(:));
        produced(produced == "") = [];

        % Keep only files that actually exist.
        existsMask = false(size(produced));
        for k = 1:numel(produced)
            existsMask(k) = isfile(produced(k));
        end
        produced = produced(existsMask);

        if isempty(produced)
            status = "ERROR";
            warningText = "SPM returned no NIfTI file.";
            rmdir(tmpDir, 's');
            return;
        end

        finalFiles = strings(numel(produced),1);

        for k = 1:numel(produced)

            oldNii = produced(k);

            if numel(produced) == 1
                newBase = string(modality);
            else
                newBase = sprintf('%s_%03d', modality, k);
            end

            newNii = fullfile(outDir, newBase + ".nii");

            if isfile(newNii)
                delete(newNii);
            end

            movefile(oldNii, newNii);
            finalFiles(k) = string(newNii);

            % If SPM created a JSON sidecar, rename it consistently.
            [oldFolder, oldBase, ~] = fileparts(oldNii);
            oldJson = fullfile(oldFolder, oldBase + ".json");
            newJson = fullfile(outDir, newBase + ".json");

            if isfile(oldJson)
                if isfile(newJson)
                    delete(newJson);
                end
                movefile(oldJson, newJson);
            end
        end

        % Remove temporary conversion folder after moving outputs.
        if isfolder(tmpDir)
            rmdir(tmpDir, 's');
        end

        nNifti = numel(finalFiles);
        generated = strjoin(finalFiles, ' || ');
        status = "OK";

        if nNifti > 1
            warningText = sprintf(['SPM produced %d NIfTI volumes. This is not ' ...
                'automatically an error, but it requires QC before downstream use.'], nNifti);
            status = "WARNING";
        end

    catch ME

        status = "ERROR";
        warningText = string(ME.message);

        if isfolder(tmpDir)
            try
                rmdir(tmpDir, 's');
            catch
            end
        end
    end
end

function files = collectReadableDICOM(folder)

    d = dir(fullfile(folder, '**', '*'));
    d = d(~[d.isdir]);

    if isempty(d)
        files = cell(0,1);
        return;
    end

    badNames = ismember(upper(string({d.name})), ...
        ["DICOMDIR","LOCKFILE","VERSION",".DS_STORE"]);
    d = d(~badNames);

    files = cell(0,1);

    for i = 1:numel(d)
        f = fullfile(d(i).folder, d(i).name);

        try
            dicominfo(f);
            files{end+1,1} = f; %#ok<AGROW>
        catch
        end
    end
end

function [files, sourceFolders, warningText] = ...
        collectPerfusionASL(sessionDir, candidateFolders)

    files = cell(0,1);
    sourceFolders = cell(0,1);
    warningText = "";

    sopUIDs = strings(0,1);
    nSourceASL = 0;
    nOtherDICOM = 0;

    for c = 1:numel(candidateFolders)

        folder = fullfile(sessionDir, candidateFolders{c});

        if ~isfolder(folder)
            continue;
        end

        candidates = collectReadableDICOM(folder);

        if isempty(candidates)
            continue;
        end

        sourceFolders{end+1,1} = folder; %#ok<AGROW>

        for i = 1:numel(candidates)

            f = candidates{i};

            try
                info = dicominfo(f);
            catch
                continue;
            end

            imageType = "";
            if isfield(info, 'ImageType')
                imageType = dicomValueToString(info.ImageType);
            end

            imageTypeUpper = upper(imageType);

            if contains(imageTypeUpper, 'PERFUSION_ASL')

                uid = "";
                if isfield(info, 'SOPInstanceUID')
                    uid = string(info.SOPInstanceUID);
                end

                if strlength(uid) > 0
                    if any(sopUIDs == uid)
                        continue;
                    end
                    sopUIDs(end+1,1) = uid; %#ok<AGROW>
                end

                files{end+1,1} = f; %#ok<AGROW>

            elseif contains(imageTypeUpper, 'ASL')
                nSourceASL = nSourceASL + 1;
            else
                nOtherDICOM = nOtherDICOM + 1;
            end
        end
    end

    if ~isempty(files)
        warningParts = strings(0,1);

        if nSourceASL > 0
            warningParts(end+1,1) = sprintf( ...
                '%d source ASL DICOM ignored; %d PERFUSION_ASL selected.', ...
                nSourceASL, numel(files)); %#ok<AGROW>
        end

        if nOtherDICOM > 0
            warningParts(end+1,1) = sprintf( ...
                '%d non-ASL DICOM in candidate folders ignored.', ...
                nOtherDICOM); %#ok<AGROW>
        end

        warningText = strjoin(warningParts, ' | ');
    end
end

function s = dicomValueToString(v)

    if ischar(v)
        s = string(v);
    elseif isstring(v)
        s = strjoin(v(:)', '\');
    elseif iscell(v)
        try
            s = strjoin(string(v(:)'), '\');
        catch
            s = "";
        end
    elseif isnumeric(v) && isscalar(v)
        s = string(v);
    else
        s = "";
    end
end

function R = makeLogRow(subject, sessionName, modality, sourceDir, outDir, ...
        status, nDicom, nNifti, generated, warningText)

    R = struct( ...
        'Subject', string(subject), ...
        'Session', string(sessionName), ...
        'Modality', string(modality), ...
        'Source', string(sourceDir), ...
        'OutputFolder', string(outDir), ...
        'Status', string(status), ...
        'SelectedDICOM', double(nDicom), ...
        'NIfTIProduced', double(nNifti), ...
        'GeneratedFiles', string(generated), ...
        'Warning', string(warningText));
end

function printResult(modality, status, nDicom, nNifti, warningText)

    fprintf('%-18s : %-8s | DICOM=%d | NIfTI=%d', ...
        modality, status, nDicom, nNifti);

    if strlength(string(warningText)) > 0
        fprintf(' | %s', string(warningText));
    end

    fprintf('\n');
end

function out = joinNonEmpty(vals, delimiter)

    vals = string(vals);
    vals = vals(strlength(strtrim(vals)) > 0);

    if isempty(vals)
        out = "";
    else
        out = strjoin(vals, delimiter);
    end
end
