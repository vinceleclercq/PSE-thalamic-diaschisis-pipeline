%% 01_PSE_DICOM_AUDIT.m
% Non-destructive audit of PACS exports for the PSE PET-MR project.
% Reads DICOM metadata only; does not move, rename, convert, or modify source files.

clear; clc;

rootDir = getenv('PSE_ROOT');
if isempty(rootDir), rootDir = fullfile(getenv('HOME'),'PSE'); end

if ~isfolder(rootDir)
    error('Dataset folder not found: %s', rootDir);
end

if isempty(which('dicominfo'))
    error(['dicominfo was not found. Install Image Processing Toolbox, ' ...
           'restart MATLAB, then run: which dicominfo']);
end

outDir = fullfile(rootDir, 'derivatives', 'audit');
if ~isfolder(outDir)
    mkdir(outDir);
end

knownModalities = { ...
    'T1','DWI','ADC','FLAIR','PET_FDG','PET_FDG_QCLEAR', ...
    'ASL_CBF','ASL_SERIES_1','ASL_SERIES_2','ASL_UNSORTED','ASL_SOURCE', ...
    'CTP_REF','CTP_CBF','CTP_CBV','CTP_TMAX','CTP_MTT','CTP_TTP','rs_fMRI'};

sessionNames = {'acute','aigu','m03','m3','3m','m12','12m'};

fprintf('\n============================================================\n');
fprintf('PSE DICOM AUDIT\n');
fprintf('Root: %s\n', rootDir);
fprintf('============================================================\n\n');

sessionDirs = findSessionDirs(rootDir, knownModalities, sessionNames);

if isempty(sessionDirs)
    error(['No subject/session folders were detected. Check the folder structure ' ...
           'under %s.', rootDir]);
end

rows = cell(0,1);
rowN = 0;

for s = 1:numel(sessionDirs)

    sessionDir = sessionDirs(s).path;
    subject    = sessionDirs(s).subject;
    session    = sessionDirs(s).session;

    fprintf('\n------------------------------------------------------------\n');
    fprintf('SUBJECT: %s    SESSION: %s\n', subject, session);
    fprintf('Folder : %s\n', sessionDir);
    fprintf('------------------------------------------------------------\n');

    for m = 1:numel(knownModalities)

        modality = knownModalities{m};
        modDir = fullfile(sessionDir, modality);

        rowN = rowN + 1;
        R = initialiseRow(subject, session, modality, modDir);

        if ~isfolder(modDir)
            R.Status = "MISSING";
            R.Warning = "Folder absent";
            fprintf('%-18s : MISSING\n', modality);
            rows{rowN} = R; %#ok<SAGROW>
            continue;
        end

        files = recursiveFiles(modDir);

        if ~isempty(files)
            badNames = ismember(upper(string({files.name})), ...
                ["DICOMDIR","LOCKFILE","VERSION",".DS_STORE"]);
            files = files(~badNames);
        end

        R.TotalFiles = numel(files);

        if isempty(files)
            R.Status = "EMPTY";
            R.Warning = "Folder present but empty";
            fprintf('%-18s : EMPTY\n', modality);
            rows{rowN} = R; %#ok<SAGROW>
            continue;
        end

        dicomCount = 0;
        perfusionASLCount = 0;
        sourceASLCount = 0;
        m0ASLCount = 0;

        seriesDesc = strings(0,1);
        imageTypes = strings(0,1);
        modalities = strings(0,1);
        seriesNums = [];
        rowsVals = [];
        colsVals = [];
        photoInterp = strings(0,1);
        samplesPerPixel = [];
        slopes = [];
        intercepts = [];
        unitsVals = strings(0,1);

        for f = 1:numel(files)

            filename = fullfile(files(f).folder, files(f).name);

            try
                info = dicominfo(filename);
            catch
                continue;
            end

            dicomCount = dicomCount + 1;

            sd = getFieldString(info, 'SeriesDescription');
            it = getFieldString(info, 'ImageType');
            mo = getFieldString(info, 'Modality');
            pi = getFieldString(info, 'PhotometricInterpretation');
            un = getFieldString(info, 'Units');

            if strlength(sd) > 0, seriesDesc(end+1,1) = sd; end %#ok<SAGROW>
            if strlength(it) > 0, imageTypes(end+1,1) = it; end %#ok<SAGROW>
            if strlength(mo) > 0, modalities(end+1,1) = mo; end %#ok<SAGROW>
            if strlength(pi) > 0, photoInterp(end+1,1) = pi; end %#ok<SAGROW>
            if strlength(un) > 0, unitsVals(end+1,1) = un; end %#ok<SAGROW>

            if isfield(info,'SeriesNumber')
                seriesNums(end+1) = double(info.SeriesNumber); %#ok<SAGROW>
            end
            if isfield(info,'Rows')
                rowsVals(end+1) = double(info.Rows); %#ok<SAGROW>
            end
            if isfield(info,'Columns')
                colsVals(end+1) = double(info.Columns); %#ok<SAGROW>
            end
            if isfield(info,'SamplesPerPixel')
                samplesPerPixel(end+1) = double(info.SamplesPerPixel); %#ok<SAGROW>
            end
            if isfield(info,'RescaleSlope')
                slopes(end+1) = double(info.RescaleSlope); %#ok<SAGROW>
            end
            if isfield(info,'RescaleIntercept')
                intercepts(end+1) = double(info.RescaleIntercept); %#ok<SAGROW>
            end

            itUpper = upper(it);

            if contains(itUpper, 'PERFUSION_ASL')
                perfusionASLCount = perfusionASLCount + 1;
            elseif contains(itUpper, 'ASL')
                sourceASLCount = sourceASLCount + 1;
            end

            flatText = upper(flattenStructText(info, 2));
            if contains(flatText, 'M_ZERO_SCAN') || contains(flatText, 'M0_SCAN')
                m0ASLCount = m0ASLCount + 1;
            end
        end

        R.DICOMFiles = dicomCount;
        R.PerfusionASLFiles = perfusionASLCount;
        R.SourceASLFiles = sourceASLCount;
        R.M0ASLFiles = m0ASLCount;

        R.SeriesDescription = joinUnique(seriesDesc);
        R.ImageType = joinUnique(imageTypes);
        R.DICOMModality = joinUnique(modalities);
        R.SeriesNumber = joinUniqueNumeric(seriesNums);
        R.Rows = joinUniqueNumeric(rowsVals);
        R.Columns = joinUniqueNumeric(colsVals);
        R.PhotometricInterpretation = joinUnique(photoInterp);
        R.SamplesPerPixel = joinUniqueNumeric(samplesPerPixel);
        R.RescaleSlope = joinUniqueNumeric(slopes);
        R.RescaleIntercept = joinUniqueNumeric(intercepts);
        R.Units = joinUnique(unitsVals);

        [status, warning] = assessSeries(R);
        R.Status = status;
        R.Warning = warning;

        fprintf('%-18s : %-8s | DICOM=%d', modality, R.Status, R.DICOMFiles);

        if contains(upper(modality),'ASL')
            fprintf(' | PERFUSION_ASL=%d | source=%d', ...
                R.PerfusionASLFiles, R.SourceASLFiles);
        end

        if strlength(R.SeriesDescription) > 0
            fprintf(' | %s', R.SeriesDescription);
        end

        if strlength(R.Warning) > 0
            fprintf(' | %s', R.Warning);
        end

        fprintf('\n');

        rows{rowN} = R; %#ok<SAGROW>
    end
end

T = struct2table(vertcat(rows{:}));

csvPath = fullfile(outDir, 'PSE_DICOM_AUDIT.csv');
writetable(T, csvPath);

matPath = fullfile(outDir, 'PSE_DICOM_AUDIT.mat');
save(matPath, 'T');

fprintf('\n============================================================\n');
fprintf('AUDIT COMPLETE\n');
fprintf('CSV: %s\n', csvPath);
fprintf('MAT: %s\n', matPath);
fprintf('============================================================\n\n');

summaryCols = {'Subject','Session','Modality','Status','DICOMFiles', ...
    'PerfusionASLFiles','SourceASLFiles','SeriesDescription','Warning'};
disp(T(:,summaryCols));

%% ========================= LOCAL FUNCTIONS ==============================

function R = initialiseRow(subject, session, modality, folder)
    R = struct( ...
        'Subject', string(subject), ...
        'Session', string(session), ...
        'Modality', string(modality), ...
        'Folder', string(folder), ...
        'TotalFiles', 0, ...
        'DICOMFiles', 0, ...
        'DICOMModality', "", ...
        'SeriesDescription', "", ...
        'SeriesNumber', "", ...
        'ImageType', "", ...
        'Rows', "", ...
        'Columns', "", ...
        'PhotometricInterpretation', "", ...
        'SamplesPerPixel', "", ...
        'RescaleSlope', "", ...
        'RescaleIntercept', "", ...
        'Units', "", ...
        'PerfusionASLFiles', 0, ...
        'SourceASLFiles', 0, ...
        'M0ASLFiles', 0, ...
        'Status', "", ...
        'Warning', "");
end

function files = recursiveFiles(folder)
    d = dir(fullfile(folder, '**', '*'));
    files = d(~[d.isdir]);
end

function sessionDirs = findSessionDirs(rootDir, knownModalities, sessionNames)

    sessionDirs = struct('path',{},'subject',{},'session',{});

    allDirs = dir(fullfile(rootDir, '**', '*'));
    allDirs = allDirs([allDirs.isdir]);

    n = 0;

    for i = 1:numel(allDirs)
        nm = lower(allDirs(i).name);

        if ismember(nm, sessionNames)
            p = fullfile(allDirs(i).folder, allDirs(i).name);
            parent = fileparts(p);
            [~, subject] = fileparts(parent);

            n = n + 1;
            sessionDirs(n).path = p; %#ok<AGROW>
            sessionDirs(n).subject = subject;
            sessionDirs(n).session = allDirs(i).name;
        end
    end

    if ~isempty(sessionDirs)
        return;
    end

    top = dir(rootDir);
    top = top([top.isdir]);

    for i = 1:numel(top)

        if startsWith(top(i).name,'.') || strcmpi(top(i).name,'derivatives')
            continue;
        end

        p = fullfile(rootDir, top(i).name);

        hasModality = false;
        for m = 1:numel(knownModalities)
            if isfolder(fullfile(p, knownModalities{m}))
                hasModality = true;
                break;
            end
        end

        if hasModality
            n = n + 1;
            sessionDirs(n).path = p; %#ok<AGROW>
            sessionDirs(n).subject = top(i).name;
            sessionDirs(n).session = 'unspecified';
        end
    end
end

function s = getFieldString(info, fieldName)
    s = "";
    if ~isfield(info, fieldName)
        return;
    end

    v = info.(fieldName);

    if ischar(v)
        s = string(v);
    elseif isstring(v)
        s = join(v(:)', '\');
    elseif iscell(v)
        try
            s = join(string(v(:)'), '\');
        catch
            s = "";
        end
    elseif isnumeric(v) && isscalar(v)
        s = string(v);
    end
end

function s = joinUnique(vals)
    if isempty(vals)
        s = "";
        return;
    end
    vals = strip(vals);
    vals(vals=="") = [];
    vals = unique(vals,'stable');
    s = join(vals,' || ');
end

function s = joinUniqueNumeric(vals)
    if isempty(vals)
        s = "";
        return;
    end
    vals = unique(vals,'stable');
    s = join(string(vals),' || ');
end

function txt = flattenStructText(x, maxDepth)
    if nargin < 2
        maxDepth = 2;
    end

    parts = strings(0,1);

    if maxDepth < 0
        txt = "";
        return;
    end

    if isstruct(x)
        f = fieldnames(x);
        for k = 1:numel(f)
            try
                val = x.(f{k});
                if ischar(val) || isstring(val)
                    parts(end+1,1) = string(val); %#ok<AGROW>
                elseif isstruct(val) && maxDepth > 0
                    parts(end+1,1) = flattenStructText(val, maxDepth-1); %#ok<AGROW>
                end
            catch
            end
        end
    end

    txt = join(parts,' ');
end

function [status, warning] = assessSeries(R)

    status = "OK";
    warnings = strings(0,1);

    if R.DICOMFiles == 0
        status = "WARNING";
        warning = "No readable DICOM file found";
        return;
    end

    modName = upper(R.Modality);
    desc = upper(R.SeriesDescription);
    photo = upper(R.PhotometricInterpretation);

    if modName == "PET_FDG"
        if contains(desc,'NAC')
            warnings(end+1) = "PET_FDG appears non-attenuation-corrected (NAC)"; %#ok<AGROW>
        end
        if contains(desc,'QCLEAR') || contains(desc,'QC70')
            warnings(end+1) = "Main PET_FDG folder appears to contain Q.Clear reconstruction"; %#ok<AGROW>
        end
    end

    if contains(modName,'ASL')
        if R.PerfusionASLFiles > 0 && R.SourceASLFiles > 0
            warnings(end+1) = "Mixed perfusion and source ASL detected; automatic sorting required"; %#ok<AGROW>
        elseif modName == "ASL_CBF" && R.PerfusionASLFiles == 0
            warnings(end+1) = "ASL_CBF folder contains no PERFUSION_ASL image"; %#ok<AGROW>
        end
    end

    if startsWith(modName,'CTP_') && modName ~= "CTP_REF"
        sppText = split(R.SamplesPerPixel,' || ');
        spp = str2double(sppText);
        hasRGBSamples = any((spp == 3) & ~isnan(spp));

        if contains(photo,'RGB') || contains(photo,'PALETTE') || hasRGBSamples
            warnings(end+1) = "Possible colour-rendered perfusion map; verify quantitative scalar map"; %#ok<AGROW>
        end
    end

    if contains(R.SeriesDescription,' || ')
        warnings(end+1) = "Multiple SeriesDescription values found in the same folder"; %#ok<AGROW>
    end

    if contains(R.SeriesNumber,' || ')
        warnings(end+1) = "Multiple SeriesNumber values found in the same folder"; %#ok<AGROW>
    end

    if ~isempty(warnings)
        status = "WARNING";
        warning = join(unique(warnings,'stable'), ' | ');
    else
        warning = "";
    end
end
