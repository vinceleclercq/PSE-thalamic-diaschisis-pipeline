%% PSE_INTENSITY_SCALING_AUDIT_01.m
% PSE PET-MR project - audit of DICOM intensity scaling metadata.
%
% Purpose:
%   Before quantitative ROI extraction, verify how PET, ASL and CTP
%   pixel values are scaled in the original DICOM files.
%
% Non-destructive: reads metadata only.
%
% Outputs:
%   <PSE_ROOT>/derivatives/scaling_audit/
%       PSE_INTENSITY_SCALING_AUDIT.csv
%
% Requirements:
%   MATLAB + Image Processing Toolbox
%
% Notes:
%   - ASL: only DICOMs with ImageType containing PERFUSION_ASL are audited.
%   - PET: main PET_FDG folder only.
%   - CTP: CBF, CBV, Tmax, MTT, TTP.
%   - The script looks for standard DICOM RescaleSlope/Intercept and
%     RealWorldValueMappingSequence when present.

clear; clc;

rootDir = getenv('PSE_ROOT');
if isempty(rootDir), rootDir = fullfile(getenv('HOME'),'PSE'); end
sessionName = 'acute';

modalities = {'PET_FDG','ASL_CBF','CTP_CBF','CTP_CBV','CTP_TMAX','CTP_MTT','CTP_TTP'};

aslCandidateFolders = {'ASL_CBF','ASL_UNSORTED','ASL_SERIES_1','ASL_SERIES_2','ASL_SOURCE'};

outDir = fullfile(rootDir,'derivatives','scaling_audit');
if ~isfolder(outDir), mkdir(outDir); end

if isempty(which('dicominfo'))
    error('dicominfo not found. Image Processing Toolbox is required.');
end

subjects = dir(fullfile(rootDir,'sub-P*'));
subjects = subjects([subjects.isdir]);

if isempty(subjects)
    error('No subject folders matching sub-P* found in %s', rootDir);
end

[~,ord] = sort({subjects.name});
subjects = subjects(ord);

rows = cell(0,1);
rowN = 0;

fprintf('\n============================================================\n');
fprintf('PSE INTENSITY SCALING AUDIT\n');
fprintf('============================================================\n\n');

for s = 1:numel(subjects)

    subject = subjects(s).name;
    sessionDir = fullfile(rootDir, subject, sessionName);

    if ~isfolder(sessionDir)
        continue;
    end

    fprintf('\n------------------------------------------------------------\n');
    fprintf('%s\n', subject);
    fprintf('------------------------------------------------------------\n');

    for m = 1:numel(modalities)

        modality = modalities{m};

        if strcmp(modality,'ASL_CBF')
            [files, sourceSummary] = collectPerfusionASL(sessionDir, aslCandidateFolders);
        else
            sourceDir = fullfile(sessionDir, modality);
            files = collectReadableDICOM(sourceDir);
            sourceSummary = sourceDir;
        end

        rowN = rowN + 1;
        R = initialiseRow(subject, sessionName, modality, sourceSummary);

        if isempty(files)
            R.Status = "MISSING";
            R.Warning = "No selected DICOM files found.";
            rows{rowN} = R; %#ok<SAGROW>
            fprintf('%-12s : MISSING\n', modality);
            continue;
        end

        slopes = [];
        intercepts = [];
        rescaleTypes = strings(0,1);
        unitsVals = strings(0,1);
        seriesDesc = strings(0,1);
        imageTypes = strings(0,1);

        rwSlopes = [];
        rwIntercepts = [];
        rwUnits = strings(0,1);

        storedMin = [];
        storedMax = [];

        for f = 1:numel(files)

            try
                info = dicominfo(files{f});
            catch
                continue;
            end

            if isfield(info,'RescaleSlope')
                slopes(end+1) = double(info.RescaleSlope); %#ok<SAGROW>
            end
            if isfield(info,'RescaleIntercept')
                intercepts(end+1) = double(info.RescaleIntercept); %#ok<SAGROW>
            end
            if isfield(info,'RescaleType')
                rescaleTypes(end+1,1) = dicomValueToString(info.RescaleType); %#ok<SAGROW>
            end
            if isfield(info,'Units')
                unitsVals(end+1,1) = dicomValueToString(info.Units); %#ok<SAGROW>
            end
            if isfield(info,'SeriesDescription')
                seriesDesc(end+1,1) = dicomValueToString(info.SeriesDescription); %#ok<SAGROW>
            end
            if isfield(info,'ImageType')
                imageTypes(end+1,1) = dicomValueToString(info.ImageType); %#ok<SAGROW>
            end
            if isfield(info,'SmallestImagePixelValue')
                storedMin(end+1) = double(info.SmallestImagePixelValue); %#ok<SAGROW>
            end
            if isfield(info,'LargestImagePixelValue')
                storedMax(end+1) = double(info.LargestImagePixelValue); %#ok<SAGROW>
            end

            [a,b,u] = extractRealWorldMapping(info);
            if ~isnan(a), rwSlopes(end+1) = a; end %#ok<SAGROW>
            if ~isnan(b), rwIntercepts(end+1) = b; end %#ok<SAGROW>
            if strlength(u)>0, rwUnits(end+1,1) = u; end %#ok<SAGROW>
        end

        R.SelectedDICOM = numel(files);
        R.SeriesDescription = joinUnique(seriesDesc);
        R.ImageType = joinUnique(imageTypes);

        R.RescaleSlope = joinUniqueNumeric(slopes);
        R.RescaleIntercept = joinUniqueNumeric(intercepts);
        R.RescaleType = joinUnique(rescaleTypes);
        R.Units = joinUnique(unitsVals);

        R.RealWorldSlope = joinUniqueNumeric(rwSlopes);
        R.RealWorldIntercept = joinUniqueNumeric(rwIntercepts);
        R.RealWorldUnits = joinUnique(rwUnits);

        R.StoredPixelMin = joinUniqueNumeric(storedMin);
        R.StoredPixelMax = joinUniqueNumeric(storedMax);

        warnings = strings(0,1);

        if isempty(slopes)
            warnings(end+1) = "No standard RescaleSlope found."; %#ok<SAGROW>
        end
        if isempty(intercepts)
            warnings(end+1) = "No standard RescaleIntercept found."; %#ok<SAGROW>
        end
        if isempty(unitsVals) && isempty(rwUnits)
            warnings(end+1) = "No explicit units found in standard fields."; %#ok<SAGROW>
        end

        % Multiple scaling factors within one selected series deserve review.
        if numel(unique(slopes)) > 1
            warnings(end+1) = "Multiple RescaleSlope values within selected DICOMs."; %#ok<SAGROW>
        end
        if numel(unique(intercepts)) > 1
            warnings(end+1) = "Multiple RescaleIntercept values within selected DICOMs."; %#ok<SAGROW>
        end
        if numel(unique(rwSlopes)) > 1
            warnings(end+1) = "Multiple RealWorldValue slopes within selected DICOMs."; %#ok<SAGROW>
        end
        if numel(unique(rwIntercepts)) > 1
            warnings(end+1) = "Multiple RealWorldValue intercepts within selected DICOMs."; %#ok<SAGROW>
        end

        if isempty(warnings)
            R.Status = "OK";
            R.Warning = "";
        else
            R.Status = "CHECK";
            R.Warning = strjoin(unique(warnings,'stable'),' | ');
        end

        rows{rowN} = R; %#ok<SAGROW>

        fprintf('%-12s : %-6s | slope=%s | intercept=%s | units=%s', ...
            modality, R.Status, R.RescaleSlope, R.RescaleIntercept, R.Units);

        if strlength(R.RealWorldSlope)>0 || strlength(R.RealWorldUnits)>0
            fprintf(' | RW slope=%s | RW units=%s', R.RealWorldSlope, R.RealWorldUnits);
        end

        if strlength(R.Warning)>0
            fprintf(' | %s', R.Warning);
        end
        fprintf('\n');
    end
end

T = struct2table(vertcat(rows{:}));

csvPath = fullfile(outDir,'PSE_INTENSITY_SCALING_AUDIT.csv');
matPath = fullfile(outDir,'PSE_INTENSITY_SCALING_AUDIT.mat');

writetable(T,csvPath);
save(matPath,'T');

fprintf('\n============================================================\n');
fprintf('SCALING AUDIT COMPLETE\n');
fprintf('CSV: %s\n', csvPath);
fprintf('============================================================\n\n');

disp(T(:,{'Subject','Modality','Status','SelectedDICOM', ...
    'RescaleSlope','RescaleIntercept','Units', ...
    'RealWorldSlope','RealWorldIntercept','RealWorldUnits','Warning'}));

%% =========================== LOCAL FUNCTIONS ============================

function R = initialiseRow(subject, sessionName, modality, sourceSummary)
    R = struct( ...
        'Subject', string(subject), ...
        'Session', string(sessionName), ...
        'Modality', string(modality), ...
        'Source', string(sourceSummary), ...
        'Status', "", ...
        'SelectedDICOM', 0, ...
        'SeriesDescription', "", ...
        'ImageType', "", ...
        'RescaleSlope', "", ...
        'RescaleIntercept', "", ...
        'RescaleType', "", ...
        'Units', "", ...
        'RealWorldSlope', "", ...
        'RealWorldIntercept', "", ...
        'RealWorldUnits', "", ...
        'StoredPixelMin', "", ...
        'StoredPixelMax', "", ...
        'Warning', "");
end

function files = collectReadableDICOM(folder)

    files = cell(0,1);

    if ~isfolder(folder)
        return;
    end

    d = dir(fullfile(folder,'**','*'));
    d = d(~[d.isdir]);

    if isempty(d)
        return;
    end

    badNames = ismember(upper(string({d.name})), ...
        ["DICOMDIR","LOCKFILE","VERSION",".DS_STORE"]);
    d = d(~badNames);

    for i = 1:numel(d)
        f = fullfile(d(i).folder,d(i).name);
        try
            dicominfo(f);
            files{end+1,1} = f; %#ok<AGROW>
        catch
        end
    end
end

function [files, sourceSummary] = collectPerfusionASL(sessionDir, candidateFolders)

    files = cell(0,1);
    sourceFolders = strings(0,1);
    sopUIDs = strings(0,1);

    for c = 1:numel(candidateFolders)

        folder = fullfile(sessionDir,candidateFolders{c});

        if ~isfolder(folder)
            continue;
        end

        candidates = collectReadableDICOM(folder);

        if isempty(candidates)
            continue;
        end

        sourceFolders(end+1,1) = string(folder); %#ok<SAGROW>

        for i = 1:numel(candidates)

            f = candidates{i};

            try
                info = dicominfo(f);
            catch
                continue;
            end

            imageType = "";
            if isfield(info,'ImageType')
                imageType = upper(dicomValueToString(info.ImageType));
            end

            if ~contains(imageType,'PERFUSION_ASL')
                continue;
            end

            uid = "";
            if isfield(info,'SOPInstanceUID')
                uid = string(info.SOPInstanceUID);
            end

            if strlength(uid)>0
                if any(sopUIDs == uid)
                    continue;
                end
                sopUIDs(end+1,1) = uid; %#ok<SAGROW>
            end

            files{end+1,1} = f; %#ok<AGROW>
        end
    end

    sourceSummary = strjoin(sourceFolders,' || ');
end

function s = dicomValueToString(v)

    if ischar(v)
        s = string(v);
    elseif isstring(v)
        s = strjoin(v(:)','\');
    elseif iscell(v)
        try
            s = strjoin(string(v(:)'),'\');
        catch
            s = "";
        end
    elseif isnumeric(v) && isscalar(v)
        s = string(v);
    else
        s = "";
    end
end

function s = joinUnique(vals)

    if isempty(vals)
        s = "";
        return;
    end

    vals = string(vals);
    vals = strip(vals);
    vals(vals=="") = [];
    vals = unique(vals,'stable');

    s = strjoin(vals,' || ');
end

function s = joinUniqueNumeric(vals)

    if isempty(vals)
        s = "";
        return;
    end

    vals = unique(double(vals),'stable');
    s = strjoin(string(vals),' || ');
end

function [slope, intercept, units] = extractRealWorldMapping(info)

    slope = NaN;
    intercept = NaN;
    units = "";

    if ~isfield(info,'RealWorldValueMappingSequence')
        return;
    end

    seq = info.RealWorldValueMappingSequence;
    itemNames = fieldnames(seq);

    if isempty(itemNames)
        return;
    end

    item = seq.(itemNames{1});

    if isfield(item,'RealWorldValueSlope')
        slope = double(item.RealWorldValueSlope);
    end

    if isfield(item,'RealWorldValueIntercept')
        intercept = double(item.RealWorldValueIntercept);
    end

    if isfield(item,'MeasurementUnitsCodeSequence')
        mus = item.MeasurementUnitsCodeSequence;
        f = fieldnames(mus);
        if ~isempty(f)
            unitItem = mus.(f{1});
            if isfield(unitItem,'CodeMeaning')
                units = string(unitItem.CodeMeaning);
            elseif isfield(unitItem,'CodeValue')
                units = string(unitItem.CodeValue);
            end
        end
    end
end
