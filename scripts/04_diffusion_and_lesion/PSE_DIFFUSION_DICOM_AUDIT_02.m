%% PSE_DIFFUSION_DICOM_AUDIT_02.m
% Audit ORIGINAL acute diffusion DICOM metadata to determine whether an ADC
% series exists and to identify the b-values / meaning of DWI_001 vs DWI_002.
%
% Default subjects:
%   sub-P003, sub-P004, sub-P005
%
% The script reads source DICOM headers only. It does NOT modify DICOMs.
% v02 stores heterogeneous dicominfo outputs in cells because different
% DICOM series may expose different metadata fields.
%
% It reports one row per detected DICOM series with:
%   Subject
%   SeriesInstanceUID
%   SeriesNumber
%   SeriesDescription
%   ProtocolName
%   SequenceName
%   ImageType
%   Manufacturer
%   Modality
%   NFiles
%   Rows / Columns
%   SliceThickness
%   SpacingBetweenSlices
%   DiffusionBValue_Min / Max / Values
%   DiffusionDirectionality
%   IsADC_Like
%   IsDWI_Like
%
% It also prints likely ADC / DWI series to the MATLAB console.
%
% SOURCE DIRECTORY EXPECTED
% -------------------------
% <PSE_ROOT>/sub-P00X/acute/
%
% If your source DICOMs are elsewhere under the subject folder, the script
% searches recursively below <PSE_ROOT>/sub-P00X/acute.
%
% OUTPUT
% ------
% <PSE_ROOT>/derivatives/coreg_diffusion/
%   PSE_DIFFUSION_DICOM_AUDIT.csv

clear; clc;

ROOT = string(getenv("PSE_ROOT"));
if strlength(ROOT)==0, ROOT = fullfile(string(getenv("HOME")),"PSE"); end
OUT_ROOT = fullfile(ROOT,"derivatives","coreg_diffusion");

subjects = ["sub-P003","sub-P004","sub-P005"];

if ~isfolder(OUT_ROOT), mkdir(OUT_ROOT); end

fprintf("============================================================\n");
fprintf("PSE ORIGINAL DICOM DIFFUSION AUDIT\n");
fprintf("============================================================\n");

allRows = {};

for is = 1:numel(subjects)

    SUB = subjects(is);
    sourceRoot = fullfile(ROOT,SUB,"acute");

    fprintf("\n[%d/%d] %s\n",is,numel(subjects),SUB);

    if ~isfolder(sourceRoot)
        fprintf("  Source directory not found: %s\n",sourceRoot);
        continue;
    end

    d = dir(fullfile(sourceRoot,"**","*"));
    d = d(~[d.isdir]);

    % Exclude obvious non-source / hidden files.
    keep = true(numel(d),1);
    for k=1:numel(d)
        nm = string(d(k).name);
        fp = string(fullfile(d(k).folder,d(k).name));
        if startsWith(nm,".") || contains(fp,filesep+"derivatives"+filesep)
            keep(k)=false;
        end
    end
    d=d(keep);

    fprintf("  Candidate files to inspect: %d\n",numel(d));

    % Collect per-file headers.
    % dicominfo structures are not guaranteed to have identical fields
    % across series/vendors, so store them in a cell array rather than a
    % MATLAB struct array.
    headers = cell(0,1);
    filePaths = strings(0,1);

    for k=1:numel(d)

        fp = string(fullfile(d(k).folder,d(k).name));

        try
            info = dicominfo(fp, "UseDictionaryVR", true);
        catch
            continue;
        end

        headers{end+1,1} = info; %#ok<AGROW>
        filePaths(end+1,1) = fp; %#ok<AGROW>
    end

    if isempty(headers)
        fprintf("  No readable DICOM files found.\n");
        continue;
    end

    fprintf("  Readable DICOM files: %d\n",numel(headers));

    % Group by SeriesInstanceUID, falling back to SeriesNumber + description.
    keys = strings(numel(headers),1);

    for k=1:numel(headers)
        Hk = headers{k};
        uid = getFieldText(Hk,"SeriesInstanceUID","");
        if strlength(uid)>0
            keys(k)=uid;
        else
            sn = getFieldNum(Hk,"SeriesNumber",NaN);
            sd = getFieldText(Hk,"SeriesDescription","");
            keys(k)="SERIES_"+string(sn)+"_"+sd;
        end
    end

    [uKeys,~,g] = unique(keys,"stable");

    for ig=1:numel(uKeys)

        idx = find(g==ig);
        H = headers(idx);
        H0 = H{1};

        seriesUID = getFieldText(H0,"SeriesInstanceUID","");
        seriesNumber = getFieldNum(H0,"SeriesNumber",NaN);
        seriesDesc = getFieldText(H0,"SeriesDescription","");
        protocolName = getFieldText(H0,"ProtocolName","");
        sequenceName = getFieldText(H0,"SequenceName","");
        manufacturer = getFieldText(H0,"Manufacturer","");
        modality = getFieldText(H0,"Modality","");
        diffusionDirectionality = getFieldText(H0,"DiffusionDirectionality","");

        imageTypeVals = strings(numel(H),1);
        bVals = nan(numel(H),1);
        rowsVals = nan(numel(H),1);
        colVals = nan(numel(H),1);
        sliceThickVals = nan(numel(H),1);
        spacingVals = nan(numel(H),1);

        for j=1:numel(H)
            Hj = H{j};
            imageTypeVals(j) = getImageTypeText(Hj);
            bVals(j) = getDiffusionBValue(Hj);
            rowsVals(j) = getFieldNum(Hj,"Rows",NaN);
            colVals(j) = getFieldNum(Hj,"Columns",NaN);
            sliceThickVals(j) = getFieldNum(Hj,"SliceThickness",NaN);
            spacingVals(j) = getFieldNum(Hj,"SpacingBetweenSlices",NaN);
        end

        imageTypeUnique = unique(imageTypeVals(imageTypeVals~=""),"stable");
        imageTypeText = strjoin(imageTypeUnique," || ");

        finiteB = bVals(isfinite(bVals));
        if isempty(finiteB)
            bMin=NaN; bMax=NaN; bText="";
        else
            ub = unique(round(finiteB,6),"stable");
            bMin=min(ub);
            bMax=max(ub);
            bText=strjoin(string(ub),";");
        end

        rowsMode = localModeFinite(rowsVals);
        colsMode = localModeFinite(colVals);
        thickMed = median(sliceThickVals(isfinite(sliceThickVals)),"omitnan");
        spacingMed = median(spacingVals(isfinite(spacingVals)),"omitnan");

        combined = lower(strjoin([seriesDesc,protocolName,sequenceName,imageTypeText]," "));

        isADC = contains(combined,"adc") || ...
                contains(combined,"apparent diffusion") || ...
                contains(combined,"apparent_diffusion") || ...
                contains(combined,"derived") && contains(combined,"diffusion") && contains(combined,"map");

        isDWI = contains(combined,"dwi") || ...
                contains(combined,"diff") || ...
                ~isempty(finiteB);

        allRows(end+1,:) = { ... %#ok<AGROW>
            char(SUB), ...
            char(seriesUID), ...
            seriesNumber, ...
            char(seriesDesc), ...
            char(protocolName), ...
            char(sequenceName), ...
            char(imageTypeText), ...
            char(manufacturer), ...
            char(modality), ...
            numel(idx), ...
            rowsMode, ...
            colsMode, ...
            thickMed, ...
            spacingMed, ...
            bMin, ...
            bMax, ...
            char(bText), ...
            char(diffusionDirectionality), ...
            isADC, ...
            isDWI, ...
            char(filePaths(idx(1)))};
    end
end

varNames = { ...
    'Subject', ...
    'SeriesInstanceUID', ...
    'SeriesNumber', ...
    'SeriesDescription', ...
    'ProtocolName', ...
    'SequenceName', ...
    'ImageType', ...
    'Manufacturer', ...
    'Modality', ...
    'NFiles', ...
    'Rows', ...
    'Columns', ...
    'SliceThickness', ...
    'SpacingBetweenSlices', ...
    'DiffusionBValue_Min', ...
    'DiffusionBValue_Max', ...
    'DiffusionBValue_Values', ...
    'DiffusionDirectionality', ...
    'IsADC_Like', ...
    'IsDWI_Like', ...
    'ExampleFile'};

if isempty(allRows)
    T = cell2table(cell(0,numel(varNames)),'VariableNames',varNames);
else
    T = cell2table(allRows,'VariableNames',varNames);
end

% Sort likely ADC/DWI first.
if ~isempty(T)
    priority = 9*ones(height(T),1);
    priority(T.IsADC_Like) = 1;
    priority(~T.IsADC_Like & T.IsDWI_Like) = 2;
    T.Priority = priority;
    T = sortrows(T,{'Subject','Priority','SeriesNumber'});
    T.Priority = [];
end

outCSV = fullfile(OUT_ROOT,"PSE_DIFFUSION_DICOM_AUDIT.csv");
writetable(T,outCSV);

fprintf("\n============================================================\n");
fprintf("DICOM DIFFUSION AUDIT COMPLETE\n");
fprintf("CSV: %s\n",outCSV);
fprintf("============================================================\n");

for SUB = subjects
    fprintf("\n--- %s : likely diffusion series ---\n",SUB);
    idx = T.Subject==SUB & (T.IsADC_Like | T.IsDWI_Like);
    disp(T(idx, {'Subject','SeriesNumber','SeriesDescription','ProtocolName', ...
        'ImageType','NFiles','Rows','Columns','SliceThickness', ...
        'DiffusionBValue_Min','DiffusionBValue_Max','DiffusionBValue_Values', ...
        'IsADC_Like','IsDWI_Like'}));
end

%% ========================================================================
% Local functions

function s = getFieldText(info,field,defaultValue)
    if isfield(info,field)
        v = info.(field);
        try
            if iscell(v)
                s = strjoin(string(v),"\");
            elseif ischar(v) || isstring(v)
                s = string(v);
            elseif isnumeric(v)
                s = string(v);
            else
                s = string(v);
            end
        catch
            s = string(defaultValue);
        end
    else
        s = string(defaultValue);
    end
end

function x = getFieldNum(info,field,defaultValue)
    if isfield(info,field)
        v=info.(field);
        if isnumeric(v) && ~isempty(v)
            x=double(v(1));
        else
            x=str2double(string(v));
            if ~isfinite(x), x=defaultValue; end
        end
    else
        x=defaultValue;
    end
end

function s = getImageTypeText(info)
    if ~isfield(info,"ImageType")
        s="";
        return;
    end
    v=info.ImageType;
    if iscell(v)
        s=strjoin(string(v),"\");
    elseif ischar(v)
        s=string(v);
    else
        s=string(v);
    end
end

function b = getDiffusionBValue(info)
    b=NaN;

    % Standard MATLAB-parsed field.
    candidateFields = ["DiffusionBValue","BValue","Diffusion_b_value"];

    for f = candidateFields
        if isfield(info,f)
            v=info.(f);
            if isnumeric(v) && ~isempty(v)
                b=double(v(1));
                return;
            else
                z=str2double(string(v));
                if isfinite(z)
                    b=z;
                    return;
                end
            end
        end
    end

    % Some vendors put b-value in nested diffusion sequence items.
    fns=fieldnames(info);
    for i=1:numel(fns)
        fn=fns{i};
        if contains(lower(fn),"diffusion") && isstruct(info.(fn))
            b=searchNestedB(info.(fn));
            if isfinite(b), return; end
        end
    end
end

function b = searchNestedB(S)
    b=NaN;

    if ~isstruct(S), return; end

    if numel(S)>1
        for k=1:numel(S)
            b=searchNestedB(S(k));
            if isfinite(b), return; end
        end
        return;
    end

    fields=fieldnames(S);

    for i=1:numel(fields)
        fn=fields{i};
        v=S.(fn);

        if contains(lower(fn),"bvalue") || contains(lower(fn),"diffusionb")
            if isnumeric(v) && ~isempty(v)
                b=double(v(1));
                return;
            else
                z=str2double(string(v));
                if isfinite(z)
                    b=z;
                    return;
                end
            end
        end

        if isstruct(v)
            b=searchNestedB(v);
            if isfinite(b), return; end
        end
    end
end

function m = localModeFinite(x)
    x=x(isfinite(x));
    if isempty(x)
        m=NaN;
    else
        m=mode(x);
    end
end
