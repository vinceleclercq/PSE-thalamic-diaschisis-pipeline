%% PSE_FREESURFER_ROI_EXTRACT_02.m
% Corrected extraction of bilateral thalamic and hippocampal volumes
% from FreeSurfer aseg.stats for the 7 acute PSE subjects.
%
% Correction vs v01:
% - eTIV extraction now matches ONLY EstimatedTotalIntraCranialVol.
%   The v01 parser could accidentally read BrainSegVol-to-eTIV first.
%
% Outputs:
%   <PSE_ROOT>/derivatives/roi_stats/
%       PSE_FREESURFER_ROI_VOLUMES_v02.csv
%       PSE_FREESURFER_ROI_VOLUMES_v02.mat
%
% Definitions:
%   Ratio_LR = Left / Right
%   AI_LR    = (Left - Right) / (Left + Right)
%
% These are descriptive LEFT/RIGHT metrics only.
% Ipsi/contra metrics will be calculated later after adding stroke side.

clear; clc;

ROOT = string(getenv("PSE_ROOT"));
if strlength(ROOT)==0, ROOT = fullfile(string(getenv("HOME")),"PSE"); end
FS_DIR = fullfile(ROOT, "derivatives", "freesurfer");
OUT_DIR = fullfile(ROOT, "derivatives", "roi_stats");

if ~isfolder(OUT_DIR)
    mkdir(OUT_DIR);
end

subjects = "sub-P" + compose("%03d", 1:7) + "_acute";
n = numel(subjects);

Subject = strings(n,1);
Left_Thalamus_mm3 = nan(n,1);
Right_Thalamus_mm3 = nan(n,1);
Left_Hippocampus_mm3 = nan(n,1);
Right_Hippocampus_mm3 = nan(n,1);
eTIV_mm3 = nan(n,1);

for i = 1:n

    Subject(i) = subjects(i);
    statsFile = fullfile(FS_DIR, subjects(i), "stats", "aseg.stats");

    fprintf("\n[%d/%d] %s\n", i, n, subjects(i));

    if ~isfile(statsFile)
        warning("Missing aseg.stats: %s", statsFile);
        continue;
    end

    txt = fileread(statsFile);
    lines = splitlines(string(txt));

    Left_Thalamus_mm3(i) = findStructVolume(lines, ...
        ["Left-Thalamus-Proper","Left-Thalamus"]);

    Right_Thalamus_mm3(i) = findStructVolume(lines, ...
        ["Right-Thalamus-Proper","Right-Thalamus"]);

    Left_Hippocampus_mm3(i) = findStructVolume(lines, ...
        "Left-Hippocampus");

    Right_Hippocampus_mm3(i) = findStructVolume(lines, ...
        "Right-Hippocampus");

    eTIV_mm3(i) = findETIV(lines);

    fprintf("  L thalamus     : %.1f mm3\n", Left_Thalamus_mm3(i));
    fprintf("  R thalamus     : %.1f mm3\n", Right_Thalamus_mm3(i));
    fprintf("  L hippocampus  : %.1f mm3\n", Left_Hippocampus_mm3(i));
    fprintf("  R hippocampus  : %.1f mm3\n", Right_Hippocampus_mm3(i));
    fprintf("  eTIV           : %.1f mm3\n", eTIV_mm3(i));
end

% Descriptive left/right metrics
Thalamus_Ratio_LR = Left_Thalamus_mm3 ./ Right_Thalamus_mm3;
Thalamus_AI_LR = (Left_Thalamus_mm3 - Right_Thalamus_mm3) ./ ...
                 (Left_Thalamus_mm3 + Right_Thalamus_mm3);

Hippocampus_Ratio_LR = Left_Hippocampus_mm3 ./ Right_Hippocampus_mm3;
Hippocampus_AI_LR = (Left_Hippocampus_mm3 - Right_Hippocampus_mm3) ./ ...
                    (Left_Hippocampus_mm3 + Right_Hippocampus_mm3);

% Normalized volumes as percentage of estimated total intracranial volume
Left_Thalamus_pct_eTIV = 100 .* Left_Thalamus_mm3 ./ eTIV_mm3;
Right_Thalamus_pct_eTIV = 100 .* Right_Thalamus_mm3 ./ eTIV_mm3;
Left_Hippocampus_pct_eTIV = 100 .* Left_Hippocampus_mm3 ./ eTIV_mm3;
Right_Hippocampus_pct_eTIV = 100 .* Right_Hippocampus_mm3 ./ eTIV_mm3;

T = table( ...
    Subject, ...
    Left_Thalamus_mm3, Right_Thalamus_mm3, ...
    Thalamus_Ratio_LR, Thalamus_AI_LR, ...
    Left_Hippocampus_mm3, Right_Hippocampus_mm3, ...
    Hippocampus_Ratio_LR, Hippocampus_AI_LR, ...
    eTIV_mm3, ...
    Left_Thalamus_pct_eTIV, Right_Thalamus_pct_eTIV, ...
    Left_Hippocampus_pct_eTIV, Right_Hippocampus_pct_eTIV);

% Sanity checks
roiVars = ["Left_Thalamus_mm3","Right_Thalamus_mm3", ...
           "Left_Hippocampus_mm3","Right_Hippocampus_mm3"];

if any(any(ismissing(T(:, roiVars))))
    warning("At least one ROI value is missing. Review the table.");
end

if any(any(table2array(T(:, roiVars)) <= 0), "all")
    warning("At least one ROI volume is <= 0. Review the table.");
end

% Physiologic plausibility check for eTIV (broad safety range)
badETIV = isnan(eTIV_mm3) | eTIV_mm3 < 800000 | eTIV_mm3 > 2500000;
if any(badETIV)
    warning("One or more eTIV values are missing or implausible. Review aseg.stats.");
end

csvFile = fullfile(OUT_DIR, "PSE_FREESURFER_ROI_VOLUMES_v02.csv");
matFile = fullfile(OUT_DIR, "PSE_FREESURFER_ROI_VOLUMES_v02.mat");

writetable(T, csvFile);
save(matFile, "T");

fprintf("\n============================================================\n");
fprintf("ROI extraction complete.\n");
fprintf("CSV: %s\n", csvFile);
fprintf("MAT: %s\n", matFile);
fprintf("============================================================\n\n");

disp(T);

%% ------------------------------------------------------------------------
function vol = findStructVolume(lines, names)

    vol = NaN;
    names = string(names);

    for k = 1:numel(lines)
        line = strtrim(lines(k));

        if strlength(line) == 0 || startsWith(line, "#")
            continue;
        end

        tok = split(line);
        tok(tok == "") = [];

        if numel(tok) < 5
            continue;
        end

        structName = tok(5);

        if any(structName == names)
            x = str2double(tok(4));
            if ~isnan(x)
                vol = x;
                return;
            end
        end
    end
end

function etiv = findETIV(lines)
% Robustly extract the specific FreeSurfer measure:
%   EstimatedTotalIntraCranialVol
%
% Do NOT match generic "eTIV", because aseg.stats also contains measures
% such as BrainSegVol-to-eTIV.

    etiv = NaN;

    for k = 1:numel(lines)
        line = strtrim(lines(k));

        if contains(line, "EstimatedTotalIntraCranialVol")
            parts = split(line, ",");
            parts = strtrim(parts);

            % Standard format:
            % # Measure EstimatedTotalIntraCranialVol, eTIV,
            % Estimated Total Intracranial Volume, VALUE, mm^3
            if numel(parts) >= 5
                x = str2double(parts(end-1));
                if ~isnan(x)
                    etiv = x;
                    return;
                end
            end
        end
    end
end
