%% PSE_BUILD_ANALYSIS_READY_ROI_DATASET_01.m
% Build a clean, analysis-ready acute multimodal ROI dataset from the
% validated v03 extraction tables.
%
% This script DOES NOT resample images and DOES NOT re-extract voxel values.
% It starts from the already generated:
%
%   <PSE_ROOT>/derivatives/roi_values/
%       PSE_ACUTE_ROI_MODALITY_VALUES_LONG_v03.csv
%       PSE_ACUTE_ROI_MODALITY_SUMMARY_LR_v03.csv
%
% and, if available:
%
%   <PSE_ROOT>/derivatives/qc_projected_rois/
%       PSE_PROJECTED_ROI_QC_STATUS.csv
%
% Outputs
% -------
%   <PSE_ROOT>/derivatives/roi_values/final/
%
%   1) PSE_ACUTE_ROI_ANALYSIS_READY_LONG.csv
%      One row per Subject x Modality x Structure.
%
%   2) PSE_ACUTE_ROI_ANALYSIS_READY_WIDE.csv
%      One row per subject, with analysis-ready variables as columns.
%
%   3) PSE_ACUTE_ROI_EXCLUSIONS_FLAGS.csv
%      Missing modalities, insufficient CTP coverage, etc.
%
%   4) PSE_STROKE_SIDE_TEMPLATE.csv
%      Fill StrokeSide with L or R, then rerun this script to automatically
%      calculate ipsilateral / contralateral values and asymmetry indices.
%
%   5) PSE_ACUTE_ROI_ANALYSIS_READY.mat
%
% Definitions
% -----------
% Left/right:
%   Ratio_LR = Left / Right
%   AI_LR    = (Left - Right) / (Left + Right)
%
% If StrokeSide is available:
%   Ratio_IpsiContra = Ipsi / Contra
%   AI_IpsiContra    = (Ipsi - Contra) / (Ipsi + Contra)
%
% IMPORTANT:
% - CTP pairs with PairQC == FAIL_COVERAGE_LT50 are excluded from primary
%   paired analysis. Their ipsi/contra and LR paired metrics are left NaN.
% - WARN_COVERAGE_LT80 is retained but flagged.
% - Missing modalities remain explicit rather than being silently dropped.

clear; clc;

%% ========================================================================
% Paths

ROOT = string(getenv("PSE_ROOT"));
if strlength(ROOT)==0, ROOT = fullfile(string(getenv("HOME")),"PSE"); end

ROI_VALUE_DIR = fullfile(ROOT, "derivatives", "roi_values");
QC_DIR = fullfile(ROOT, "derivatives", "qc_projected_rois");
OUT_DIR = fullfile(ROI_VALUE_DIR, "final");

if ~isfolder(OUT_DIR)
    mkdir(OUT_DIR);
end

LONG_IN = fullfile(ROI_VALUE_DIR, ...
    "PSE_ACUTE_ROI_MODALITY_VALUES_LONG_v03.csv");

SUMMARY_IN = fullfile(ROI_VALUE_DIR, ...
    "PSE_ACUTE_ROI_MODALITY_SUMMARY_LR_v03.csv");

QC_IN = fullfile(QC_DIR, ...
    "PSE_PROJECTED_ROI_QC_STATUS.csv");

STROKE_SIDE_FILE = fullfile(OUT_DIR, ...
    "PSE_STROKE_SIDE_TEMPLATE.csv");

LONG_OUT = fullfile(OUT_DIR, ...
    "PSE_ACUTE_ROI_ANALYSIS_READY_LONG.csv");

WIDE_OUT = fullfile(OUT_DIR, ...
    "PSE_ACUTE_ROI_ANALYSIS_READY_WIDE.csv");

FLAGS_OUT = fullfile(OUT_DIR, ...
    "PSE_ACUTE_ROI_EXCLUSIONS_FLAGS.csv");

MAT_OUT = fullfile(OUT_DIR, ...
    "PSE_ACUTE_ROI_ANALYSIS_READY.mat");

%% ========================================================================
% Check inputs

if ~isfile(LONG_IN)
    error("Missing input file: %s", LONG_IN);
end

if ~isfile(SUMMARY_IN)
    error("Missing input file: %s", SUMMARY_IN);
end

Tvoxel = readtable(LONG_IN, "TextType", "string");
T = readtable(SUMMARY_IN, "TextType", "string");

fprintf("============================================================\n");
fprintf("PSE ANALYSIS-READY ROI DATASET\n");
fprintf("============================================================\n");
fprintf("Rows in long voxel summary: %d\n", height(Tvoxel));
fprintf("Rows in paired LR summary:  %d\n", height(T));

%% ========================================================================
% Normalize key columns to strings

T.Subject = string(T.Subject);
T.Modality = string(T.Modality);
T.Structure = string(T.Structure);
T.PairQC = string(T.PairQC);

%% ========================================================================
% Create / read stroke-side file

subjects = unique(T.Subject, "stable");

if ~isfile(STROKE_SIDE_FILE)

    Tside = table(subjects, repmat("", numel(subjects), 1), ...
        'VariableNames', {'Subject','StrokeSide'});

    writetable(Tside, STROKE_SIDE_FILE);

    fprintf("\nCreated stroke-side template:\n  %s\n", STROKE_SIDE_FILE);
    fprintf("Fill StrokeSide with L or R and rerun to calculate ipsi/contra metrics.\n");

else

    Tside = readtable(STROKE_SIDE_FILE, "TextType", "string");

    if ~all(ismember(["Subject","StrokeSide"], string(Tside.Properties.VariableNames)))
        error("Stroke-side file must contain columns Subject and StrokeSide.");
    end

    Tside.Subject = string(Tside.Subject);
    Tside.StrokeSide = upper(strtrim(string(Tside.StrokeSide)));
end

%% ========================================================================
% Merge StrokeSide onto paired summary

T.StrokeSide = repmat("", height(T), 1);

for i = 1:height(T)
    j = find(Tside.Subject == T.Subject(i), 1);

    if ~isempty(j)
        side = Tside.StrokeSide(j);

        if side == "L" || side == "R"
            T.StrokeSide(i) = side;
        end
    end
end

%% ========================================================================
% Primary usability flags

T.PrimaryPairUsable = true(height(T),1);
T.AnalysisFlag = repmat("OK", height(T),1);

for i = 1:height(T)

    qc = T.PairQC(i);

    if qc == "FAIL_COVERAGE_LT50"
        T.PrimaryPairUsable(i) = false;
        T.AnalysisFlag(i) = "EXCLUDE_PRIMARY_PAIR_LOW_COVERAGE";

    elseif qc == "WARN_COVERAGE_LT80"
        T.AnalysisFlag(i) = "RETAIN_WITH_COVERAGE_WARNING";

    elseif qc == "OK"
        T.AnalysisFlag(i) = "OK";

    else
        T.AnalysisFlag(i) = "REVIEW_QC";
    end
end

%% ========================================================================
% Ipsi / contra variables
%
% These remain NaN until StrokeSide has been entered as L/R.

n = height(T);

IpsiMean = nan(n,1);
ContraMean = nan(n,1);
Mean_Ratio_IpsiContra = nan(n,1);
Mean_AI_IpsiContra = nan(n,1);

IpsiMedian = nan(n,1);
ContraMedian = nan(n,1);
Median_Ratio_IpsiContra = nan(n,1);
Median_AI_IpsiContra = nan(n,1);

Ipsi_StatsCoverageFraction = nan(n,1);
Contra_StatsCoverageFraction = nan(n,1);

for i = 1:n

    if ~T.PrimaryPairUsable(i)
        continue;
    end

    side = T.StrokeSide(i);

    if side == "L"

        IpsiMean(i) = T.LeftMean(i);
        ContraMean(i) = T.RightMean(i);

        IpsiMedian(i) = T.LeftMedian(i);
        ContraMedian(i) = T.RightMedian(i);

        Ipsi_StatsCoverageFraction(i) = T.Left_StatsCoverageFraction(i);
        Contra_StatsCoverageFraction(i) = T.Right_StatsCoverageFraction(i);

    elseif side == "R"

        IpsiMean(i) = T.RightMean(i);
        ContraMean(i) = T.LeftMean(i);

        IpsiMedian(i) = T.RightMedian(i);
        ContraMedian(i) = T.LeftMedian(i);

        Ipsi_StatsCoverageFraction(i) = T.Right_StatsCoverageFraction(i);
        Contra_StatsCoverageFraction(i) = T.Left_StatsCoverageFraction(i);

    else
        continue;
    end

    Mean_Ratio_IpsiContra(i) = safeDivide(IpsiMean(i), ContraMean(i));
    Mean_AI_IpsiContra(i) = safeAI(IpsiMean(i), ContraMean(i));

    Median_Ratio_IpsiContra(i) = safeDivide(IpsiMedian(i), ContraMedian(i));
    Median_AI_IpsiContra(i) = safeAI(IpsiMedian(i), ContraMedian(i));
end

T.IpsiMean = IpsiMean;
T.ContraMean = ContraMean;
T.Mean_Ratio_IpsiContra = Mean_Ratio_IpsiContra;
T.Mean_AI_IpsiContra = Mean_AI_IpsiContra;

T.IpsiMedian = IpsiMedian;
T.ContraMedian = ContraMedian;
T.Median_Ratio_IpsiContra = Median_Ratio_IpsiContra;
T.Median_AI_IpsiContra = Median_AI_IpsiContra;

T.Ipsi_StatsCoverageFraction = Ipsi_StatsCoverageFraction;
T.Contra_StatsCoverageFraction = Contra_StatsCoverageFraction;

%% ========================================================================
% Add projected-ROI visual QC status, if available

T.ProjectedROIVisualQC = repmat("NOT_AVAILABLE", height(T),1);

if isfile(QC_IN)

    Tqc = readtable(QC_IN, "TextType", "string");

    Tqc.Subject = string(Tqc.Subject);
    Tqc.Modality = string(Tqc.Modality);
    Tqc.Status = string(Tqc.Status);

    for i = 1:height(T)

        j = find( ...
            Tqc.Subject == T.Subject(i) & ...
            Tqc.Modality == T.Modality(i), ...
            1);

        if ~isempty(j)
            T.ProjectedROIVisualQC(i) = Tqc.Status(j);
        end
    end
end

%% ========================================================================
% Ensure failed paired CTP comparisons cannot accidentally be analyzed

bad = ~T.PrimaryPairUsable;

pairedVarsToBlank = [ ...
    "Mean_Ratio_LR", ...
    "Mean_AI_LR", ...
    "Median_Ratio_LR", ...
    "Median_AI_LR", ...
    "IpsiMean", ...
    "ContraMean", ...
    "Mean_Ratio_IpsiContra", ...
    "Mean_AI_IpsiContra", ...
    "IpsiMedian", ...
    "ContraMedian", ...
    "Median_Ratio_IpsiContra", ...
    "Median_AI_IpsiContra"];

for v = pairedVarsToBlank
    if ismember(v, string(T.Properties.VariableNames))
        T.(v)(bad) = NaN;
    end
end

%% ========================================================================
% Create explicit exclusions / warnings table

flagRows = {};

% Pair-level flags from quantitative coverage
for i = 1:height(T)

    if T.AnalysisFlag(i) ~= "OK"
        flagRows(end+1,:) = { ... %#ok<AGROW>
            char(T.Subject(i)), ...
            char(T.Modality(i)), ...
            char(T.Structure(i)), ...
            char(T.AnalysisFlag(i)), ...
            char(T.PairQC(i))};
    end
end

% Missing modalities from visual-QC status table
if exist("Tqc", "var")

    miss = startsWith(Tqc.Status, "MISSING");

    for i = find(miss)'

        flagRows(end+1,:) = { ... %#ok<AGROW>
            char(Tqc.Subject(i)), ...
            char(Tqc.Modality(i)), ...
            "ALL", ...
            "MISSING_MODALITY", ...
            char(Tqc.Message(i))};
    end
end

flagVars = {'Subject','Modality','Structure','Flag','Detail'};

if isempty(flagRows)
    Tflags = cell2table(cell(0,5), "VariableNames", flagVars);
else
    Tflags = cell2table(flagRows, "VariableNames", flagVars);
end

%% ========================================================================
% Save clean LONG table

writetable(T, LONG_OUT);

%% ========================================================================
% Build WIDE table: one row per subject
%
% Only analysis-relevant fields are pivoted.
% Column example:
%   PET_FDG_Thalamus_LeftMean
%   ASL_CBF_Hippocampus_Mean_AI_IpsiContra
%   CTP_CBF_Thalamus_Mean_AI_LR

wideVars = [ ...
    "LeftMean", ...
    "RightMean", ...
    "Mean_Ratio_LR", ...
    "Mean_AI_LR", ...
    "LeftMedian", ...
    "RightMedian", ...
    "Median_Ratio_LR", ...
    "Median_AI_LR", ...
    "Left_StatsCoverageFraction", ...
    "Right_StatsCoverageFraction", ...
    "IpsiMean", ...
    "ContraMean", ...
    "Mean_Ratio_IpsiContra", ...
    "Mean_AI_IpsiContra", ...
    "IpsiMedian", ...
    "ContraMedian", ...
    "Median_Ratio_IpsiContra", ...
    "Median_AI_IpsiContra", ...
    "PrimaryPairUsable"];

Twide = table(subjects, 'VariableNames', {'Subject'});

% Add subject-level stroke side once
strokeSideWide = repmat("", numel(subjects),1);

for is = 1:numel(subjects)
    j = find(Tside.Subject == subjects(is), 1);

    if ~isempty(j)
        strokeSideWide(is) = Tside.StrokeSide(j);
    end
end

Twide.StrokeSide = strokeSideWide;

for i = 1:height(T)

    sub = T.Subject(i);

    rowIdx = find(Twide.Subject == sub, 1);

    prefix = matlab.lang.makeValidName( ...
        T.Modality(i) + "_" + T.Structure(i));

    for v = wideVars

        if ~ismember(v, string(T.Properties.VariableNames))
            continue;
        end

        newVar = matlab.lang.makeValidName(prefix + "_" + v);

        if ~ismember(newVar, string(Twide.Properties.VariableNames))

            if islogical(T.(v))
                Twide.(newVar) = false(height(Twide),1);
            else
                Twide.(newVar) = nan(height(Twide),1);
            end
        end

        Twide.(newVar)(rowIdx) = T.(v)(i);
    end
end

writetable(Twide, WIDE_OUT);
writetable(Tflags, FLAGS_OUT);

save(MAT_OUT, "T", "Twide", "Tflags", "Tvoxel", "Tside");

%% ========================================================================
% Summary

fprintf("\n============================================================\n");
fprintf("ANALYSIS-READY DATASET CREATED\n");
fprintf("============================================================\n");
fprintf("Long table:\n  %s\n", LONG_OUT);
fprintf("Wide table:\n  %s\n", WIDE_OUT);
fprintf("Flags/exclusions:\n  %s\n", FLAGS_OUT);
fprintf("Stroke-side template:\n  %s\n", STROKE_SIDE_FILE);
fprintf("MAT file:\n  %s\n", MAT_OUT);

fprintf("\nPair QC summary:\n");
disp(groupsummary(T, "AnalysisFlag"));

fprintf("\nStroke side currently available for %d / %d subjects.\n", ...
    nnz(Twide.StrokeSide=="L" | Twide.StrokeSide=="R"), height(Twide));

if all(Twide.StrokeSide=="")
    fprintf("=> Fill PSE_STROKE_SIDE_TEMPLATE.csv, then rerun this script.\n");
end

fprintf("============================================================\n");

%% ========================================================================
% Local functions

function y = safeDivide(a,b)

    if isempty(a) || isempty(b) || ...
       ~isfinite(a) || ~isfinite(b) || b == 0
        y = NaN;
    else
        y = a / b;
    end
end

function y = safeAI(ipsi,contra)

    d = ipsi + contra;

    if isempty(ipsi) || isempty(contra) || ...
       ~isfinite(ipsi) || ~isfinite(contra) || d == 0
        y = NaN;
    else
        y = (ipsi - contra) / d;
    end
end
