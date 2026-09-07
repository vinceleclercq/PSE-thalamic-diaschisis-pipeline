%% PSE_ACUTE_PILOT_ANALYSIS_01.m
% Acute pilot analysis of ipsilateral vs contralateral thalamic and
% hippocampal multimodal ROI measures in the PSE project.
%
% INPUT
% -----
% <PSE_ROOT>/derivatives/roi_values/final/
%   PSE_ACUTE_ROI_ANALYSIS_READY_LONG.csv
%
% EXPECTED VARIABLES
% ------------------
% Subject, Modality, Structure, StrokeSide,
% LeftMean, RightMean, Mean_Ratio_LR, Mean_AI_LR,
% LeftMedian, RightMedian, Median_Ratio_LR, Median_AI_LR,
% IpsiMean, ContraMean, Mean_Ratio_IpsiContra, Mean_AI_IpsiContra,
% IpsiMedian, ContraMedian, Median_Ratio_IpsiContra, Median_AI_IpsiContra,
% PrimaryPairUsable, PairQC, AnalysisFlag, ...
%
% PRIMARY ANALYSIS CHOICE
% -----------------------
% Median ROI values are used as the primary descriptive measure because
% they are more robust to residual outliers / partial-volume effects.
% Means are retained in the source dataset and can be switched on below.
%
% PRIMARY MODALITIES
% ------------------
% PET_FDG, ASL_CBF, CTP_CBF
%
% SECONDARY CTP PARAMETERS
% ------------------------
% CTP_CBV, CTP_MTT, CTP_TMAX, CTP_TTP
%
% IMPORTANT
% ---------
% This is a descriptive PILOT analysis (n=7). The script deliberately does
% not perform null-hypothesis significance testing. Pairwise modality
% concordance is summarized with Spearman rho and N only.
%
% AI definition:
%   AI_ipsi-contra = (Ipsi - Contra) / (Ipsi + Contra)
%
% Interpretation for PET_FDG / ASL_CBF / CTP_CBF:
%   AI < 0 : lower value ipsilateral to stroke
%   AI = 0 : symmetry
%   AI > 0 : higher value ipsilateral to stroke
%
% OUTPUT ROOT
% -----------
% <PSE_ROOT>/derivatives/pilot_analysis/
%
% CSV outputs:
%   PSE_ACUTE_PILOT_PRIMARY_LONG.csv
%   PSE_ACUTE_PILOT_PRIMARY_WIDE.csv
%   PSE_ACUTE_PILOT_DESCRIPTIVE_SUMMARY.csv
%   PSE_ACUTE_PILOT_MODALITY_CONCORDANCE.csv
%   PSE_ACUTE_PILOT_RANKING.csv
%   PSE_ACUTE_PILOT_SECONDARY_CTP.csv
%
% Figures:
%   heatmaps/
%   paired/
%   scatter/
%   ranking/
%
% Requires only base MATLAB.

clear; clc; close all;

%% ========================================================================
% Configuration

ROOT = string(getenv("PSE_ROOT"));
if strlength(ROOT)==0, ROOT = fullfile(string(getenv("HOME")),"PSE"); end

INFILE = fullfile(ROOT, "derivatives", "roi_values", "final", ...
    "PSE_ACUTE_ROI_ANALYSIS_READY_LONG.csv");

OUTROOT = fullfile(ROOT, "derivatives", "pilot_analysis");
CSVROOT = fullfile(OUTROOT, "csv");
FIGROOT = fullfile(OUTROOT, "figures");

HEATROOT = fullfile(FIGROOT, "heatmaps");
PAIREDROOT = fullfile(FIGROOT, "paired");
SCATTERROOT = fullfile(FIGROOT, "scatter");
RANKROOT = fullfile(FIGROOT, "ranking");

dirs = [OUTROOT, CSVROOT, FIGROOT, HEATROOT, PAIREDROOT, SCATTERROOT, RANKROOT];
for d = dirs
    if ~isfolder(d)
        mkdir(d);
    end
end

% Primary statistic for plots and summary:
PRIMARY_STAT = "Median";     % "Median" or "Mean"

PRIMARY_MODALITIES = ["PET_FDG", "ASL_CBF", "CTP_CBF"];
SECONDARY_CTP = ["CTP_CBV", "CTP_MTT", "CTP_TMAX", "CTP_TTP"];
STRUCTURES = ["Thalamus", "Hippocampus"];

% AI display range for heatmaps. Values outside are visually clipped only.
AI_DISPLAY_LIMIT = 0.25;

fprintf("============================================================\n");
fprintf("PSE ACUTE PILOT ANALYSIS\n");
fprintf("Primary statistic: %s\n", PRIMARY_STAT);
fprintf("============================================================\n");

%% ========================================================================
% Load and validate

if ~isfile(INFILE)
    error("Input file not found: %s", INFILE);
end

T = readtable(INFILE, "TextType", "string");

requiredVars = [ ...
    "Subject","Modality","Structure","StrokeSide", ...
    "IpsiMean","ContraMean","Mean_AI_IpsiContra", ...
    "IpsiMedian","ContraMedian","Median_AI_IpsiContra", ...
    "PrimaryPairUsable","PairQC"];

missingVars = setdiff(requiredVars, string(T.Properties.VariableNames));

if ~isempty(missingVars)
    error("Missing required variable(s): %s", strjoin(missingVars, ", "));
end

T.Subject = string(T.Subject);
T.Modality = string(T.Modality);
T.Structure = string(T.Structure);
T.StrokeSide = string(T.StrokeSide);
T.PairQC = string(T.PairQC);

% Normalize logical if CSV import returned strings / numeric.
if ~islogical(T.PrimaryPairUsable)
    if isnumeric(T.PrimaryPairUsable)
        T.PrimaryPairUsable = T.PrimaryPairUsable ~= 0;
    else
        x = lower(strtrim(string(T.PrimaryPairUsable)));
        T.PrimaryPairUsable = x=="true" | x=="1";
    end
end

%% ========================================================================
% Choose primary columns

if PRIMARY_STAT == "Median"
    ipsiVar = "IpsiMedian";
    contraVar = "ContraMedian";
    aiVar = "Median_AI_IpsiContra";
    ratioVar = "Median_Ratio_IpsiContra";
elseif PRIMARY_STAT == "Mean"
    ipsiVar = "IpsiMean";
    contraVar = "ContraMean";
    aiVar = "Mean_AI_IpsiContra";
    ratioVar = "Mean_Ratio_IpsiContra";
else
    error("PRIMARY_STAT must be Median or Mean.");
end

%% ========================================================================
% Primary analysis table

isPrimary = ismember(T.Modality, PRIMARY_MODALITIES) & ...
            ismember(T.Structure, STRUCTURES);

Tp = T(isPrimary, :);

% Keep explicit non-usable rows in the dataset but blank paired metrics.
bad = ~Tp.PrimaryPairUsable;

Tp.PrimaryIpsi = Tp.(ipsiVar);
Tp.PrimaryContra = Tp.(contraVar);
Tp.PrimaryAI = Tp.(aiVar);
Tp.PrimaryRatio = Tp.(ratioVar);

Tp.PrimaryIpsi(bad) = NaN;
Tp.PrimaryContra(bad) = NaN;
Tp.PrimaryAI(bad) = NaN;
Tp.PrimaryRatio(bad) = NaN;

primaryOut = fullfile(CSVROOT, "PSE_ACUTE_PILOT_PRIMARY_LONG.csv");
writetable(Tp, primaryOut);

%% ========================================================================
% Build compact wide table: one row per subject

subjects = unique(T.Subject, "stable");

Twide = table(subjects, 'VariableNames', {'Subject'});

strokeSide = repmat("", numel(subjects), 1);

for is = 1:numel(subjects)
    j = find(T.Subject == subjects(is) & T.StrokeSide ~= "", 1);
    if ~isempty(j)
        strokeSide(is) = T.StrokeSide(j);
    end
end

Twide.StrokeSide = strokeSide;

for i = 1:height(Tp)

    row = find(Twide.Subject == Tp.Subject(i), 1);

    prefix = matlab.lang.makeValidName( ...
        Tp.Modality(i) + "_" + Tp.Structure(i));

    variablesToPivot = [ ...
        "PrimaryIpsi", ...
        "PrimaryContra", ...
        "PrimaryAI", ...
        "PrimaryRatio"];

    for v = variablesToPivot

        newName = matlab.lang.makeValidName(prefix + "_" + v);

        if ~ismember(newName, string(Twide.Properties.VariableNames))
            Twide.(newName) = nan(height(Twide),1);
        end

        Twide.(newName)(row) = Tp.(v)(i);
    end
end

wideOut = fullfile(CSVROOT, "PSE_ACUTE_PILOT_PRIMARY_WIDE.csv");
writetable(Twide, wideOut);

%% ========================================================================
% Descriptive summary by modality and structure
%
% No hypothesis testing. Report:
% N, median AI, IQR, min, max, N(AI<0), proportion AI<0,
% median ipsi, median contra.

summaryRows = {};

for im = 1:numel(PRIMARY_MODALITIES)
    MOD = PRIMARY_MODALITIES(im);

    for istr = 1:numel(STRUCTURES)
        STR = STRUCTURES(istr);

        idx = Tp.Modality == MOD & Tp.Structure == STR & ...
              Tp.PrimaryPairUsable & isfinite(Tp.PrimaryAI);

        ai = Tp.PrimaryAI(idx);
        ipsi = Tp.PrimaryIpsi(idx);
        contra = Tp.PrimaryContra(idx);

        n = numel(ai);

        if n == 0
            medAI = NaN; q1 = NaN; q3 = NaN;
            minAI = NaN; maxAI = NaN;
            nNeg = 0; propNeg = NaN;
            medIpsi = NaN; medContra = NaN;
        else
            medAI = median(ai);
            q1 = localQuantile(ai, 0.25);
            q3 = localQuantile(ai, 0.75);
            minAI = min(ai);
            maxAI = max(ai);
            nNeg = nnz(ai < 0);
            propNeg = nNeg / n;
            medIpsi = median(ipsi, "omitnan");
            medContra = median(contra, "omitnan");
        end

        summaryRows(end+1,:) = { ... %#ok<AGROW>
            char(MOD), char(STR), n, ...
            medAI, q1, q3, q3-q1, minAI, maxAI, ...
            nNeg, propNeg, medIpsi, medContra};
    end
end

summaryVars = { ...
    'Modality','Structure','N', ...
    'MedianAI','Q1_AI','Q3_AI','IQR_AI','MinAI','MaxAI', ...
    'N_AI_Negative','Proportion_AI_Negative', ...
    'MedianIpsi','MedianContra'};

Tsummary = cell2table(summaryRows, 'VariableNames', summaryVars);

summaryOut = fullfile(CSVROOT, "PSE_ACUTE_PILOT_DESCRIPTIVE_SUMMARY.csv");
writetable(Tsummary, summaryOut);

%% ========================================================================
% Ranking by most negative AI (strongest ipsilateral reduction first)

rankingRows = {};

for im = 1:numel(PRIMARY_MODALITIES)
    MOD = PRIMARY_MODALITIES(im);

    for istr = 1:numel(STRUCTURES)
        STR = STRUCTURES(istr);

        idx = Tp.Modality == MOD & Tp.Structure == STR & ...
              Tp.PrimaryPairUsable & isfinite(Tp.PrimaryAI);

        Ts = Tp(idx, {'Subject','StrokeSide','PrimaryIpsi', ...
                      'PrimaryContra','PrimaryAI','PrimaryRatio'});

        if isempty(Ts)
            continue;
        end

        Ts = sortrows(Ts, "PrimaryAI", "ascend");

        for r = 1:height(Ts)
            rankingRows(end+1,:) = { ... %#ok<AGROW>
                char(MOD), char(STR), r, ...
                char(Ts.Subject(r)), char(Ts.StrokeSide(r)), ...
                Ts.PrimaryIpsi(r), Ts.PrimaryContra(r), ...
                Ts.PrimaryAI(r), Ts.PrimaryRatio(r)};
        end
    end
end

rankingVars = { ...
    'Modality','Structure','RankMostNegativeAI', ...
    'Subject','StrokeSide','Ipsi','Contra','AI','Ratio_IpsiContra'};

Tranking = cell2table(rankingRows, 'VariableNames', rankingVars);

rankingOut = fullfile(CSVROOT, "PSE_ACUTE_PILOT_RANKING.csv");
writetable(Tranking, rankingOut);

%% ========================================================================
% Modality concordance: pairwise Spearman rho on AI
%
% Pairs:
% PET_FDG vs ASL_CBF
% PET_FDG vs CTP_CBF
% ASL_CBF vs CTP_CBF
%
% Descriptive rho only; no p-values.

pairList = [ ...
    "PET_FDG", "ASL_CBF"; ...
    "PET_FDG", "CTP_CBF"; ...
    "ASL_CBF", "CTP_CBF"];

corrRows = {};

for istr = 1:numel(STRUCTURES)

    STR = STRUCTURES(istr);

    for ip = 1:size(pairList,1)

        M1 = pairList(ip,1);
        M2 = pairList(ip,2);

        [subCommon, x, y] = getPairedAI(Tp, STR, M1, M2);

        n = numel(x);

        if n >= 2
            rho = localSpearman(x, y);
        else
            rho = NaN;
        end

        corrRows(end+1,:) = { ... %#ok<AGROW>
            char(STR), char(M1), char(M2), n, rho, ...
            strjoin(subCommon, ";")};
    end
end

corrVars = { ...
    'Structure','Modality1','Modality2','N','SpearmanRho','Subjects'};

Tcorr = cell2table(corrRows, 'VariableNames', corrVars);

corrOut = fullfile(CSVROOT, "PSE_ACUTE_PILOT_MODALITY_CONCORDANCE.csv");
writetable(Tcorr, corrOut);

%% ========================================================================
% Secondary CTP table

isSecondary = ismember(T.Modality, SECONDARY_CTP) & ...
              ismember(T.Structure, STRUCTURES);

Tsec = T(isSecondary, :);

Tsec.PrimaryIpsi = Tsec.(ipsiVar);
Tsec.PrimaryContra = Tsec.(contraVar);
Tsec.PrimaryAI = Tsec.(aiVar);
Tsec.PrimaryRatio = Tsec.(ratioVar);

badSec = ~Tsec.PrimaryPairUsable;
Tsec.PrimaryIpsi(badSec) = NaN;
Tsec.PrimaryContra(badSec) = NaN;
Tsec.PrimaryAI(badSec) = NaN;
Tsec.PrimaryRatio(badSec) = NaN;

secOut = fullfile(CSVROOT, "PSE_ACUTE_PILOT_SECONDARY_CTP.csv");
writetable(Tsec, secOut);

%% ========================================================================
% FIGURE 1: AI heatmaps, one per structure

for istr = 1:numel(STRUCTURES)

    STR = STRUCTURES(istr);

    A = nan(numel(subjects), numel(PRIMARY_MODALITIES));

    for is = 1:numel(subjects)
        for im = 1:numel(PRIMARY_MODALITIES)

            idx = Tp.Subject == subjects(is) & ...
                  Tp.Structure == STR & ...
                  Tp.Modality == PRIMARY_MODALITIES(im) & ...
                  Tp.PrimaryPairUsable;

            if nnz(idx) == 1
                A(is,im) = Tp.PrimaryAI(idx);
            end
        end
    end

    fig = figure( ...
        "Visible","off", ...
        "Color","w", ...
        "Position",[100 100 900 650]);

    ax = axes(fig);

    % NaNs display as white.
    h = imagesc(ax, A);
    h.AlphaData = ~isnan(A);

    set(ax, "Color", [1 1 1]);
    caxis(ax, [-AI_DISPLAY_LIMIT AI_DISPLAY_LIMIT]);

    colormap(ax, localDivergingMap(256));
    cb = colorbar(ax);
    cb.Label.String = "AI = (Ipsi - Contra) / (Ipsi + Contra)";

    xticks(ax, 1:numel(PRIMARY_MODALITIES));
    xticklabels(ax, strrep(PRIMARY_MODALITIES, "_", " "));

    yticks(ax, 1:numel(subjects));
    yticklabels(ax, subjects);

    title(ax, sprintf("%s — acute ipsi/contra asymmetry", STR), ...
        "Interpreter","none");

    xlabel(ax, "Modality");
    ylabel(ax, "Subject");

    % Add numeric values.
    for r = 1:size(A,1)
        for c = 1:size(A,2)
            if isfinite(A(r,c))
                text(ax, c, r, sprintf("%.3f", A(r,c)), ...
                    "HorizontalAlignment","center", ...
                    "VerticalAlignment","middle", ...
                    "FontSize",9, ...
                    "FontWeight","bold");
            else
                text(ax, c, r, "NA", ...
                    "HorizontalAlignment","center", ...
                    "VerticalAlignment","middle", ...
                    "FontSize",9);
            end
        end
    end

    outPng = fullfile(HEATROOT, ...
        sprintf("PSE_%s_AI_heatmap.png", STR));

    exportgraphics(fig, outPng, "Resolution", 220);
    close(fig);
end

%% ========================================================================
% FIGURE 2: paired ipsi vs contra plots

for im = 1:numel(PRIMARY_MODALITIES)

    MOD = PRIMARY_MODALITIES(im);

    for istr = 1:numel(STRUCTURES)

        STR = STRUCTURES(istr);

        idx = Tp.Modality == MOD & ...
              Tp.Structure == STR & ...
              Tp.PrimaryPairUsable & ...
              isfinite(Tp.PrimaryIpsi) & ...
              isfinite(Tp.PrimaryContra);

        Ts = Tp(idx, :);

        if isempty(Ts)
            continue;
        end

        fig = figure( ...
            "Visible","off", ...
            "Color","w", ...
            "Position",[100 100 750 650]);

        ax = axes(fig);
        hold(ax, "on");

        for i = 1:height(Ts)
            plot(ax, [1 2], ...
                [Ts.PrimaryContra(i), Ts.PrimaryIpsi(i)], ...
                "-o", ...
                "LineWidth", 1.2, ...
                "MarkerSize", 6);

            text(ax, 2.04, Ts.PrimaryIpsi(i), Ts.Subject(i), ...
                "FontSize",8, ...
                "Interpreter","none");
        end

        xlim(ax, [0.7 2.45]);
        xticks(ax, [1 2]);
        xticklabels(ax, ["Contralateral","Ipsilateral"]);

        ylabel(ax, sprintf("%s ROI value", PRIMARY_STAT));
        title(ax, sprintf("%s — %s — acute paired values", MOD, STR), ...
            "Interpreter","none");

        grid(ax, "on");
        box(ax, "off");

        outPng = fullfile(PAIREDROOT, ...
            sprintf("PSE_%s_%s_paired_%s.png", MOD, STR, PRIMARY_STAT));

        exportgraphics(fig, outPng, "Resolution", 220);
        close(fig);
    end
end

%% ========================================================================
% FIGURE 3: pairwise modality AI scatterplots

for istr = 1:numel(STRUCTURES)

    STR = STRUCTURES(istr);

    for ip = 1:size(pairList,1)

        M1 = pairList(ip,1);
        M2 = pairList(ip,2);

        [subCommon, x, y] = getPairedAI(Tp, STR, M1, M2);

        if numel(x) < 2
            continue;
        end

        rho = localSpearman(x, y);

        fig = figure( ...
            "Visible","off", ...
            "Color","w", ...
            "Position",[100 100 700 650]);

        ax = axes(fig);
        hold(ax, "on");

        scatter(ax, x, y, 55, "filled");

        for i = 1:numel(x)
            text(ax, x(i), y(i), "  " + subCommon(i), ...
                "FontSize",8, ...
                "Interpreter","none");
        end

        % Zero reference axes.
        xline(ax, 0, "--");
        yline(ax, 0, "--");

        xlabel(ax, sprintf("%s AI", strrep(M1,"_"," ")));
        ylabel(ax, sprintf("%s AI", strrep(M2,"_"," ")));

        title(ax, sprintf("%s — modality concordance", STR), ...
            "Interpreter","none");

        subtitle(ax, sprintf("Spearman rho = %.3f | N = %d", rho, numel(x)));

        grid(ax, "on");
        box(ax, "off");
        axis(ax, "square");

        outPng = fullfile(SCATTERROOT, ...
            sprintf("PSE_%s_%s_vs_%s_AI_scatter.png", STR, M1, M2));

        exportgraphics(fig, outPng, "Resolution", 220);
        close(fig);
    end
end

%% ========================================================================
% FIGURE 4: ranking plots (most negative AI first)

for im = 1:numel(PRIMARY_MODALITIES)

    MOD = PRIMARY_MODALITIES(im);

    for istr = 1:numel(STRUCTURES)

        STR = STRUCTURES(istr);

        idx = Tranking.Modality == MOD & Tranking.Structure == STR;
        Ts = Tranking(idx, :);

        if isempty(Ts)
            continue;
        end

        fig = figure( ...
            "Visible","off", ...
            "Color","w", ...
            "Position",[100 100 850 600]);

        ax = axes(fig);

        barh(ax, 1:height(Ts), Ts.AI);

        yticklabels(ax, Ts.Subject);
        yticks(ax, 1:height(Ts));

        set(ax, "YDir", "reverse");

        xline(ax, 0, "--");

        xlabel(ax, "AI = (Ipsi - Contra) / (Ipsi + Contra)");
        ylabel(ax, "Subject");

        title(ax, sprintf("%s — %s — ranked asymmetry", MOD, STR), ...
            "Interpreter","none");

        grid(ax, "on");
        box(ax, "off");

        outPng = fullfile(RANKROOT, ...
            sprintf("PSE_%s_%s_AI_ranking.png", MOD, STR));

        exportgraphics(fig, outPng, "Resolution", 220);
        close(fig);
    end
end

%% ========================================================================
% Save MAT bundle

MATOUT = fullfile(OUTROOT, "PSE_ACUTE_PILOT_ANALYSIS_01.mat");

save(MATOUT, ...
    "T", "Tp", "Twide", "Tsummary", "Tcorr", ...
    "Tranking", "Tsec", ...
    "PRIMARY_STAT", "PRIMARY_MODALITIES", "SECONDARY_CTP", "STRUCTURES");

%% ========================================================================
% Console summary

fprintf("\n============================================================\n");
fprintf("PILOT ANALYSIS COMPLETE\n");
fprintf("============================================================\n");
fprintf("Primary long CSV:\n  %s\n", primaryOut);
fprintf("Primary wide CSV:\n  %s\n", wideOut);
fprintf("Descriptive summary:\n  %s\n", summaryOut);
fprintf("Modality concordance:\n  %s\n", corrOut);
fprintf("Ranking:\n  %s\n", rankingOut);
fprintf("Secondary CTP:\n  %s\n", secOut);
fprintf("Figures:\n  %s\n", FIGROOT);
fprintf("MAT bundle:\n  %s\n", MATOUT);

fprintf("\nDescriptive AI summary:\n");
disp(Tsummary);

fprintf("\nPairwise modality concordance:\n");
disp(Tcorr);

fprintf("============================================================\n");

%% ========================================================================
% Local functions
% ========================================================================

function q = localQuantile(x, p)
% Simple linear-interpolated quantile, base MATLAB only.

    x = sort(x(isfinite(x)));

    n = numel(x);

    if n == 0
        q = NaN;
        return;
    elseif n == 1
        q = x;
        return;
    end

    pos = 1 + (n-1)*p;
    lo = floor(pos);
    hi = ceil(pos);

    if lo == hi
        q = x(lo);
    else
        w = pos - lo;
        q = (1-w)*x(lo) + w*x(hi);
    end
end

function rho = localSpearman(x, y)
% Spearman rank correlation without Statistics Toolbox.

    x = x(:);
    y = y(:);

    keep = isfinite(x) & isfinite(y);
    x = x(keep);
    y = y(keep);

    if numel(x) < 2
        rho = NaN;
        return;
    end

    rx = localTiedRank(x);
    ry = localTiedRank(y);

    C = corrcoef(rx, ry);

    if numel(C) >= 4
        rho = C(1,2);
    else
        rho = NaN;
    end
end

function r = localTiedRank(x)
% Average ranks for ties, base MATLAB only.

    x = x(:);
    n = numel(x);

    [sx, order] = sort(x);
    rSorted = nan(n,1);

    i = 1;

    while i <= n
        j = i;

        while j < n && sx(j+1) == sx(i)
            j = j + 1;
        end

        avgRank = mean(i:j);
        rSorted(i:j) = avgRank;

        i = j + 1;
    end

    r = nan(n,1);
    r(order) = rSorted;
end

function [subjectsCommon, x, y] = getPairedAI(Tp, structure, modality1, modality2)

    T1 = Tp( ...
        Tp.Structure == structure & ...
        Tp.Modality == modality1 & ...
        Tp.PrimaryPairUsable & ...
        isfinite(Tp.PrimaryAI), ...
        {'Subject','PrimaryAI'});

    T2 = Tp( ...
        Tp.Structure == structure & ...
        Tp.Modality == modality2 & ...
        Tp.PrimaryPairUsable & ...
        isfinite(Tp.PrimaryAI), ...
        {'Subject','PrimaryAI'});

    subjectsCommon = intersect(T1.Subject, T2.Subject, "stable");

    x = nan(numel(subjectsCommon),1);
    y = nan(numel(subjectsCommon),1);

    for i = 1:numel(subjectsCommon)

        j1 = find(T1.Subject == subjectsCommon(i), 1);
        j2 = find(T2.Subject == subjectsCommon(i), 1);

        x(i) = T1.PrimaryAI(j1);
        y(i) = T2.PrimaryAI(j2);
    end
end

function cmap = localDivergingMap(n)
% Blue -> white -> red diverging colormap centered on zero.

    if nargin < 1
        n = 256;
    end

    n1 = floor(n/2);
    n2 = n - n1;

    blue = [0.15 0.35 0.80];
    white = [1.00 1.00 1.00];
    red = [0.80 0.20 0.20];

    a = linspace(0,1,n1)';
    b = linspace(0,1,n2)';

    left = blue.*(1-a) + white.*a;
    right = white.*(1-b) + red.*b;

    cmap = [left; right];
end
