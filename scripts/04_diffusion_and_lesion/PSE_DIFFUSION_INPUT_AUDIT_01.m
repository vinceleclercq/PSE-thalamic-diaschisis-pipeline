%% PSE_DIFFUSION_INPUT_AUDIT_01.m
% Audit acute diffusion-related NIfTI files.
% No files are modified.

clear; clc;

ROOT = string(getenv("PSE_ROOT"));
if strlength(ROOT)==0, ROOT = fullfile(string(getenv("HOME")),"PSE"); end
NIFTI_ROOT = fullfile(ROOT, "derivatives", "nifti");
OUT_ROOT = fullfile(ROOT, "derivatives", "coreg_diffusion");

SPM_DIR = string(getenv("SPM_DIR"));
if strlength(SPM_DIR)>0
    SPM_CANDIDATES = SPM_DIR;
else
    SPM_CANDIDATES = [fullfile(string(getenv("HOME")),"spm25"), ...
                      fullfile(string(getenv("HOME")),"spm")];
end
subjects = discoverSubjectDirs(NIFTI_ROOT);
if isempty(subjects)
    error("No subject NIfTI directories were found under: %s", NIFTI_ROOT);
end

if ~isfolder(OUT_ROOT), mkdir(OUT_ROOT); end

if exist('spm_vol','file') ~= 2
    for i=1:numel(SPM_CANDIDATES)
        if isfolder(SPM_CANDIDATES(i))
            addpath(char(SPM_CANDIDATES(i)));
            break;
        end
    end
end

if exist('spm_vol','file') ~= 2
    error("SPM not found.");
end

rows = {};

fprintf("============================================================\n");
fprintf("PSE DIFFUSION INPUT AUDIT\n");
fprintf("============================================================\n");

for is = 1:numel(subjects)

    SUB = subjects(is);
    acuteDir = fullfile(NIFTI_ROOT,SUB,"acute");

    fprintf("\n%s\n",SUB);

    if ~isfolder(acuteDir)
        fprintf("  Missing acute directory.\n");
        continue;
    end

    d1 = dir(fullfile(acuteDir,"**","*.nii"));
    d2 = dir(fullfile(acuteDir,"**","*.nii.gz"));
    d = [d1; d2];

    for i=1:numel(d)

        p = string(fullfile(d(i).folder,d(i).name));
        low = lower(p);

        if contains(low,"adc")
            cls = "ADC_EXACT";
        elseif contains(low,"apparent") || contains(low,"coefficient") || ...
               contains(low,"eadc") || contains(low,"trace")
            cls = "ADC_LIKE";
        elseif contains(low,"dwi") || contains(low,"diffus")
            cls = "DWI_OR_DIFFUSION";
        elseif contains(low,"t1")
            cls = "T1";
        else
            cls = "OTHER";
        end

        pRead = p;

        try
            if endsWith(lower(p),".nii.gz")
                tmpDir = fullfile(OUT_ROOT,"_audit_tmp",SUB);
                if ~isfolder(tmpDir), mkdir(tmpDir); end
                files = gunzip(p,tmpDir);
                pRead = string(files{1});
            end

            V = spm_vol(char(pRead));
            dim = V(1).dim;
            vox = sqrt(sum(V(1).mat(1:3,1:3).^2,1));
            nVol = numel(V);

            dimTxt = sprintf("%dx%dx%d",dim(1),dim(2),dim(3));
            voxTxt = sprintf("%.4gx%.4gx%.4g",vox(1),vox(2),vox(3));

        catch
            dimTxt = "";
            voxTxt = "";
            nVol = NaN;
            cls = cls + "_READ_ERROR";
        end

        rows(end+1,:) = {char(SUB),char(cls),char(p),char(string(d(i).name)), ...
                         char(dimTxt),nVol,char(voxTxt)}; %#ok<AGROW>
    end
end

T = cell2table(rows,'VariableNames', ...
    {'Subject','CandidateClass','FullPath','Filename','Dimensions','NVolumes','VoxelSize_mm'});

priority = 9*ones(height(T),1);
priority(T.CandidateClass=="ADC_EXACT") = 1;
priority(T.CandidateClass=="ADC_LIKE") = 2;
priority(T.CandidateClass=="DWI_OR_DIFFUSION") = 3;
priority(T.CandidateClass=="T1") = 4;

T.Priority = priority;
T = sortrows(T,{'Subject','Priority','FullPath'});
T.Priority = [];

outCSV = fullfile(OUT_ROOT,"PSE_DIFFUSION_INPUT_AUDIT.csv");
writetable(T,outCSV);

fprintf("\n============================================================\n");
fprintf("AUDIT COMPLETE\n");
fprintf("CSV: %s\n",outCSV);
fprintf("============================================================\n");

for SUB = subjects
    fprintf("\n--- %s ---\n",SUB);
    disp(T(T.Subject==SUB & T.CandidateClass~="OTHER",:));
end

function subjects = discoverSubjectDirs(parentDir)
% Discover subject directories dynamically.

    d = dir(fullfile(parentDir, "sub-*"));
    d = d([d.isdir]);
    subjects = sort(string({d.name}));
end
