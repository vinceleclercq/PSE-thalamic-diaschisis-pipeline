%% PSE_COREG_VISUAL_QC_01.m
% Interactive visual QC of acute rigid coregistration using SPM Check Reg.
%
% Non-destructive: opens images only.
%
% For each subject it displays, when available:
%   1. T1 native reference
%   2. PET_FDG coregistered copy
%   3. ASL_CBF coregistered copy
%   4. CTP_REF coregistered copy
%
% Use the linked crosshair in SPM to judge anatomical alignment.
% Press ENTER in the MATLAB Command Window to move to the next subject.
% Type q then ENTER to stop.

clear; clc;

rootDir = getenv('PSE_ROOT');
if isempty(rootDir), rootDir = fullfile(getenv('HOME'),'PSE'); end
niftiRoot = fullfile(rootDir,'derivatives','nifti');
coregRoot = fullfile(rootDir,'derivatives','coreg');

if isempty(which('spm_check_registration'))
    error('SPM is not available on the MATLAB path.');
end

spm('defaults','fmri');

subjects = dir(fullfile(niftiRoot,'sub-P*'));
subjects = subjects([subjects.isdir]);
[~,ord] = sort({subjects.name});
subjects = subjects(ord);

fprintf('\n============================================================\n');
fprintf('PSE COREGISTRATION VISUAL QC\n');
fprintf('Pay particular attention to:\n');
fprintf('  - PET_FDG in sub-P004\n');
fprintf('  - ASL_CBF in sub-P001 and sub-P004\n');
fprintf('  - CTP_REF in ALL subjects\n');
fprintf('============================================================\n\n');

for s = 1:numel(subjects)

    subject = subjects(s).name;

    t1 = findSingleNifti(fullfile(niftiRoot,subject,'acute','T1'));

    pet = string(fullfile(coregRoot,subject,'acute','PET_FDG','PET_FDG_coreg.nii'));
    asl = string(fullfile(coregRoot,subject,'acute','ASL_CBF','ASL_CBF_coreg.nii'));
    ctp = string(fullfile(coregRoot,subject,'acute','CTP','CTP_REF_coreg.nii'));

    paths = strings(0,1);
    labels = strings(0,1);

    if strlength(t1)>0 && isfile(t1)
        paths(end+1,1) = t1; %#ok<SAGROW>
        labels(end+1,1) = "T1"; %#ok<SAGROW>
    end
    if isfile(pet)
        paths(end+1,1) = pet; %#ok<SAGROW>
        labels(end+1,1) = "PET_FDG"; %#ok<SAGROW>
    end
    if isfile(asl)
        paths(end+1,1) = asl; %#ok<SAGROW>
        labels(end+1,1) = "ASL_CBF"; %#ok<SAGROW>
    end
    if isfile(ctp)
        paths(end+1,1) = ctp; %#ok<SAGROW>
        labels(end+1,1) = "CTP_REF"; %#ok<SAGROW>
    end

    fprintf('\n------------------------------------------------------------\n');
    fprintf('%s\n',subject);
    fprintf('Displayed order: %s\n',strjoin(labels,' | '));
    fprintf('------------------------------------------------------------\n');

    if numel(paths)<2
        fprintf('Not enough images available for registration QC.\n');
        continue;
    end

    spm_check_registration(char(paths));

    response = input(sprintf( ...
        '%s: inspect alignment in SPM, then press ENTER for next subject (q to quit): ', ...
        subject),'s');

    if strcmpi(strtrim(response),'q')
        break;
    end
end

fprintf('\nVisual QC session finished.\n');

function nii = findSingleNifti(folder)
    nii = "";
    if ~isfolder(folder)
        return;
    end
    d = dir(fullfile(folder,'*.nii'));
    if numel(d)==1
        nii = string(fullfile(d(1).folder,d(1).name));
    end
end
