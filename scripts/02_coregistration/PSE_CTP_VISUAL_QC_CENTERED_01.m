%% PSE_CTP_VISUAL_QC_CENTERED_01.m
% Interactive visual QC of corrected CTP -> T1 coregistration.
%
% Displays, when available:
%   1. T1 native reference
%   2. corrected CTP_REF
%   3. corrected CTP_CBF
%
% Non-destructive. Press ENTER in MATLAB Command Window to advance.
% Type q then ENTER to stop.

clear; clc;

rootDir = getenv('PSE_ROOT');
if isempty(rootDir), rootDir = fullfile(getenv('HOME'),'PSE'); end
niftiRoot = fullfile(rootDir,'derivatives','nifti');
ctpRoot = fullfile(rootDir,'derivatives','coreg_ctp_centered');

if isempty(which('spm_check_registration'))
    error('SPM is not available on the MATLAB path.');
end

spm('defaults','fmri');

subjects = dir(fullfile(niftiRoot,'sub-P*'));
subjects = subjects([subjects.isdir]);
[~,ord] = sort({subjects.name});
subjects = subjects(ord);

fprintf('\n============================================================\n');
fprintf('PSE CTP CENTERED COREGISTRATION - VISUAL QC\n');
fprintf('Displayed order: T1 | corrected CTP_REF | corrected CTP_CBF\n');
fprintf('============================================================\n\n');

for s = 1:numel(subjects)

    subject = subjects(s).name;

    t1 = findSingleNifti(fullfile(niftiRoot,subject,'acute','T1'));
    ctpRef = string(fullfile(ctpRoot,subject,'acute','CTP','CTP_REF_coreg_centered.nii'));
    ctpCbf = string(fullfile(ctpRoot,subject,'acute','CTP','CTP_CBF_coreg_centered.nii'));

    paths = strings(0,1);
    labels = strings(0,1);

    if strlength(t1)>0 && isfile(t1)
        paths(end+1,1)=t1; %#ok<SAGROW>
        labels(end+1,1)="T1"; %#ok<SAGROW>
    end

    if isfile(ctpRef)
        paths(end+1,1)=ctpRef; %#ok<SAGROW>
        labels(end+1,1)="CTP_REF"; %#ok<SAGROW>
    end

    if isfile(ctpCbf)
        paths(end+1,1)=ctpCbf; %#ok<SAGROW>
        labels(end+1,1)="CTP_CBF"; %#ok<SAGROW>
    end

    fprintf('\n------------------------------------------------------------\n');
    fprintf('%s\n',subject);
    fprintf('Displayed: %s\n',strjoin(labels,' | '));
    fprintf('------------------------------------------------------------\n');

    if numel(paths)<2
        fprintf('Not enough images available for QC.\n');
        continue;
    end

    spm_check_registration(char(paths));

    response = input(sprintf( ...
        '%s: inspect alignment, press ENTER for next subject (q to quit): ', ...
        subject),'s');

    if strcmpi(strtrim(response),'q')
        break;
    end
end

fprintf('\nVisual QC finished.\n');

function nii = findSingleNifti(folder)
    nii = "";
    if ~isfolder(folder), return; end
    d = dir(fullfile(folder,'*.nii'));
    if numel(d)==1
        nii = string(fullfile(d(1).folder,d(1).name));
    end
end
