%% PSE_DIFFUSION_COREG_T1_02.m
% Revised diffusion -> T1 coregistration.
%
% Changes from v01
% ----------------
% 1) Uses original ADC when available.
% 2) Falls back to a reconstructed ADC when a vendor ADC is absent.
% 3) Explicitly identifies the high-b DWI from DWI_001 / DWI_002 using the
%    lower positive-signal median, instead of selecting an arbitrary DWI.
% 4) Records ADC source and ADC/DWI center distance before registration.
%
% Strategy
% --------
% - T1 = reference
% - ADC = source for rigid NMI estimate
% - The same rigid transform is applied to copies of ADC and high-b DWI
%   headers, preserving their native voxel grids and mutual scanner-space
%   relationship.
%
% OUTPUT
% ------
% <PSE_ROOT>/derivatives/coreg_diffusion_v02/
%
% Requires SPM.

clear; clc;

ROOT = string(getenv("PSE_ROOT"));
if strlength(ROOT)==0, ROOT = fullfile(string(getenv("HOME")),"PSE"); end
NIFTI_ROOT = fullfile(ROOT,"derivatives","nifti");
ADC_RECON_ROOT = fullfile(ROOT,"derivatives","adc_reconstructed");
OUT_ROOT = fullfile(ROOT,"derivatives","coreg_diffusion_v02");

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

if exist('spm_coreg','file') ~= 2
    for i=1:numel(SPM_CANDIDATES)
        if isfolder(SPM_CANDIDATES(i))
            addpath(char(SPM_CANDIDATES(i)));
            break;
        end
    end
end
if exist('spm_coreg','file') ~= 2, error("SPM not found."); end

spm('defaults','PET');

flags = struct( ...
    'sep',[4 2], ...
    'cost_fun','nmi', ...
    'tol',[0.0200 0.0200 0.0200 0.0010 0.0010 0.0010 ...
           0.0100 0.0100 0.0100 0.0010 0.0010 0.0010], ...
    'fwhm',[7 7]);

rows={};

fprintf("============================================================\n");
fprintf("PSE DIFFUSION -> T1 COREGISTRATION v02\n");
fprintf("SPM: %s\n",spm('Ver'));
fprintf("============================================================\n");

for is=1:numel(subjects)

    SUB=subjects(is);
    acuteDir=fullfile(NIFTI_ROOT,SUB,"acute");

    fprintf("\n[%d/%d] %s\n",is,numel(subjects),SUB);

    t1File=fullfile(acuteDir,"T1","T1.nii");
    adcOriginal=fullfile(acuteDir,"ADC","ADC.nii");
    adcRecon=fullfile(ADC_RECON_ROOT,SUB,"acute","ADC_reconstructed_mm2_s.nii");
    dwiDir=fullfile(acuteDir,"DWI");

    if ~isfile(t1File)
        rows(end+1,:)={char(SUB),"MISSING_T1","","","","",NaN,""}; %#ok<AGROW>
        continue;
    end

    if isfile(adcOriginal)
        adcFile=string(adcOriginal);
        adcSource="ORIGINAL";
    elseif isfile(adcRecon)
        adcFile=string(adcRecon);
        adcSource="RECONSTRUCTED";
    else
        rows(end+1,:)={char(SUB),"MISSING_ADC",char(t1File),"","","",NaN,""}; %#ok<AGROW>
        fprintf("  Missing both original and reconstructed ADC.\n");
        continue;
    end

    [dwiHigh,dwiMsg]=findHighBDWI(dwiDir);

    if strlength(dwiHigh)==0
        rows(end+1,:)={char(SUB),"MISSING_DWI",char(t1File),char(adcFile),char(adcSource),"",NaN,char(dwiMsg)}; %#ok<AGROW>
        continue;
    end

    outDir=fullfile(OUT_ROOT,SUB,"acute");
    if ~isfolder(outDir), mkdir(outDir); end

    try
        T1out=copyfileNii(t1File,fullfile(outDir,"T1_reference.nii"));
        ADCout=copyfileNii(adcFile,fullfile(outDir,"ADC_coreg_T1.nii"));
        DWIout=copyfileNii(dwiHigh,fullfile(outDir,"DWI_b1000_coreg_T1.nii"));

        Vref=spm_vol(char(T1out)); Vref=Vref(1);
        Vadc=spm_vol(char(ADCout)); Vadc=Vadc(1);
        Vdwi=spm_vol(char(DWIout)); Vdwi=Vdwi(1);

        centerADC=worldCenter(Vadc);
        centerDWI=worldCenter(Vdwi);
        centerDist=norm(centerADC-centerDWI);

        x=spm_coreg(Vref,Vadc,flags);
        M=spm_matrix(x);

        applyHeaderTransform(ADCout,M);
        applyHeaderTransform(DWIout,M);

        if centerDist>10
            status="PROCESSED_WARN_ADC_DWI_CENTER";
        else
            status="PROCESSED";
        end

        msg=sprintf("ADC=%s; %s; x=[%s]",adcSource,dwiMsg,sprintf(' %.6g',x));

        rows(end+1,:)={char(SUB),char(status),char(T1out),char(ADCout), ...
            char(adcSource),char(DWIout),centerDist,char(msg)}; %#ok<AGROW>

        fprintf("  ADC source: %s\n",adcSource);
        fprintf("  High-b DWI: %s\n",dwiHigh);
        fprintf("  ADC/DWI center distance before transform: %.2f mm\n",centerDist);
        fprintf("  x=[%s]\n",sprintf(' %.6g',x));

    catch ME
        rows(end+1,:)={char(SUB),"ERROR",char(t1File),char(adcFile), ...
            char(adcSource),char(dwiHigh),NaN,char(ME.message)}; %#ok<AGROW>
        fprintf("  ERROR: %s\n",ME.message);
    end
end

Tstatus=cell2table(rows,'VariableNames', ...
    {'Subject','Status','T1File','ADCFile','ADCSource','DWIFile', ...
     'ADC_DWI_CenterDistance_mm','Message'});

outCSV=fullfile(OUT_ROOT,"PSE_DIFFUSION_COREG_T1_V02_STATUS.csv");
writetable(Tstatus,outCSV);

fprintf("\n============================================================\n");
fprintf("DIFFUSION COREG v02 COMPLETE\n");
fprintf("Status: %s\n",outCSV);
fprintf("============================================================\n");
disp(Tstatus);

%% ------------------------------------------------------------------------
function [f,msg]=findHighBDWI(dwiDir)
    f=""; msg="";
    d=dir(fullfile(dwiDir,"DWI_*.nii"));
    if isempty(d)
        msg="no DWI_*.nii found"; return;
    end

    paths=strings(numel(d),1);
    meds=nan(numel(d),1);

    for k=1:numel(d)
        paths(k)=string(fullfile(d(k).folder,d(k).name));
        V=spm_vol(char(paths(k)));
        if numel(V)~=1
            continue;
        end
        Y=spm_read_vols(V);
        vals=Y(isfinite(Y) & Y>0);
        if ~isempty(vals), meds(k)=median(vals); end
    end

    good=isfinite(meds);
    if ~any(good)
        msg="no readable positive-signal DWI"; return;
    end

    idx=find(good);
    [~,j]=min(meds(idx));
    sel=idx(j);
    f=paths(sel);

    if nnz(good)>=2
        sorted=sort(meds(good));
        msg=sprintf("lower-median DWI selected as high-b; medians=[%s]", ...
            sprintf(' %.4g',sorted));
    else
        msg="only one readable DWI candidate";
    end
end

function out=copyfileNii(inFile,outFile)
    copyfile(char(inFile),char(outFile),'f');
    out=string(outFile);
end

function c=worldCenter(V)
    v=(V.dim(:)+1)/2;
    w=V.mat*[v;1];
    c=w(1:3);
end

function applyHeaderTransform(niftiFile,M)
    niftiFile=char(niftiFile);
    V=spm_vol(niftiFile);
    for k=1:numel(V)
        spec=sprintf('%s,%d',niftiFile,k);
        oldMat=spm_get_space(spec);
        spm_get_space(spec,M\oldMat);
    end
end

function subjects = discoverSubjectDirs(parentDir)
% Discover subject directories dynamically.

    d = dir(fullfile(parentDir, "sub-*"));
    d = d([d.isdir]);
    subjects = sort(string({d.name}));
end
