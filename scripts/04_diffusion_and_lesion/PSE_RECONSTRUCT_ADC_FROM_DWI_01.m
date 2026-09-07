%% PSE_RECONSTRUCT_ADC_FROM_DWI_01.m
% Reconstruct ADC maps for acute subjects without a vendor-provided ADC map.
%
% This reproduces the reconstruction used in the PSE acute pilot for
% sub-P003, sub-P004 and sub-P005 after source-DICOM audit confirmed that:
%   DWI_002.nii = b=0 s/mm^2
%   DWI_001.nii = b=1000 s/mm^2
%
% Formula:
%   ADC = -log(S_b1000 / S_b0) / 1000
%
% ADC is stored in mm^2/s. Invalid voxels (non-finite or non-positive source
% signal) are set to zero. Negative finite ADC values are retained here; the
% final NVAUTO preprocessing later clips negative intensities to zero.
%
% IMPORTANT: do not reuse the DWI_001/DWI_002 mapping on another dataset
% without first verifying the b-values from source DICOM metadata.
%
% Requirements: MATLAB + SPM.

clear; clc;

ROOT = string(getenv("PSE_ROOT"));
if strlength(ROOT)==0, ROOT = fullfile(string(getenv("HOME")),"PSE"); end

SPM_DIR = string(getenv("SPM_DIR"));
if strlength(SPM_DIR)>0 && isfolder(SPM_DIR)
    addpath(char(SPM_DIR));
end
if exist('spm_vol','file') ~= 2
    error('SPM not found. Set SPM_DIR or add SPM to the MATLAB path.');
end

NIFTI_ROOT = fullfile(ROOT,"derivatives","nifti");
OUT_ROOT = fullfile(ROOT,"derivatives","adc_reconstructed");
if ~isfolder(OUT_ROOT), mkdir(OUT_ROOT); end

subjects = ["sub-P003","sub-P004","sub-P005"];
bValue = 1000; % s/mm^2
rows = {};

fprintf("============================================================\n");
fprintf("PSE ADC RECONSTRUCTION FROM b0/b1000 DWI\n");
fprintf("============================================================\n");

for i = 1:numel(subjects)
    SUB = subjects(i);
    dwiDir = fullfile(NIFTI_ROOT,SUB,"acute","DWI");
    b0File = fullfile(dwiDir,"DWI_002.nii");
    b1000File = fullfile(dwiDir,"DWI_001.nii");
    outDir = fullfile(OUT_ROOT,SUB,"acute");
    if ~isfolder(outDir), mkdir(outDir); end
    outFile = fullfile(outDir,"ADC_reconstructed_mm2_s.nii");
    qcFile = fullfile(outDir,SUB + "_ADC_reconstruction_QC.png");

    fprintf("\n[%d/%d] %s\n",i,numel(subjects),SUB);

    try
        if ~isfile(b0File) || ~isfile(b1000File)
            error('Missing DWI_001.nii or DWI_002.nii.');
        end

        V0 = spm_vol(char(b0File)); V0 = V0(1);
        Vb = spm_vol(char(b1000File)); Vb = Vb(1);
        if any(V0.dim ~= Vb.dim) || max(abs(V0.mat(:)-Vb.mat(:))) > 1e-4
            error('b0 and b1000 NIfTI grids differ; inspect before reconstruction.');
        end

        S0 = spm_read_vols(V0);
        Sb = spm_read_vols(Vb);

        valid = isfinite(S0) & isfinite(Sb) & S0>0 & Sb>0;
        ADC = zeros(size(S0),'single');
        ADC(valid) = single(-log(Sb(valid)./S0(valid))/bValue);
        ADC(~isfinite(ADC)) = 0;

        Vout = V0;
        Vout.fname = char(outFile);
        Vout.dt = [16 0]; % float32
        Vout.pinfo = [1;0;0];
        spm_write_vol(Vout,double(ADC));

        s0Vals = S0(valid); sbVals = Sb(valid);
        adcPos = double(ADC(valid & ADC>0));
        med0 = median(s0Vals,'omitnan');
        medb = median(sbVals,'omitnan');
        medADC = median(adcPos,'omitnan');
        p05 = prctile(adcPos,5);
        p95 = prctile(adcPos,95);

        makeQC(S0,Sb,double(ADC),qcFile,SUB);

        rows(end+1,:) = {char(SUB),'RECONSTRUCTED',char(outFile), ... %#ok<AGROW>
            char(b0File),char(b1000File),med0,medb,med0/medb, ...
            medADC,p05,p95,'b0=DWI_002.nii; b1000=DWI_001.nii'};

        fprintf('  Median reconstructed ADC: %.9f mm^2/s\n',medADC);
        fprintf('  Output: %s\n',outFile);

    catch ME
        rows(end+1,:) = {char(SUB),'ERROR',char(outFile),char(b0File), ... %#ok<AGROW>
            char(b1000File),NaN,NaN,NaN,NaN,NaN,NaN,char(ME.message)};
        fprintf('  ERROR: %s\n',ME.message);
    end
end

T = cell2table(rows,'VariableNames',{
    'Subject','Status','ADCFile','B0File','B1000File', ...
    'MedianSignal_B0','MedianSignal_B1000','MedianSignalRatio_B0_to_B1000', ...
    'MedianADC_mm2_s','P05_ADC_mm2_s','P95_ADC_mm2_s','Message'});

outCSV = fullfile(OUT_ROOT,'PSE_ADC_RECONSTRUCTION_STATUS.csv');
writetable(T,outCSV);

disp(T);
fprintf('\nStatus CSV: %s\n',outCSV);

%% ------------------------------------------------------------------------
function makeQC(S0,Sb,ADC,outFile,SUB)
    mask = ADC>0 & isfinite(ADC);
    z = find(squeeze(any(any(mask,1),2)));
    if isempty(z)
        zlist = round(linspace(1,size(ADC,3),5));
    else
        zlist = unique(round(linspace(min(z),max(z),5)));
    end

    f = figure('Visible','off','Color','w','Position',[100 100 1500 720]);
    tl = tiledlayout(3,numel(zlist),'Padding','compact','TileSpacing','compact');
    for j=1:numel(zlist)
        zz=zlist(j);
        nexttile; imagesc(rot90(S0(:,:,zz))); axis image off; colormap gray;
        title(sprintf('b0 z=%d',zz));
        nexttile; imagesc(rot90(Sb(:,:,zz))); axis image off; colormap gray;
        title('b1000');
        nexttile; imagesc(rot90(ADC(:,:,zz)),[0 0.003]); axis image off; colormap gray;
        title('ADC (mm^2/s)');
    end
    title(tl,SUB + " reconstructed ADC");
    exportgraphics(f,char(outFile),'Resolution',180);
    close(f);
end
