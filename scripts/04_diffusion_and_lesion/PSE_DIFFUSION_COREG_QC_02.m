%% PSE_DIFFUSION_COREG_QC_02.m
% Visual QC for diffusion -> T1 coregistration v02.
%
% Reads:
% <PSE_ROOT>/derivatives/coreg_diffusion_v02/
%   sub-P00X/acute/
%       T1_reference.nii
%       ADC_coreg_T1.nii
%       DWI_b1000_coreg_T1.nii
%
% For each subject, the DWI grid is used as the QC grid.
% T1 and ADC are sampled into the DWI grid in world coordinates ONLY FOR
% DISPLAY. No files are modified or written except PNG QC images.
%
% Figure rows:
%   1) T1 resampled to DWI grid
%   2) b1000 DWI
%   3) ADC resampled to DWI grid
%   4) T1/DWI RGB overlay (T1=red, DWI=green)
%
% Output:
% <PSE_ROOT>/derivatives/coreg_diffusion_v02_qc/
%   sub-P00X_diffusion_T1_coreg_v02_QC.png
%   PSE_DIFFUSION_COREG_V02_QC_STATUS.csv
%
% Requires SPM.

clear; clc; close all;

ROOT = string(getenv("PSE_ROOT"));
if strlength(ROOT)==0, ROOT = fullfile(string(getenv("HOME")),"PSE"); end
IN_ROOT = fullfile(ROOT,"derivatives","coreg_diffusion_v02");
OUT_ROOT = fullfile(ROOT,"derivatives","coreg_diffusion_v02_qc");

SPM_DIR = string(getenv("SPM_DIR"));
if strlength(SPM_DIR)>0
    SPM_CANDIDATES = SPM_DIR;
else
    SPM_CANDIDATES = [fullfile(string(getenv("HOME")),"spm25"), ...
                      fullfile(string(getenv("HOME")),"spm")];
end
subjects = discoverSubjectDirs(IN_ROOT);
if isempty(subjects)
    error("No diffusion-coregistration subject directories were found under: %s", IN_ROOT);
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
if exist('spm_vol','file') ~= 2, error("SPM not found."); end

rows={};

fprintf("============================================================\n");
fprintf("PSE DIFFUSION -> T1 COREG v02 VISUAL QC\n");
fprintf("============================================================\n");

for is=1:numel(subjects)

    SUB=subjects(is);
    inDir=fullfile(IN_ROOT,SUB,"acute");

    t1File=fullfile(inDir,"T1_reference.nii");
    adcFile=fullfile(inDir,"ADC_coreg_T1.nii");
    dwiFile=fullfile(inDir,"DWI_b1000_coreg_T1.nii");

    fprintf("\n[%d/%d] %s\n",is,numel(subjects),SUB);

    if ~isfile(t1File) || ~isfile(adcFile) || ~isfile(dwiFile)
        rows(end+1,:)={char(SUB),"MISSING_INPUT",""}; %#ok<AGROW>
        fprintf("  Missing input.\n");
        continue;
    end

    try
        Vt=spm_vol(char(t1File)); Vt=Vt(1);
        Va=spm_vol(char(adcFile)); Va=Va(1);
        Vd=spm_vol(char(dwiFile)); Vd=Vd(1);

        D=spm_read_vols(Vd);
        T1onDWI=resampleToReference(Vt,Vd,1);
        ADConDWI=resampleToReference(Va,Vd,1);

        % Choose slices with largest DWI positive support.
        support=isfinite(D) & D>0;
        sliceCounts=squeeze(sum(sum(support,1),2));

        good=find(sliceCounts > 0.20*max(sliceCounts));

        if numel(good)>=6
            zList=unique(round(linspace(min(good),max(good),6)));
        else
            zList=unique(round(linspace(1,size(D,3),6)));
        end

        fig=figure('Visible','off','Color','w','Position',[100 100 1700 1150]);
        tiledlayout(fig,4,numel(zList),'Padding','compact','TileSpacing','compact');

        wt=robustWindow(T1onDWI);
        wd=robustWindow(D);
        wa=robustWindow(ADConDWI);

        for iz=1:numel(zList)
            z=zList(iz);

            t=rot90(T1onDWI(:,:,z));
            d=rot90(D(:,:,z));
            a=rot90(ADConDWI(:,:,z));

            nexttile(iz);
            imagesc(t,wt); axis image off; colormap gray;
            title(sprintf('T1 on DWI | z=%d',z));

            nexttile(numel(zList)+iz);
            imagesc(d,wd); axis image off; colormap gray;
            title(sprintf('b1000 DWI | z=%d',z));

            nexttile(2*numel(zList)+iz);
            imagesc(a,wa); axis image off; colormap gray;
            title(sprintf('ADC on DWI | z=%d',z));

            nexttile(3*numel(zList)+iz);
            rgb=makeRGOverlay(t,d,wt,wd);
            image(rgb); axis image off;
            title('Overlay T1=R / DWI=G');
        end

        sgtitle(sprintf('%s | diffusion -> T1 coregistration v02 QC',SUB), ...
            'Interpreter','none');

        outPng=fullfile(OUT_ROOT,sprintf('%s_diffusion_T1_coreg_v02_QC.png',SUB));
        exportgraphics(fig,outPng,'Resolution',180);
        close(fig);

        rows(end+1,:)={char(SUB),"QC_CREATED",char(outPng)}; %#ok<AGROW>
        fprintf("  QC: %s\n",outPng);

    catch ME
        if exist('fig','var') && isgraphics(fig), close(fig); end
        rows(end+1,:)={char(SUB),"ERROR",char(ME.message)}; %#ok<AGROW>
        fprintf("  ERROR: %s\n",ME.message);
    end
end

Tstatus=cell2table(rows,'VariableNames',{'Subject','Status','MessageOrFile'});
outCSV=fullfile(OUT_ROOT,'PSE_DIFFUSION_COREG_V02_QC_STATUS.csv');
writetable(Tstatus,outCSV);

fprintf("\n============================================================\n");
fprintf("COREG v02 QC COMPLETE\n");
fprintf("Status: %s\n",outCSV);
fprintf("============================================================\n");
disp(Tstatus);

%% ------------------------------------------------------------------------
function Yout=resampleToReference(Vsrc,Vref,interp)
    dim=Vref.dim;
    Yout=nan(dim);
    invSrc=inv(Vsrc.mat);

    for z=1:dim(3)
        [x,y]=ndgrid(1:dim(1),1:dim(2));
        zz=z*ones(size(x));

        voxRef=[x(:)';y(:)';zz(:)';ones(1,numel(x))];
        world=Vref.mat*voxRef;
        voxSrc=invSrc*world;

        vals=spm_sample_vol(Vsrc,voxSrc(1,:),voxSrc(2,:),voxSrc(3,:),interp);
        Yout(:,:,z)=reshape(vals,dim(1),dim(2));
    end
end

function w=robustWindow(X)
    v=X(isfinite(X));
    if isempty(v), w=[0 1]; return; end
    v=v(v~=0);
    if isempty(v), w=[0 1]; return; end
    v=sort(v);
    n=numel(v);
    lo=v(max(1,round(0.02*n)));
    hi=v(min(n,round(0.98*n)));
    if hi<=lo, hi=lo+1; end
    w=[lo hi];
end

function rgb=makeRGOverlay(T,D,wt,wd)
    tn=(T-wt(1))/(wt(2)-wt(1));
    dn=(D-wd(1))/(wd(2)-wd(1));

    tn=max(0,min(1,tn));
    dn=max(0,min(1,dn));

    rgb=zeros([size(T),3]);
    rgb(:,:,1)=tn;
    rgb(:,:,2)=dn;
    rgb(:,:,3)=0.15*(tn+dn);

    bad=~isfinite(T) & ~isfinite(D);
    for c=1:3
        tmp=rgb(:,:,c);
        tmp(bad)=0;
        rgb(:,:,c)=tmp;
    end
end

function subjects = discoverSubjectDirs(parentDir)
% Discover subject directories dynamically from the processing output.

    d = dir(fullfile(parentDir, "sub-*"));
    d = d([d.isdir]);
    subjects = sort(string({d.name}));
end
