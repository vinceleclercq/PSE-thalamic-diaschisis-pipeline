%% PSE_CTP_OVERLAY_QC_01.m
% PSE PET-MR project - overlay QC for corrected CTP -> T1 registration.
%
% Purpose:
%   Create TEMPORARY resliced copies of CTP_REF in T1 space ONLY for visual
%   quality control. Quantitative CTP maps remain untouched in native space.
%
% For each subject with CTP_REF:
%   1) copy corrected CTP_REF to derivatives/qc_ctp_overlay/
%   2) reslice that COPY to the T1 grid with SPM
%   3) create a PNG with T1 background + CT edge contours in 3 orthogonal views
%
% Outputs:
%   <PSE_ROOT>/derivatives/qc_ctp_overlay/
%       sub-Pxxx/
%           CTP_REF_resliced_QC.nii
%           sub-Pxxx_CTP_T1_overlay.png
%
% Requirements:
%   MATLAB + Image Processing Toolbox + SPM

clear; clc;

%% ----------------------------- SETTINGS --------------------------------

rootDir = getenv('PSE_ROOT');
if isempty(rootDir), rootDir = fullfile(getenv('HOME'),'PSE'); end

niftiRoot = fullfile(rootDir,'derivatives','nifti');
ctpRoot   = fullfile(rootDir,'derivatives','coreg_ctp_centered');
qcRoot    = fullfile(rootDir,'derivatives','qc_ctp_overlay');

if isempty(which('spm')) || isempty(which('spm_reslice'))
    error('SPM is not available on the MATLAB path.');
end

if isempty(which('edge'))
    error('Image Processing Toolbox is required for edge overlays.');
end

if ~isfolder(qcRoot)
    mkdir(qcRoot);
end

spm('defaults','fmri');
try
    spm_get_defaults('cmdline',true);
catch
end

subjects = dir(fullfile(niftiRoot,'sub-P*'));
subjects = subjects([subjects.isdir]);

[~,ord] = sort({subjects.name});
subjects = subjects(ord);

fprintf('\n============================================================\n');
fprintf('PSE CTP -> T1 OVERLAY QC\n');
fprintf('Temporary CT reslicing for QC only.\n');
fprintf('Quantitative CTP maps are NOT resampled.\n');
fprintf('============================================================\n\n');

for s = 1:numel(subjects)

    subject = subjects(s).name;

    t1 = findSingleNifti(fullfile(niftiRoot,subject,'acute','T1'));
    ctp = string(fullfile(ctpRoot,subject,'acute','CTP','CTP_REF_coreg_centered.nii'));

    if strlength(t1)==0 || ~isfile(ctp)
        fprintf('%s : skipped (T1 or corrected CTP_REF missing)\n',subject);
        continue;
    end

    subjectQcDir = fullfile(qcRoot,subject);
    if ~isfolder(subjectQcDir)
        mkdir(subjectQcDir);
    end

    % Copy CTP_REF into QC folder so that reslicing NEVER touches the
    % corrected registration output used by the pipeline.
    ctpQcCopy = fullfile(subjectQcDir,'CTP_REF_QCcopy.nii');
    if isfile(ctpQcCopy)
        delete(ctpQcCopy);
    end
    copyfile(ctp,ctpQcCopy);

    % SPM reslice: first image = reference grid (T1), second = image to reslice.
    flags = struct( ...
        'interp',1, ...      % trilinear; QC structural image only
        'wrap',[0 0 0], ...
        'mask',0, ...
        'mean',0, ...
        'which',1, ...
        'prefix','r');

    P = char([string(t1); string(ctpQcCopy)]);
    spm_reslice(P,flags);

    reslicedDefault = fullfile(subjectQcDir,'rCTP_REF_QCcopy.nii');

    if ~isfile(reslicedDefault)
        fprintf('%s : ERROR - SPM did not create resliced CT\n',subject);
        continue;
    end

    finalResliced = fullfile(subjectQcDir,'CTP_REF_resliced_QC.nii');
    if isfile(finalResliced)
        delete(finalResliced);
    end
    movefile(reslicedDefault,finalResliced);

    % Read both images in exactly the same T1 grid.
    Vt1 = spm_vol(char(t1));
    Vct = spm_vol(char(finalResliced));

    T1 = spm_read_vols(Vt1);
    CT = spm_read_vols(Vct);

    if ~isequal(size(T1),size(CT))
        fprintf('%s : ERROR - grids do not match after reslicing\n',subject);
        continue;
    end

    % Select slices near the geometric centre of the common grid.
    sz = size(T1);
    ix = round(sz(1)/2);
    iy = round(sz(2)/2);
    iz = round(sz(3)/2);

    % Generate orthogonal views.
    t1Sag = squeeze(T1(ix,:,:))';
    ctSag = squeeze(CT(ix,:,:))';

    t1Cor = squeeze(T1(:,iy,:))';
    ctCor = squeeze(CT(:,iy,:))';

    t1Axi = squeeze(T1(:,:,iz))';
    ctAxi = squeeze(CT(:,:,iz))';

    % Standard radiological-like display orientation for easy viewing.
    t1Sag = flipud(t1Sag); ctSag = flipud(ctSag);
    t1Cor = flipud(t1Cor); ctCor = flipud(ctCor);
    t1Axi = flipud(t1Axi); ctAxi = flipud(ctAxi);

    pngPath = fullfile(subjectQcDir,subject + "_CTP_T1_overlay.png");

    fig = figure('Visible','off','Color','w','Position',[100 100 1500 500]);

    makeOverlayPanel(t1Sag,ctSag,1,'Sagittal');
    makeOverlayPanel(t1Cor,ctCor,2,'Coronal');
    makeOverlayPanel(t1Axi,ctAxi,3,'Axial');

    sgtitle(sprintf('%s - T1 background with CTP_REF edges (QC only)',subject), ...
        'Interpreter','none');

    exportgraphics(fig,pngPath,'Resolution',160);
    close(fig);

    fprintf('%s : created %s\n',subject,pngPath);
end

fprintf('\n============================================================\n');
fprintf('OVERLAY QC COMPLETE\n');
fprintf('Open PNG files under:\n%s\n',qcRoot);
fprintf('============================================================\n');

%% ============================= FUNCTIONS ================================

function nii = findSingleNifti(folder)
    nii = "";
    if ~isfolder(folder), return; end
    d = dir(fullfile(folder,'*.nii'));
    if numel(d)==1
        nii = string(fullfile(d(1).folder,d(1).name));
    end
end

function makeOverlayPanel(t1Slice,ctSlice,panelNumber,panelTitle)

    subplot(1,3,panelNumber);

    t1Norm = robustNormalize(t1Slice);

    imagesc(t1Norm);
    axis image off;
    colormap(gca,gray);
    hold on;

    % CT brain/skull edges. Use robust CT window before edge detection.
    ctNorm = robustNormalize(ctSlice);

    % Ignore completely empty background.
    nonzero = abs(ctSlice) > 0;
    ctNorm(~nonzero) = 0;

    bw = edge(ctNorm,'Canny');

    % Remove image-border artefacts.
    bw([1 end],:) = false;
    bw(:,[1 end]) = false;

    contour(bw,[0.5 0.5],'r','LineWidth',0.8);

    title(panelTitle);
    hold off;
end

function y = robustNormalize(x)

    x = double(x);
    vals = x(isfinite(x));

    if isempty(vals)
        y = zeros(size(x));
        return;
    end

    lo = localPercentile(vals,1);
    hi = localPercentile(vals,99);

    if ~isfinite(lo) || ~isfinite(hi) || hi<=lo
        lo = min(vals);
        hi = max(vals);
    end

    if hi<=lo
        y = zeros(size(x));
        return;
    end

    x(~isfinite(x)) = lo;
    y = (x-lo)./(hi-lo);
    y = min(max(y,0),1);
end

function p = localPercentile(x,pct)

    x = sort(x(:));
    n = numel(x);

    if n==0
        p = NaN;
        return;
    end

    idx = 1 + (n-1)*pct/100;
    i0 = floor(idx);
    i1 = ceil(idx);

    if i0==i1
        p = x(i0);
    else
        p = x(i0) + (idx-i0)*(x(i1)-x(i0));
    end
end
