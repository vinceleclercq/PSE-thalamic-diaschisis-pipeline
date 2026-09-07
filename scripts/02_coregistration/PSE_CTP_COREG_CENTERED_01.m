%% PSE_CTP_COREG_CENTERED_01.m
% PSE PET-MR project - corrected CT perfusion coregistration to T1.
%
% Purpose:
%   The first CT->T1 registration attempt failed because CT and T1 did not
%   sufficiently overlap in world coordinates. This version first aligns
%   the geometric centres of the CT and T1 volumes, then runs rigid NMI
%   coregistration with SPM.
%
% IMPORTANT:
%   - Source NIfTI files are NEVER modified.
%   - Quantitative CTP maps are NOT resampled.
%   - The same rigid transform estimated from CTP_REF is applied to the
%     affine headers of CTP_CBF / CTP_CBV / CTP_TMAX / CTP_MTT / CTP_TTP.
%   - Outputs go to a NEW folder:
%       derivatives/coreg_ctp_centered/
%
% Requirements:
%   MATLAB + SPM

clear; clc;

%% ----------------------------- SETTINGS --------------------------------

rootDir = getenv('PSE_ROOT');
if isempty(rootDir), rootDir = fullfile(getenv('HOME'),'PSE'); end
sessionName = 'acute';

niftiRoot = fullfile(rootDir,'derivatives','nifti');
outRoot   = fullfile(rootDir,'derivatives','coreg_ctp_centered');
logRoot   = fullfile(rootDir,'derivatives','coreg_ctp_centered_logs');

overwriteExisting = true;

ctpMaps = {'CTP_CBF','CTP_CBV','CTP_TMAX','CTP_MTT','CTP_TTP'};

%% ---------------------------- PREFLIGHT --------------------------------

if isempty(which('spm')) || isempty(which('spm_coreg')) || isempty(which('spm_get_space'))
    error('SPM is not available on the MATLAB path.');
end

if ~isfolder(niftiRoot)
    error('NIfTI root not found: %s', niftiRoot);
end

if ~isfolder(outRoot), mkdir(outRoot); end
if ~isfolder(logRoot), mkdir(logRoot); end

spm('defaults','fmri');
try
    spm_get_defaults('cmdline',true);
catch
end

flags = spm_get_defaults('coreg.estimate');
flags.cost_fun = 'nmi';
flags.sep      = [8 4 2];
flags.fwhm     = [7 7];
flags.graphics = 0;

fprintf('\n============================================================\n');
fprintf('PSE CTP -> T1 COREGISTRATION (CENTERED INITIALISATION)\n');
fprintf('No quantitative CTP map will be resampled.\n');
fprintf('============================================================\n\n');

subjects = dir(fullfile(niftiRoot,'sub-P*'));
subjects = subjects([subjects.isdir]);

if isempty(subjects)
    error('No sub-P* folders found under %s', niftiRoot);
end

[~,ord] = sort({subjects.name});
subjects = subjects(ord);

rows = cell(0,1);
rowN = 0;

%% ------------------------------- LOOP -----------------------------------

for s = 1:numel(subjects)

    subject = subjects(s).name;

    inSession  = fullfile(niftiRoot,subject,sessionName);
    outSession = fullfile(outRoot,subject,sessionName,'CTP');

    fprintf('\n------------------------------------------------------------\n');
    fprintf('%s\n',subject);
    fprintf('------------------------------------------------------------\n');

    t1 = findSingleNifti(fullfile(inSession,'T1'));
    refSrc = findSingleNifti(fullfile(inSession,'CTP_REF'));

    rowN = rowN + 1;
    R = initRow(subject,sessionName);

    if strlength(t1)==0
        R.Status = "ERROR";
        R.Warning = "No unique T1 NIfTI found.";
        rows{rowN} = R;
        printRow(R);
        continue;
    end

    if strlength(refSrc)==0
        R.Status = "MISSING";
        R.Warning = "No CTP_REF NIfTI found.";
        rows{rowN} = R;
        printRow(R);
        continue;
    end

    if ~isfolder(outSession)
        mkdir(outSession);
    end

    refDst = fullfile(outSession,'CTP_REF_coreg_centered.nii');

    try
        copyNiftiAndJson(refSrc,refDst,overwriteExisting);

        % Copy available parametric maps before changing any affine headers.
        mapInfo = struct('modality',{},'src',{},'dst',{}, ...
                         'oldMat',{},'preMat',{},'finalMat',{});

        for k = 1:numel(ctpMaps)

            mod = ctpMaps{k};
            mapSrc = findSingleNifti(fullfile(inSession,mod));

            if strlength(mapSrc)==0
                continue;
            end

            mapDst = fullfile(outSession,mod + "_coreg_centered.nii");
            copyNiftiAndJson(mapSrc,mapDst,overwriteExisting);

            mapInfo(end+1).modality = mod; %#ok<AGROW>
            mapInfo(end).src = mapSrc;
            mapInfo(end).dst = string(mapDst);
            mapInfo(end).oldMat = spm_get_space(char(mapDst));
            mapInfo(end).preMat = [];
            mapInfo(end).finalMat = [];
        end

        % Read current geometry.
        VG = spm_vol(char(t1));
        VF = spm_vol(char(refDst));

        t1Mat  = VG.mat;
        ctpMat = spm_get_space(char(refDst));

        % Geometric centres in world (mm) coordinates.
        t1Centre  = volumeCentreWorld(VG.dim,t1Mat);
        ctpCentre = volumeCentreWorld(VF.dim,ctpMat);

        delta = t1Centre - ctpCentre;

        % Pre-align by translating the CTP volume centre onto the T1 centre.
        Tcentre = eye(4);
        Tcentre(1:3,4) = delta(1:3);

        refPreMat = Tcentre * ctpMat;
        spm_get_space(char(refDst),refPreMat);

        % Apply exactly the same initial translation to every CTP map.
        for k = 1:numel(mapInfo)
            mapInfo(k).preMat = Tcentre * mapInfo(k).oldMat;
            spm_get_space(char(mapInfo(k).dst),mapInfo(k).preMat);
        end

        % Run rigid NMI registration after geometric centre alignment.
        VG = spm_vol(char(t1));
        VF = spm_vol(char(refDst));

        x = spm_coreg(VG,VF,flags);
        M = spm_matrix(x);

        refFinalMat = M \ refPreMat;
        spm_get_space(char(refDst),refFinalMat);

        for k = 1:numel(mapInfo)
            mapInfo(k).finalMat = M \ mapInfo(k).preMat;
            spm_get_space(char(mapInfo(k).dst),mapInfo(k).finalMat);
        end

        R.Status = "OK";
        R.Output = string(refDst);

        R.PreShiftX_mm = delta(1);
        R.PreShiftY_mm = delta(2);
        R.PreShiftZ_mm = delta(3);

        R.Tx_mm = x(1);
        R.Ty_mm = x(2);
        R.Tz_mm = x(3);
        R.Rx_deg = rad2deg(x(4));
        R.Ry_deg = rad2deg(x(5));
        R.Rz_deg = rad2deg(x(6));

        warningText = transformMagnitudeWarning(x);

        if strlength(warningText)>0
            R.Status = "WARNING";
            R.Warning = warningText;
        end

        transform = struct( ...
            'subject',string(subject), ...
            'reference',string(t1), ...
            'source',string(refSrc), ...
            'output',string(refDst), ...
            't1CentreWorld',t1Centre, ...
            'ctpCentreWorld',ctpCentre, ...
            'centreTranslation',delta, ...
            'Tcentre',Tcentre, ...
            'x_spm_coreg',x, ...
            'M_spm_matrix',M, ...
            'ctpOriginalMat',ctpMat, ...
            'ctpPreMat',refPreMat, ...
            'ctpFinalMat',refFinalMat, ...
            'maps',mapInfo);

        save(fullfile(outSession,'CTP_to_T1_centered_transform.mat'), ...
            '-struct','transform');

    catch ME
        R.Status = "ERROR";
        R.Warning = string(ME.message);
    end

    rows{rowN} = R;
    printRow(R);
end

%% ------------------------------- LOG ------------------------------------

Log = struct2table(vertcat(rows{:}));

csvPath = fullfile(logRoot,'PSE_CTP_COREG_CENTERED_LOG.csv');
matPath = fullfile(logRoot,'PSE_CTP_COREG_CENTERED_LOG.mat');

writetable(Log,csvPath);
save(matPath,'Log');

fprintf('\n============================================================\n');
fprintf('CTP CENTERED COREGISTRATION COMPLETE\n');
fprintf('Output: %s\n',outRoot);
fprintf('Log   : %s\n',csvPath);
fprintf('============================================================\n\n');

disp(Log(:,{'Subject','Status', ...
    'PreShiftX_mm','PreShiftY_mm','PreShiftZ_mm', ...
    'Tx_mm','Ty_mm','Tz_mm','Rx_deg','Ry_deg','Rz_deg','Warning'}));

%% ============================== FUNCTIONS ===============================

function c = volumeCentreWorld(dim,mat)
    v = [(double(dim(1))+1)/2; ...
         (double(dim(2))+1)/2; ...
         (double(dim(3))+1)/2; ...
         1];
    w = mat * v;
    c = w(1:3);
end

function R = initRow(subject,sessionName)
    R = struct( ...
        'Subject',string(subject), ...
        'Session',string(sessionName), ...
        'Status',"", ...
        'PreShiftX_mm',NaN, ...
        'PreShiftY_mm',NaN, ...
        'PreShiftZ_mm',NaN, ...
        'Tx_mm',NaN, ...
        'Ty_mm',NaN, ...
        'Tz_mm',NaN, ...
        'Rx_deg',NaN, ...
        'Ry_deg',NaN, ...
        'Rz_deg',NaN, ...
        'Output',"", ...
        'Warning',"");
end

function warningText = transformMagnitudeWarning(x)

    translation = norm(x(1:3));
    rotationDeg = norm(rad2deg(x(4:6)));

    warnings = strings(0,1);

    if translation > 50
        warnings(end+1,1) = sprintf( ...
            'Large residual translation after centre alignment: %.1f mm.',translation); %#ok<AGROW>
    end

    if rotationDeg > 35
        warnings(end+1,1) = sprintf( ...
            'Large residual rotation after centre alignment: %.1f deg.',rotationDeg); %#ok<AGROW>
    end

    if isempty(warnings)
        warningText = "";
    else
        warningText = strjoin(warnings,' | ');
    end
end

function nii = findSingleNifti(folder)

    nii = "";

    if ~isfolder(folder)
        return;
    end

    d = dir(fullfile(folder,'*.nii'));

    if numel(d)~=1
        return;
    end

    nii = string(fullfile(d(1).folder,d(1).name));
end

function copyNiftiAndJson(src,dst,overwriteExisting)

    src = string(src);
    dst = string(dst);

    if isfile(dst)
        if overwriteExisting
            delete(dst);
        else
            error('Output already exists: %s',dst);
        end
    end

    copyfile(src,dst);

    [srcFolder,srcBase,~] = fileparts(src);
    [dstFolder,dstBase,~] = fileparts(dst);

    srcJson = fullfile(srcFolder,srcBase + ".json");
    dstJson = fullfile(dstFolder,dstBase + ".json");

    if isfile(srcJson)
        if isfile(dstJson) && overwriteExisting
            delete(dstJson);
        end
        copyfile(srcJson,dstJson);
    end
end

function printRow(R)

    fprintf('CTP_REF : %-8s',R.Status);

    if isfinite(R.PreShiftX_mm)
        fprintf(' | pre-shift=[%.1f %.1f %.1f] mm', ...
            R.PreShiftX_mm,R.PreShiftY_mm,R.PreShiftZ_mm);
    end

    if isfinite(R.Tx_mm)
        fprintf(' | residual T=[%.1f %.1f %.1f] mm', ...
            R.Tx_mm,R.Ty_mm,R.Tz_mm);
        fprintf(' | R=[%.1f %.1f %.1f] deg', ...
            R.Rx_deg,R.Ry_deg,R.Rz_deg);
    end

    if strlength(R.Warning)>0
        fprintf(' | %s',R.Warning);
    end

    fprintf('\n');
end
