function [ details, paths ] = mn_subjects( id, options )
%MNKET_SUBJECTS Function that sets all filenames and paths
%   IN:     EITHER (for quering both general and subject-specific paths:
%           id                  - the subject number as a string, e.g. '0001'
%           options (optional)  - the struct that holds all analysis options
%           OR (for only quering general paths & files):
%           options - the struct that holds all analysis options 
%   OUT:    details     - a struct that holds all filenames and paths for
%                       the subject with subject number id
%           paths       - a struct that holds all paths and config files
%                       that are not subject-specific

%-- check input -----------------------------------------------------------------------------------%
if isstruct(id)
    options = id;
    id = 'dummy';
elseif ischar(id) && nargin < 2
    options = mn_set_analysis_options;
end

%-- general paths and files -----------------------------------------------------------------------%
paths.confroot      = fullfile(options.workdir, 'config');
paths.tonesroot     = fullfile(options.workdir, 'tones');
paths.erproot       = fullfile(options.workdir, 'erp');
paths.statroot      = fullfile(options.workdir, 'stats_model');
paths.erpstatroot   = fullfile(options.workdir, 'stats_erp');
paths.dcmroot   =   fullfile(options.workdir, 'dcm', options.condition);
paths.grouproot =  fullfile(options.workdir, 'groupdiff');
paths.tf        =  fullfile(options.workdir, 'tf');
paths.qualityroot = fullfile(options.workdir, 'quality');
% config files
paths.paradigm      = fullfile(paths.tonesroot, 'paradigm.mat');
paths.elec          = fullfile(paths.confroot, 'mnket_sensors.sfp');
paths.channeldef    = fullfile(paths.confroot, 'mnket_chandef.mat');
paths.montage       = fullfile(paths.confroot, 'mnket_avref_eog.mat');
paths.trialdef      = fullfile(paths.confroot, ['mnket_trialdef_' options.preproc.trialdef '.mat']);

% erp analysis folders
paths.erpfold       = fullfile(paths.erproot, options.erp.type, options.condition);
paths.erpdiffold    = fullfile(paths.erproot, options.erp.type, 'drugdiff');
paths.erpgafold     = fullfile(paths.erpfold, 'GA');
paths.erpspmfold    = fullfile(paths.erpfold, 'SPM');

% erp analysis files
for iCh = 1: numel(options.erp.channels)
    paths.erpgafiles{iCh} = fullfile(paths.erpgafold,['tnu_ga_' options.condition ...
        '_' options.erp.type '_' options.erp.channels{iCh} '.mat']);
    paths.diffgafiles{iCh} = fullfile(paths.erpdiffold, ['tnu_ga_diffwaves_' options.erp.type ...
        '_' options.erp.channels{iCh} '.mat']);                
end
paths.spmganame     = ['spm_GA_' options.condition '_' options.erp.type '.mat'];  
paths.spmgafile     = fullfile(paths.erpspmfold, paths.spmganame);    
paths.dcmfile   = fullfile(paths.dcmroot, ['dcm_' options.erp.type '.mat']);

% model stats folders
paths.statfold      = fullfile(paths.statroot, options.condition);
paths.statdifffold  = fullfile(paths.statroot, 'drugdiff');

% Groupxcond stats folders
%paths.groupfold = fullfile(paths.grouproot, options.condition);
paths.groupfold = fullfile(paths.grouproot);
% erp stats folders
paths.erpstatfold      = fullfile(paths.erpstatroot, options.erp.type, options.condition);
paths.erpstatdifffold  = fullfile(paths.erpstatroot, options.erp.type, 'drugdiff');

% erp stats files
paths.erpspmfile    = fullfile(paths.erpstatfold, 'SPM.mat');

% tf group files
paths.psd_std_pla = fullfile(paths.tf, ['psd_std_pla_' options.erp.type '_'...
                                        options.tf.mode '_' ...
                                        options.tf.nTrials '.mat']);
paths.psd_std_ket = fullfile(paths.tf, ['psd_std_ket_' options.erp.type '_',...
                                        options.tf.mode '_' ...
                                        options.tf.nTrials '.mat']);
paths.psd_dev_pla = fullfile(paths.tf, ['psd_dev_pla_' options.erp.type '_',...
                                        options.tf.mode '_' ...
                                        options.tf.nTrials '.mat']);
paths.psd_dev_ket = fullfile(paths.tf, ['psd_dev_ket_' options.erp.type '_',...
                                        options.tf.mode '_' ...
                                        options.tf.nTrials '.mat']);

paths.psd_std_pla2 = fullfile(paths.tf, ['psd_std_pla2_' options.erp.type '_'...
                                        options.tf.mode '_' ...
                                        options.tf.nTrials '.mat']);
paths.psd_std_psi = fullfile(paths.tf, ['psd_std_psi_' options.erp.type '_',...
                                        options.tf.mode '_' ...
                                        options.tf.nTrials '.mat']);
paths.psd_dev_pla2 = fullfile(paths.tf, ['psd_dev_pla2_' options.erp.type '_',...
                                        options.tf.mode '_' ...
                                        options.tf.nTrials '.mat']);
paths.psd_dev_psi = fullfile(paths.tf, ['psd_dev_psi_' options.erp.type '_',...
                                        options.tf.mode '_' ...
                                        options.tf.nTrials '.mat']);
paths.subj_list = fullfile(paths.tf, 'subj_list.mat');

% logging 2nd level analyses and quality check
paths.logfile       = fullfile(options.workdir, 'secondlevel.log');
paths.qualityroot   = fullfile(options.workdir, 'quality');
paths.qualityfold   = fullfile(paths.qualityroot, options.preproc.eyeblinktreatment, ...
    options.condition);
paths.trialstatstab = fullfile(paths.qualityfold, [options.condition, '_' ...
    options.preproc.eyeblinktreatment, '_table_trial_stats.mat']);
paths.trialstatsfig = fullfile(paths.qualityfold, [options.condition, '_' ...
    options.preproc.eyeblinktreatment, '_overview_trial_stats']);

paths.dcmFigParam = fullfile(paths.dcmroot, ['dcm_BParameters' options.erp.type '.fig']);
paths.dcmPredictedSimulated...
                    = fullfile(paths.dcmroot, ['dcm_PredictedSimulated' options.erp.type '.fig']);
paths.dcmScalpMaps...
                    = fullfile(paths.dcmroot, ['dcm_ScalpMaps' options.erp.type '.fig']);

paths.trialstats_summary_all = fullfile(paths.qualityroot, options.preproc.eyeblinktreatment, ...
    'summary_trial_stats_across_drugs.mat');
% in case we exclude participant 4497: recomputed trial statistics
paths.trialstats_summary_red = fullfile(paths.qualityroot, options.preproc.eyeblinktreatment, ...
    'summary_trial_stats_across_drugs_without_4497.mat');


%-- subject-specific options, paths and files -----------------------------------------------------%
% swap channels F1 and C2: which datasets are affected
switch options.condition
    case 'placebo'
        switch id
            case {'4520', '4534'}
                details.swapchannels = true;
            otherwise
                details.swapchannels = false;
        end
    case 'ketamine'
        switch id
            case {'4422', '4488', '4520'}
                details.swapchannels = true;
            otherwise
                details.swapchannels = false;
        end

    case 'psilocybin'
        switch id
            case {''}
                details.swapchannels = true;
             otherwise
                details.swapchannels = false;
        end
end

        
% EB detection threshold
switch id
    case {'4497', '4447', '4478'}
        details.eyeblinkthreshold = 2;
    case {'4499', '4532'}
        details.eyeblinkthreshold = 5;   
    otherwise
        details.eyeblinkthreshold = options.preproc.eyeblinkthreshold;
end

% index of first event of interest for EB detection
switch options.condition
    case 'placebo'
        switch id
            case {'4447', '4459', '4487', '4494', '4497', '4500', '4532', '4548'}
                details.windowForEyeblinkdetection = 4;
            otherwise
                details.windowForEyeblinkdetection = 3;
        end
    case 'ketamine'
        switch id
            case {'4446', '4459', '4478', '4482', '4487', '4494', '4499', '4507'}
                details.windowForEyeblinkdetection = 4;
            otherwise
                details.windowForEyeblinkdetection = 3;
        end

    case 'psilocybin'
        switch id
            case {''}
                details.windowForEyeblinkdetection = 4;
            otherwise
                details.windowForEyeblinkdetection = 3;
        end

end

% bad channels before EB confound estimation (only needed for SSP eyeblink correction)
switch options.condition
    case 'placebo'
        switch id
            case '4507'
                details.preclean.badchannels = 39; % F2
            otherwise
                details.preclean.badchannels = [];
        end
    case 'ketamine'
        switch id
            case '4422'
                details.preclean.badchannels = [28 59]; % Iz, P6
            case '4494'
                details.preclean.badchannels = 9; % FC5
            otherwise
                details.preclean.badchannels = [];
        end
end
    
% raw file names
switch options.condition
    case 'placebo'
        if strcmp(options.analysis,'MNKET')
            rawsuffix = '_1_pla';
            if strcmp(id, '4532') || strcmp(id, '4534')
                rawsuffix = '_pla';
            end
        else
            if strcmp(options.analysis,'MNPSI')
                rawsuffix = '_2_pla';
            end
            if strcmp(options.analysis,'group_analysis')
                rawsuffix = '';
            end
        end
    case 'ketamine'
        rawsuffix = '_1_ket';
        if strcmp(id, '4532') || strcmp(id, '4534')
            rawsuffix = '_ket';
        end
    case 'psilocybin'
        rawsuffix = '_2_psi';
end

% names
details.subjectname  = ['MMN_' id];

details.rawfilename  = [details.subjectname rawsuffix];
details.prepfilename = [details.subjectname '_prep'];
details.erpfilename  = [details.subjectname '_' options.erp.type '_erp'];
details.mergfilename = [details.subjectname '_' options.erp.type '_erp_merged'];

% directories
details.subjectroot = fullfile(options.workdir, 'subjects', options.condition, details.subjectname);

details.rawroot     = fullfile(details.subjectroot, 'eeg');  
details.tonesroot   = fullfile(details.subjectroot, 'tones');  
details.preproot    = fullfile(details.subjectroot, 'spm_prep', options.preproc.eyeblinktreatment);
details.erproot     = fullfile(details.subjectroot, 'spm_erp', options.erp.type);
details.tfroot      = fullfile(details.subjectroot, 'tf');
% Added by Milad S.
if ~exist(fullfile(details.subjectroot,'tf'), 'dir')
    mkdir(details.subjectroot,'tf');
end
% files
details.logfile     = fullfile(details.subjectroot, [details.subjectname '.log']);
details.tonestxt    = fullfile(paths.tonesroot, 'textfiles', options.condition, ['sub' id '.txt']);
details.tonesfile   = fullfile(details.tonesroot, 'tones.mat');
details.eegtones    = fullfile(details.tonesroot, 'eegtones.mat');


details.rawfile     = fullfile(details.rawroot, [details.rawfilename '.bdf']);
details.rawfile_edf = fullfile(details.rawroot, [details.rawfilename '.edf']);
details.prepfile    = fullfile(details.preproot, [details.prepfilename '.mat']);
details.prepfile_modelbased = fullfile(details.preproot, [details.prepfilename '_modelbased.mat']);
details.ebfile      = fullfile(details.preproot, ['fEBbfdfMspmeeg_' details.rawfilename '.mat']);

details.trialstats  = fullfile(details.preproot, [details.subjectname '_trialStats.mat']);
details.artefacts   = fullfile(details.preproot, [details.subjectname '_nArtefacts.mat']);
details.eyeblinks   = fullfile(details.preproot, [details.subjectname '_nEyeblinks.mat']);

details.redeffile   = fullfile(details.erproot, ['redef_' details.subjectname '.mat']);
details.avgfile     = fullfile(details.erproot, ['avg_' details.subjectname '.mat']);
details.erpfile     = fullfile(details.erproot, [details.erpfilename '.mat']);
details.difffile    = fullfile(details.erproot, ['diff_' details.erpfilename '.mat']);
details.mergfile    = fullfile(details.erproot, [details.mergfilename]);

% Added by Milad S.:------------
details.prepfile_replaced = fullfile(details.tfroot,[details.subjectname '_prepfile_replaced.mat']);
details.fullStd     = fullfile(details.tfroot,[details.subjectname '_prepfile_stds.mat']);
details.fullDev     = fullfile(details.tfroot,[details.subjectname '_prepfile_devs_' options.tf.mode '.mat']);
details.fullStd_pwr = fullfile(details.tfroot,[details.subjectname '_PSpect_stds_' options.tf.mode '.mat']);
details.fullDev_pwr = fullfile(details.tfroot,[details.subjectname '_PSpect_devs_' options.tf.mode '.mat']);
details.trialPwr    = fullfile(details.tfroot,[options.erp.type details.subjectname '_trialPwr_' options.tf.mode '.mat']);
details.trialTF     = fullfile(details.tfroot,[options.erp.type details.subjectname '_trialTF_' options.tf.mode '.mat']);
details.CSD         = fullfile(details.tfroot,[options.erp.type details.subjectname '_CSD_' options.tf.mode '.mat']);
details.PSD         = fullfile(details.tfroot,[options.erp.type details.subjectname '_PSD_' options.tf.mode '.mat']);

details.waveletCoefDev = fullfile(details.tfroot,[details.subjectname '_waveletCoefDev.mat']);
details.subbandsDev = fullfile(details.tfroot,[details.subjectname '_subbandsDev.mat']);
details.waveletCoefStd = fullfile(details.tfroot,[details.subjectname '_waveletCoefStd.mat']);
details.subbandsStd = fullfile(details.tfroot,[details.subjectname '_subbandsStd.mat']);

details.PSpect      = fullfile(details.tfroot,[details.subjectname '_PSpect_' options.tf.mode '.fig']);
details.denoised    = fullfile(details.tfroot,[details.subjectname '_denoised.fig']);
details.trialPowerFig  = fullfile(details.tfroot,[options.erp.type details.subjectname '_Power_' options.tf.mode '.fig']);
details.trialTFFig     = fullfile(details.tfroot,[options.erp.type details.subjectname '_trialTF_' options.tf.mode '.fig']);
%-------------------------------

% fiducials
details.fid.labels  = {'NAS'; 'LPA'; 'RPA'};
details.fid.data    = [1, 85, -41; -83, -20, -65; 83, -20, -65];

% figures
details.ebdetectfig     = fullfile(details.preproot, [details.subjectname '_eyeblinkDetection.fig']);
details.ebspatialfig    = fullfile(details.preproot, [details.subjectname '_eyeblinkConfounds.fig']);
% only needed for EB rejection:
details.eboverlapfig      = fullfile(details.preproot, [details.subjectname '_eyeblinkTrialOverlap.fig']);
% only needed for EB correction:
details.ebcorrectfig    = fullfile(details.preproot, [details.subjectname '_eyeblinkCorrection.fig']);
details.coregdatafig    = fullfile(details.preproot, [details.subjectname '_coregistration_data.fig']);
details.coregmeshfig    = fullfile(details.preproot, [details.subjectname '_coregistration_mesh.fig']);


details.erpfigure         = fullfile(details.erproot, [details.subjectname '_ERP_' ...
    options.erp.electrode '.fig']);

end
