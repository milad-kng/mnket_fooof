function options = mn_set_analysis_options
%MNKET_SET_ANALYSIS_OPTIONS Analysis options function for MNKET project
%   IN:     -
%   OUT:    options     - a struct with all analysis parameters for the
%                       different steps of the analysis.

%-- where to find the data -----------------------------------------------%
[~, uid] = unix('whoami'); 
switch uid(1: end-1)     
    case 'laptop-b8rt1ssp\milad'
        options.maindir    = 'D:\Science\Coding\mnket_fooof_paper';
        options.preprocdir = fullfile(options.maindir, 'results'); 
    otherwise
        error(['Undefined user. Please specify a user in mn_set_analysis_options ' ...
            'and provide the path to the data']);
end
options.analysis = 'MNPSI' ;%'MNKET';% 'MNPSI', 'group_analysis' 
if strcmp(options.analysis, 'MNKET')
    options.workdir = fullfile(options.preprocdir,'test_mnket'); % 'test_mnpsi'
else
    if strcmp(options.analysis, 'MNPSI')
        options.workdir = fullfile(options.preprocdir,'test_mnpsi'); % 'test_mnpsi'
    end
    if strcmp(options.analysis, 'group_analysis')
        options.workdir = fullfile(options.preprocdir,'group_analysis'); % 'test_mnpsi'
    end
end
options.rawdir  = 'E:\Educational\5.Toronto\Data\mn_test\raw';
options.codedir = fullfile(options.maindir, 'mnket_fooof', 'preprocessing', 'code');
%% Specify default option functions --------------------------------------%
options.funs.details = @mn_subjects; % Specify paths
options.funs.subjects = @mn_set_subject_groups; % Specify subjects groups 
options.funs.eeg = @mn_prepare_eeg; % Specify eeg options
%% Evaluate option functions
 options= feval(options.funs.subjects, options); % Get subject list
 options= feval(options.funs.eeg, options);

%-- condition info -------------------------------------------------------% 
options.condition   = 'placebo'; % 'placebo', 'ketamine','psilocybin','drugdiff'
if strcmp(options.analysis, 'MNKET')
    options.conditions  = {'placebo','ketamine'};
else
    if strcmp(options.analysis, 'MNPSI')
            options.conditions  = {'placebo','psilocybin'};
    end
end
%-- preparation ----------------------------------------------------------%
options.prepare.subjectIDs  = options.subjects.all; % data preparation (tone sequences)
options.prepare.overwrite   = 0; % whether to overwrite any previous prep
                           
%-- modeling -------------------------------------------------------------%
options.model.subjectIDs    = options.subjects.all; % modeling with the HGF
options.model.overwrite     = 0; % whether to overwrite any previous model

%-- preprocessing --------------------------------------------------------%
options.preproc.subjectIDs      = options.subjects.all;
options.preproc.overwrite       = 0; % whether to overwrite any prev. prepr
options.preproc.keep            = 1;  % whether to keep intermediate data

% swap channel option: this is where we decide what to do about the data sets with apparently
% swapped electrodes C2 and F1 (see readme and swapchannels.pdf for details)
options.preproc.swapchannels    = 'reject'; % 'swap'(swap channels back), 'reject'(mark as bad), 
                                            % '' (do nothing about it)

options.preproc.rereferencing   = 'average';
options.preproc.keepotherchannels = 1;
options.preproc.lowpassfreq     = 45;
options.preproc.highpassfreq    = 0.5; 
options.preproc.downsamplefreq  = 256;

options.preproc.trialdef            = 'tone'; % 'MMN', 'tone' - choose 'tone' for modelbased analysis
options.preproc.epochwin            = [-100 400];
options.preproc.baselinecorrection  = 1;

options.preproc.eyeblinktreatment   = 'reject'; % 'reject', 'ssp'
options.preproc.mrifile             = 'template';
options.preproc.eyeblinkchannels    = {'VEOG'};
options.preproc.windowForEyeblinkdetection = 3; % first event of interest (and optionally last)
% NOTE: This sets the default index of the first even of interest in the EEG file, however, this 
% will be adjusted individually for subjects if their EEG file requires a different value. For all
% adjustments, see mnket_subjects.
options.preproc.eyeblinkthreshold   = 3; % for SD thresholding: in standard deviations, for amp in uV
% NOTE: This sets the default threshold for detecting eye blink events in the EOG, however, this 
% will be adjusted individually for subjects if their EOG requires a different threshold. For all
% adjustments, see mnket_subjects.
options.preproc.eyeconfoundcomps    = 1;
options.preproc.eyeblinkmode        = 'eventbased'; % uses EEG triggers for trial onsets
options.preproc.eyeblinkwindow      = 0.5; % in s around blink events
options.preproc.eyeblinktrialoffset = 0.1; % in s: EBs won't hurt <100ms after tone onset
options.preproc.eyeblinkEOGchannel  = 'VEOG'; % EOG channel (name/idx) to plot
options.preproc.eyebadchanthresh    = 0.4; % prop of bad trials due to EBs
options.preproc.eyecorrectionchans  = {'Fp1', 'Fz', 'AF8', 'T7', 'Oz'};

% in case we use the SSP eye blink correction, this section defines the amount of pre-cleaning
options.preproc.preclean.doFilter           = true;
options.preproc.preclean.lowPassFilterFreq  = 10;
options.preproc.preclean.doBadChannels      = true;
options.preproc.preclean.doRejection        = true;
options.preproc.preclean.badtrialthresh     = 500;
options.preproc.preclean.badchanthresh      = 0.5;
options.preproc.preclean.rejectPrefix       = 'cleaned_';

options.preproc.badchanthresh   = 0.2; % proportion of bad trials
options.preproc.badtrialthresh  = 80; % in microVolt

%-- erp ------------------------------------------------------------------%
options.erp.subjectIDs  = options.subjects.all;                        
options.erp.overwrite   = 0; % whether to overwrite any previous erp

options.erp.type        = 'roving';%'roving';  % roving (sta=6, dev>=5), mmnad (sta=6, dev=1), 
                            % tone (nothing), memory (t6, t8, t10), 
                            % repetitions (t1, t2, ..., t11)
options.erp.electrode   = 'C1';
options.erp.averaging   = 's'; % s (standard), r (robust)
switch options.erp.averaging
    case 'r'
        options.erp.addfilter = 'f';
    case 's'
        options.erp.addfilter = '';
end
options.erp.percentPE = 20; % for ERP type lowhighPE: how many percent PEs

options.erp.contrastWeighting   = 1;
options.erp.contrastPrefix      = 'diff_';
options.erp.contrastName        = 'mmn';
% these channels are set after inspecting the main results of the modelbased statistics: we plot all
% channels' grand averages + 1.96 * standard error, which lie in the peak of significant clusters of
% 2nd level contrasts
options.erp.channels            = {'C4','C3', 'C1', 'C2','Cz', ...
                                    'FC1', 'FC2', 'FCz', ...
                                    'F1', 'F2', 'Fz', ...
                                    'P7', 'P8', 'P9', 'P10', ...
                                    'TP7'};


%-- T-F analysis ---------------------------------------------------------%
options.tf.mode = 'trial_series_csd'; %options: 'whole signal', 'trial', 'trial series csd'
options.tf.spectrum.overwrite = 1;
options.tf.nTrials = 'all';
options.tf.subtract_erp = 1;              
end

