function D = mnket_preprocessing_reject(id, options)
%MNKET_PREPROCESSING_REJECT Preprocesses one subject from the MNKET study using EB rejection.
%   IN:     id          - subject identifier, e.g '0001'
%           optionally:
%           options     - the struct that holds all analysis options
%   OUT:    D           - preprocessed data set

% general analysis options
if nargin < 2
    options = mn_set_analysis_options;
end

options.preproc.eyeblinktreatment = 'reject';

% paths and files
[details, paths] = mn_subjects(id, options);

% record what we're doing
diary(details.logfile);
tnueeg_display_analysis_step_header('preprocessing_reject', 'mnket', id, options.preproc);

% check destination folder
if ~exist(details.preproot, 'dir')
    mkdir(details.preproot);
end
cd(details.preproot);

try
    % check for previous preprocessing
    D = spm_eeg_load(details.prepfile);
    disp(['Subject ' id ' has been preprocessed before.']);
    if options.preproc.overwrite
        clear D;
        disp('Overwriting...');
        error('Continue to preprocessing script');
    else
        disp('Nothing is being done.');
    end
catch
    disp(['Preprocessing subject ' id ' ...']);
    
    %-- preparation -------------------------------------------------------------------------------%
    spm('defaults', 'eeg');
    spm_jobman('initcfg');

    % convert bdf
    D = tnueeg_convert(details.rawfile);
    fprintf('\nConversion done.\n\n');
    
    % set channel types (EEG, EOG) 
    if ~exist(paths.channeldef, 'file')
        chandef = mnket_channel_definition(paths);
    else
        load(paths.channeldef);
    end
    D = tnueeg_set_channeltypes(D, chandef, options);
    fprintf('\nChanneltypes done.\n\n');

    % This is where we swap channels F1 and C2 for the affected subjects, if indicated in the
    % options to do so. 
    if strcmp(options.preproc.swapchannels, 'swap')
        % for some subjects, swap channels
        if details.swapchannels
            montage = mnket_montage_swap_channels;
                D = tnueeg_do_invisible_montage(D, montage, options);
                fprintf('\nChannel swapping done.\n\n');
        end
    end

    % do montage (rereference to average reference, and create EOG channels)
    if ~exist(paths.montage, 'file')
        montage = mnket_montage(options);
        save(paths.montage, 'montage');
    else
        load(paths.montage);
    end
    D = tnueeg_do_montage(D, montage, options);
    fprintf('\nMontage done.\n\n');
    
    %-- filtering ---------------------------------------------------------------------------------%
    D = tnueeg_filter(D, 'high', options);
    if options.preproc.lowpassfreq > 48
        D = tnueeg_filter(D, 'stop', options); 
    end
    D = tnueeg_downsample(D, options);
    D = tnueeg_filter(D, 'low', options); 
    fprintf('\nFilters & Downsampling done.\n\n');       
    
    %-- trigger event correction ------------------------------------------------------------------%
    % manage special case of MMN_4458 during ketamine by removing the first trigger event (for
    % details, see readme_special_case2_4458)
    switch id
        case '4458'
            switch options.condition
                case 'ketamine'
                    D = mnket_remove_incorrect_first_trial(D);
            end
    end
        
    %-- eye blink detection -----------------------------------------------------------------------%
    % subject-specific EB detection thresholds (are always applied to both sessions of the same 
    % participant to ensure compatibility).
    options.preproc.eyeblinkthreshold = details.eyeblinkthreshold;
    % subject-and-session-specific first events of interest
    options.preproc.windowForEyeblinkdetection = details.windowForEyeblinkdetection;
    
    [Dm, trialStats.numEyeblinks] = tnueeg_eyeblink_detection_spm(D, options);
    savefig(details.ebdetectfig);
    fprintf('\nEye blink detection done.\n\n');

    %--trial definition --------------------------------------------------------------------------%
    if ~exist(paths.trialdef, 'file')
        trialdef = mnket_trial_definition(options, paths);
    else
        load(paths.trialdef);
    end
    fprintf('\nTrial definition present.\n\n');
    
    %-- eye blink rejection -----------------------------------------------------------------------%
    [Dc, ~, trialStats.numEyeartefacts, trialStats.idxEyeartefacts, trialStats.nTrialsNoBlinks, ...
        fh, nTrialsInitial] = ...
        tnueeg_eyeblink_rejection_on_continuous_eeg(Dm, trialdef, options);
    savefig(fh, details.eboverlapfig)
    trialStats.nTrialsInitial = nTrialsInitial.tone;
    fprintf('\nEye blink rejection done.\n\n');
    
    %-- experimental epoching ---------------------------------------------------------------------%
    D = tnueeg_epoch_experimental(Dc, trialdef, options);
    fprintf('\nExperimental epoching done.\n\n');
    
    %-- reject remaining artefacts ----------------------------------------------------------------%
    if strcmp(options.preproc.swapchannels, 'reject')
        % for some subjects, reject swapped channels
        if details.swapchannels
            D = badchannels(D, [4, 49], [1 1]);
        end
    end
    [D, trialStats.numArtefacts, trialStats.idxArtefacts, trialStats.badChannels] = ...
        tnueeg_artefact_rejection_threshold(D, options);
    fprintf('\nArtefact rejection done.\n\n');
    
    %-- finish ------------------------------------------------------------------------------------%
    D = move(D, details.prepfilename);
    trialStats.nTrialsFinal = tnueeg_count_good_trials(D);
    save(details.trialstats, 'trialStats');

    disp('   ');
    disp(['Detected ' num2str(trialStats.numEyeblinks) ' eye blinks for subject ' id]);
    disp(['Excluded ' num2str(trialStats.numEyeartefacts.all) ' trials due to eye blinks.']);
    disp(['Rejected ' num2str(trialStats.numArtefacts) ' additional bad trials.']);
    disp(['Marked ' num2str(trialStats.badChannels.numBadChannels) ' channels as bad.']);
    disp([num2str(trialStats.nTrialsFinal.all) ' remaining good trials in D.']);
    fprintf('\nPreprocessing done: subject %s in condition %s.\n', id, options.condition);
    disp('   ');
    disp('*----------------------------------------------------*');
    disp('   ');
   
    

%------tone sequence extraction---------------------------------------%
%      mnCHR_binary_trialDef(D,details);

cd(options.workdir);
close all

diary OFF
end
