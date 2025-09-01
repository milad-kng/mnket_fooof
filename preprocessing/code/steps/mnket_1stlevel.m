function mnket_1stlevel(id, options)
%MNKET_1STLEVEL Computes the first level contrast images for modelbased 
%statistics for one subject from the MNKET study.
%   IN:     id      - subject identifier, e.g '0001'
%           options - the struct that holds all analysis options
%   OUT:    --

% general analysis options
if nargin < 2
    options = mn_set_analysis_options;
end

% paths and files
[details, ~] = mn_subjects(id, options);

% record what we're doing
diary(details.logfile);
mnket_display_analysis_step_header('firstlevel stats', id, options.stats);

try
    % check for previous statistics
    load(details.spmfile);
    disp(['1st level stats for subject ' id ' have been computed before.']);
    if options.stats.overwrite
        delete(details.spmfile);
        disp('Overwriting...');
        error('Continue to 1st level stats step');
    else
        disp('Nothing is being done.');
    end
catch
    disp(['Computing 1st level stats for subject ' id ' ...']);
    
    % first, correct the regressors for trials that were rejected during
    % preprocessing
    if strcmp(id, '4534') && strcmp(options.condition, 'ketamine') 
        % special case for subject 4534 in ketamine condition
        design = mnket_correct_regressors_for_EEG_artefacts_special_case(details, options);
    else
        design = mnket_correct_regressors_for_EEG_artefacts(details, options);
    end
    save(details.eegdesign, 'design');
    disp(['Corrected regressors for subject ' id]);
    
    % make sure we have a results directory
    facdir = details.statroot;
    if ~exist(facdir, 'dir')
        mkdir(facdir);
    end
    
    % smoothed images of final preprocessed EEG file serve as input for the
    % statistics - one image per trial, but here, we only indicate the path
    % to and basename of the images
    impath = details.smoofile{1};
    
    % compute the regression of the EEG signal onto our single-trial
    % regressors
    tnueeg_singletrial_1stlevel_stats(facdir, impath, design);
    disp(['Computed statistics for subject ' id]);
end
cd(options.workdir);

diary OFF
end



