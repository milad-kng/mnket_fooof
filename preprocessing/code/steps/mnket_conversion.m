function mnket_conversion(id, options)
%MNKET_CONVERSION Converts preprocessed epoched EEG data to 3D images and 
%smoothes them for one subject from the MNKET study.
%   IN:     id                  - subject identifier, e.g '0001'
%   OUT:    --

% general analysis options
if nargin < 2
    options = mn_set_analysis_options;
end

% paths and files
[details, ~] = mn_subjects(id, options);

% record what we're doing
diary(details.logfile);
mnket_display_analysis_step_header('conversion', id, options.conversion);

% determine the preprocessed file to use for conversion
switch options.conversion.mode
    case 'diffWaves'
        prepfile = details.difffile;
    case 'modelbased'
        prepfile = details.prepfile_modelbased;
    case 'ERPs'
        prepfile = details.erpfile;
    case 'mERPs'
        prepfile = details.mergfile;
end

try
    % check for previous smoothing
    im = spm_vol(details.smoofile{1});
    disp(['Images for subject ' id ' have been converted and smoothed before.']);
    if options.conversion.overwrite
        clear im;
        disp('Overwriting...');
        error('Continue to conversion step');
    else
        disp('Nothing is being done.');
    end
catch
    disp(['Converting subject ' id ' ...']);
    
    switch options.conversion.mode
        case 'modelbased'
            % in the modelbased analysis, we have to remove the first EEG trial because the model 
            % only defines PEs starting from the 2nd tone (1st observed transition)            

            if strcmp(id, '4534') && strcmp(options.condition, 'ketamine')
                % special case for subject 4534 in ketamine condition: first 29 trials were not 
                % recorded in the EEG, therefore we don't need to correct the EEG for the 1st trial. 
                D = spm_eeg_load(details.prepfile);
                D = copy(D, details.prepfile_modelbased);
%                 clear D;
                prepfile = details.prepfile_modelbased;
            else
                prepfile = mnket_remove_first_EEG_trial(details, options);
            end
    end
    
    % try whether this file exists
    try
        D = spm_eeg_load(prepfile);
    catch
        disp('This file: ')
        disp(prepfile)
        disp('could not been found.')
        error('No final preprocessed EEG file')
    end

    % convert EEG data
    [images, ~] = tnueeg_convert2images(D, options);
    disp(['Converted EEG data for subject ' id ' in mode ' options.conversion.mode]);

    % and smooth the resulting images
    tnueeg_smooth_images(images, options);
    disp(['Smoothed images for subject ' id ' in mode ' options.conversion.mode]);
end

cd(options.workdir);

diary OFF
end
