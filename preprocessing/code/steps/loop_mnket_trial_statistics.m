function loop_mnket_trial_statistics
%LOOP_MNKET_TRIAL_STATISTICS Loops over conditions and preprocessing options for a summary of trial
%statistics under the given settings.
%   IN:     -
%   OUT:    -

options = mnket_set_analysis_options;

options.preproc.eyeblinktreatment = 'reject';
mnket_trial_statistics(options);

options.preproc.eyeblinktreatment = 'ssp';
mnket_trial_statistics(options);

end

