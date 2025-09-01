function mnket_quality_check( options )
%MNKET_QUALITY_CHECK Performs all quality checks for the MNKET study.
%   IN:     optionally:
%           options         - the struct that contains all analysis options
%   OUT:    -

if nargin < 1
    options = mn_set_analysis_options;
end

[~, paths] = mn_subjects(options);


if ~exist(paths.qualityroot, 'dir')
    mkdir(paths.qualityroot);
end

% Modelling
% no checks implemented so far...

% Collect trial numbers
mnket_trial_statistics(options);

% Eye-blink treatment
mnket_summarize_step('eyeblinkdetection', options);
switch lower(options.preproc.eyeblinktreatment)
    case 'ssp'
        mnket_summarize_step('eyeblinkconfounds', options);
        mnket_report_step_posthoc('eyeblinkcorrection', options);
        mnket_summarize_step('eyeblinkcorrection', options);
    case 'berg'
        mnket_summarize_step('eyeblinkconfounds', options);
        mnket_report_step_posthoc('eyeblinkcorrection', options);
        mnket_summarize_step('eyeblinkcorrection', options);
        mnket_report_step_posthoc('coregistration', options);
        mnket_summarize_step('coregistration', options);
    case 'reject'
        mnket_summarize_step('eyeblinkoverlap', options);
end

% First level stats
mnket_report_step_posthoc('regressors', options);
mnket_summarize_step('regressors', options);
mnket_report_step_posthoc('firstlevelmask', options);
mnket_summarize_step('firstlevelmask', options);
%mnket_summarize_step('firstlevelregressors', options); %TODO
%mnket_summarize_step('firstleveldesign', options); %TODO




end
