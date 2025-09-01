function mnket_analyze_subject( id, options )
%MNKET_ANALYZE_SUBJECT Performs all analysis steps for one subject of the MNKET study (up until 
%first level modelbased statistics)
%   IN:     id  - subject identifier string, e.g. '0001'
%   OUT:    --

if nargin < 2
    options = mn_set_analysis_options;
end

% loop over both sessions for this subject
for optionsCell = options.conditions
    options.condition = char(optionsCell);
    fprintf('\n\n --------- Working on: %s session ---------\n\n', upper(options.condition));
    
    %% Preparation and Modeling
      mnket_data_preparation(id, options);
    
    %% Pre-processing: reject eyeblinks
    options.preproc.eyeblinktreatment = 'reject';
    fprintf('\n\n --- Subject analysis using: %s method ---\n\n', upper(options.preproc.eyeblinktreatment));
    mnket_preprocessing_reject(id, options);
     
    %% ERP analysis: tone definition

    options.erp.type = 'tone';
    mnket_erp(id, options);

    options.erp.type = 'roving';
    mnket_erp(id, options);

    mnket_time_freq_analysis(id, options)
 
end
    
    
    
  



