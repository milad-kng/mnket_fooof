function loop_mnket_subject_analysis( options )
%LOOP_MNKET_SUBJECT_ANALYSIS Loops over all subjects in the MNKET study and executes all analysis 
%steps
%   IN:     options     - the struct that holds all analysis options
%   OUT:    -

if nargin < 1 
    options = mn_set_analysis_options;
end


for idCell = options.subjects.all
    id = char(idCell);

    mnket_analyze_subject(id, options);

end

