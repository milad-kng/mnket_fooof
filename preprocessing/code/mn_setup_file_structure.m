function options = mn_setup_file_structure(options)
%--------------------------------------------------------------------------
% Sets up results structure.
% modified from mnCHR_setup_file_structure by COLLEEN CHARLTON
% GABRIELLE ALLOHVERDI
%--------------------------------------------------------------------------
%% Create some roots
% You can simply add new root directories here
options.roots.results               = fullfile(options.workdir, options.analysis);

%% Create folder tree
if ~exist(options.workdir, 'dir')
    mkdir(options.workdir);
end
cd(options.workdir);

folders = struct2cell(options.roots);
for i = 1:size(folders,1)
    if ~(exist(folders{i},'dir'))
        mkdir(folders{i})
    end
end
fprintf('Created analysis folder tree at %s.\n', options.roots.results);

end