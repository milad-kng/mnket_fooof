function mnket_setup_paths
%MNKET_SETUP_PATHS Sets up all Matlab paths needed for the MNKET analysis. Either takes the
%environmental variable 'MNKETROOT' as project path, or, if that is not defined on the current
%machine, simply takes the path to this m-file as project path.

%% remove all other toolboxes 
restoredefaultpath; 

%% get MNKET project code
% if strcmp(getenv('MNKETROOT'), '')
%     warning('MNKETROOT not defined');
%     pathProject = fullfile(fileparts(mfilename('fullpath')), 'code');
%     % add project path with all sub-paths 
%     addpath(genpath(pathProject));
% else
%    if ~isdir(getenv('MNKETROOT'))
%        warning('MNKETROOT is not a directory');
%        pathProject = fullfile(fileparts(mfilename('fullpath')), 'code'); 
%        % add project path with all sub-paths 
%        addpath(genpath(pathProject));
%    else
%        pathProject = fullfile(getenv('MNKETROOT'));
%        % add project path with all sub-paths 
%        addpath(genpath(pathProject));
%    end
% end
pathProject = fileparts(mfilename('fullpath'));

% remove all other toolboxes
restoredefaultpath;

% to remove obvious path warnings for now
warning off; 

% add project path with all sub-paths
addpath(genpath(pathProject));


%% remove SPM subfolder paths 
% NOTE: NEVER add SPM with subfolders to your path, since it creates
% conflicts with Matlab core functions, e.g., uint16

pathSpm = fileparts(which('spm'));
% remove subfolders of SPM, since it is recommended,
% and fieldtrip creates conflicts with Matlab functions otherwise
rmpath(genpath(pathSpm));
addpath(pathSpm);

modality = 'EEG';

% init batch editor
spm_jobman('initcfg');

% to initialize field trip etc.
spm('defaults', modality);

% large mat files (>2GB) can only be saved with 7.3
spm_get_defaults('mat.format', '-v7.3');

% this should speed up GLM estimation:
spm_get_defaults('stats.maxmem', 2^31); % 2 GB RAM
spm_get_defaults('stats.resmem', true); % store GLM temp files in memory, not on disk

% change to true for cluster job where no graphics can be output
spm_get_defaults('cmdline', false);%% remove SPM subfolder paths 


end