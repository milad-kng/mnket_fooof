% Matlab startup file for MNKET

%% remove all other toolboxes 
restoredefaultpath; 

%% get MNKET project code
if strcmp(getenv('MNKETROOT'), '')
    warning('MNKETROOT not defined');
    pathProject = fullfile(fileparts(mfilename('fullpath')), 'code');
    % add project path with all sub-paths 
    addpath(genpath(pathProject));

else
   if ~isdir(getenv('MNKETROOT'))
       warning('MNKETROOT is not a directory');
       pathProject = fullfile(fileparts(mfilename('fullpath')), 'code'); 
       % add project path with all sub-paths 
       addpath(genpath(pathProject));
   else
       pathProject = fullfile(getenv('MNKETROOT'));
       addpath(genpath(pathProject));
   end
end

%% remove SPM subfolder paths 
% NOTE: NEVER add SPM with subfolders to your path, since it creates 
% conflicts with Matlab core functions, e.g., uint16 
pathSpm = fileparts(which('spm')); 
% remove subfolders of SPM, since it is recommended, 
% and fieldtrip creates conflicts with Matlab functions otherwise 
rmpath(genpath(pathSpm)); 
addpath(pathSpm);