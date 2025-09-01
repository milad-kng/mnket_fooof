function design = mnket_model(id, options)
%MNKET_MODEL Simulates the beliefs of one subject from the MNKET study and saves the trajectories
%for modelbased analysis of EEG data.
%   IN:     id          - subject identifier, e.g '0001'
%           optionally:
%           options     - the struct that holds all analysis options
%   OUT:    design      - the design file which holds the modelbased regressors for this subject

% general analysis options
if nargin < 2
    options = mn_set_analysis_options;
end

% paths and files
details = mn_subjects(id, options);

% record what we're doing
diary(details.logfile);
tnueeg_display_analysis_step_header('model', 'mnket', id, options.model);

% check destination folder
if ~exist(details.modelroot, 'dir')
    mkdir(details.modelroot);
end

try
    % check for previous preprocessing
    load(details.design);
    disp(['Subject ' id ' has been modeled before.']);
    if options.model.overwrite
        clear design;
        disp('Overwriting...');
        error('Continue to modeling step.');
    else
        disp('Nothing is being done.');
    end
catch
    disp(['Modeling subject ' id ' ...']);
    
    %-- simulate beliefs --------------------------------------------------------------------------%
    tones = getfield(load(details.tonesfile), 'tones');
    sim = mnket_simulate_beliefs(tones);
    save(details.simfilepre, 'sim');
    fprintf('\nSimulation done.\n\n');
    
    %-- calculate trial-by-trial PEs --------------------------------------------------------------%
    sim = mnket_calculate_regressors(sim);
    % Special case for subject 4421 in the psilocybin condition. The first two trials are were not recorded, so we
    % begin at the 3rd trajectory in the sequence
    switch id
        case {'4421'}
            switch options.condition
                case {'psilocybin'}
                    switch options.erp.type
                        case 'precision'
                            sim.reg.pihat2 = [sim.reg.pihat2(3:end)];
                            sim.reg.pihat3 = [sim.reg.pihat3(3:end)];
                            sim.reg.pihat1 = [sim.reg.pihat1(3:end)];
                        case {'epsilon', 'lowhighEpsi2', 'lowhighEpsi3'}
                            sim.reg.epsi2 = [sim.reg.epsi2(3:end)];
                            sim.reg.epsi3 = [sim.reg.epsi3(3:end)];
                    end
                    
            end
    end
    
                    
    save(details.simfilepost, 'sim');
    fprintf('\nPE calculation done.\n\n');
    
    %-- create design files for EEG analysis ------------------------------------------------------%
    design = mnket_create_design_file(sim, options.stats.design);
    % We have a special case for subject MMN_4534 in condition ketamine: EEG recording started later
    % than audio stimulation, therefore we have to skip the first 28 regressor entries for EEG.
    % For longer explanation, see contents of the textfile readme_4534_special_case.txt.
    switch id
        case {'4534'}
            switch options.condition
                case {'ketamine'}
                    condNames = fieldnames(design);
                    for iCon = 1: numel(condNames)
                        design.(condNames{iCon})(1: 28) = [];
                    end
                    fprintf(...
                        '\nRemoved the first 28 trials from design file for MMN_4534 in ketamine');
            end
    end
    save(details.design, 'design');
    fprintf('\nDesign file has been created.\n\n');
    
    fprintf('\nModeling done: subject %s in condition %s.\n', id, options.condition);
    disp('   ');
    disp('*----------------------------------------------------*');
    disp('   ');
    
end

cd(options.workdir);
close all

diary OFF
end
