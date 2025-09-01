% --- Analysis script for MnketAnalysis EEG dataset --- %
% Adapted from mnketAnalysis repo by Dr. Lilian A. Weber and Dr. Andreea O. Diaconescu
% Author: Gabrielle Allohverdi

%% Setup
% % set up Matlab environment
mnket_setup_paths;

% set all analysis options and provide the path to the data
options = mn_set_analysis_options;

% create the folder structure needed for the full analysis, and fill with necessary raw data
mnket_setup_analysis_folder(options);

%% run the full first-level analysis
% includes: data preparation, EEG preprocessing, ERPs, conversion to images, 1st level statistics
loop_mnket_subject_analysis(options); 
%% summarize quality of preprocessing and trial statistics
loop_mnket_quality_check(options);

%% Group level t-f analysis
options.tf.spectrum.overwrite = 0;
options.tf.collectMMN = 1;
options.datatype = 'psi';
close all
options.erp.type = 'roving';%'roving';
mnket_tf_group(options)
close all



