# -*- coding: utf-8 -*-
"""
The pipeline to analyze the power spectral density features of ketamine and 
psilocybin datasets obtained using fooof (SpecParam) model.

@author: Milad Soltanzadeh
"""
# =============================================================================
# Importage
# =============================================================================
import numpy as np
import scipy.io as sio
import pickle
# import os
import pandas as pd
import matplotlib.pyplot as plt
# import matplotlib
# from   matplotlib import cm, colors, colorbar
# import matplotlib.patches as patches
import seaborn as sns
# from matplotlib.colors import TwoSlopeNorm

# from sklearn.metrics import r2_score
# from scipy.stats import ttest_ind
import os
from os.path import isfile, exists
from aux_functions import extract_params_fg,check_nans,\
                          get_analysis_options
from fooof_qc import plot_retention_summary, plot_group_barplot,\
                     get_good_channels_by_r2, plot_retention_boxplots,\
                         compute_common_retained_channels, qc_filter,\
                         qc_plot_detectable, filter_common_channels_between_drugs
from fooof_pls_loops import run_fooof, load_fooof, loop_param_extract, run_plot_mpls,\
                            prep_ket_psi_param_diff, run_plot_bpls
from plot_helpers import plot_violin_distributions
from scipy.stats import pearsonr
# from scipy.stats import ttest_ind
# from statsmodels.multivariate.manova import MANOVA

# =============================================================================
# Set global plotting attributes
# =============================================================================
plt.rcParams['axes.edgecolor'] = 'black'
plt.rcParams['axes.linewidth'] = 1.5  
# font = {'family': 'Times New Roman',
#         'weight': 'normal'}  

# # Apply font properties
# plt.rc('font', **font)
# sns.set(font='Times New Roman') 

# =============================================================================
# Setup Paths
# =============================================================================
results_path = 'results'
file_path = results_path + '/files/'
fig_path = results_path + '/figs/'
# data_path = 'psd_data/fourtyFiveHzWithoutERP/Newer/'
data_path = 'psd_data/working_files/'
# data_path = 'psd_data/'

erp_path = 'erp_data/'
if not exists (file_path):
    os.mkdir('results')
    os.mkdir(file_path)
if not exists (fig_path):
    os.mkdir(fig_path)
paths = {'results': results_path,
         'data'   : data_path,
         'erp'    : erp_path}
# =============================================================================
# Exatraction of a/periodic parameters using FOOOF model
# =============================================================================
good_channels = {}
shared_good_channel_indices = {}
param_list = {}
param_list_power = {}

excluded_subjects_dict = {}
retained_channels = {}
common_retained_std = {}
common_retained_dev = {}
# cond_to_analyze = 'both'

detect_stats = {} 
detectable_percent = {} 
fit_metric = {}

for drug in ['ket','psi']:
    # =============================================================================
    # Define EEG data features and variables
    # =============================================================================

    options = get_analysis_options(drug)

    # =============================================================================
    #     Loop FOOOF
    # =============================================================================
    filename = file_path+'fgs'+'_'+str(options['ntrials'])+'_'+drug
    if (options['overwrite_fooof']):
        print('Running FOOOF ...')
        fgs_fit = run_fooof(drug, paths, options, do_plot=True)
        
    else:
        print('Loading FOOOF ...')
        fgs_fit = load_fooof(filename, options, do_plot=True)
    
    # =============================================================================
    # # Extract aperiodic and band-specific parameters
    # =============================================================================
    print('Extracting parameters ...')
    param_list[drug], detect_stats[drug], detectable_percent[drug],\
                fit_metric[drug] =  loop_param_extract(fgs_fit, options, paths, do_plot=True)
    
    
    # =============================================================================
    #     FOOOF Model QC
    # =============================================================================
    print('Running qulaity control ...')
    df, df_retention = qc_plot_detectable(detect_stats[drug], detectable_percent[drug],
                                          options, paths, drug, do_plot=True)
    
    param_list[drug], retained_channels[drug], param_labels_filtered, good_channels[drug] = \
                                   qc_filter(param_list[drug], detectable_percent[drug], 
                                            fit_metric[drug], options)
    # =========================================================================
    # PLS analysis on FOOOF model's parameters -- Number of channels detectable
    # =========================================================================
    # Nchan = len(shared_good_channels[drug])
    # filename = file_path+"mpls_fooof_params_meancentre_"+str(mc)+"_"+str(options['ntrials'])+'_'+drug+".pkl"
    
    file_prefix = f'mpls_fooof_params_meancentre_{options["ntrials"]}_{drug}'
    print('Mean-centered PLS within subjects ...')
    mpls_fooof_dict = run_plot_mpls(param_list[drug], detect_stats[drug], 
                                    paths, param_labels_filtered, options,
                                    file_prefix=file_prefix,
                                    exclude_bad_fits=True, 
                                    good_channels=good_channels[drug], 
                                    only_retained=False, 
                                    common_retained=retained_channels[drug][options['cond_to_analyze']], 
                                    mc=2, do_plot=True, 
                                    analyze_percentage=True, 
                                    analyze_power=True,
                                    analyze_aperiodic=True,
                                    analyze_all=True, 
                                    do_topo=False,
                                    detectable_percent=detectable_percent[drug], drug=drug)
    
    print('behavioral PLS within subjects ...')
    bpls_fooof_dict = run_plot_bpls(param_list[drug], param_labels_filtered, options,
                      paths, drug=drug, do_predict=True, exclude_bad_fits=True, 
                      good_channels=good_channels[drug], only_retained=False,
                      common_retained=retained_channels[drug], 
                      only_frontal=True, do_plot=True)   

# =============================================================================
# Ketamine-Psilocybin comparison
# =============================================================================
print('Running QC between drugs ...')
common_detect_ketpsi, common_good_ketpsi = filter_common_channels_between_drugs(param_list, 
                                                                                good_channels,
                                                                                retained_channels, 
                                                                                options)
# options['cond_to_analyze'] = 'both'
# options['interact_mode'] = 'original'
param_ketpsi2 = prep_ket_psi_param_diff(param_list, options)
plot_violin_distributions(param_list, fig_path, param_labels_filtered, options)
file_prefix = f'mpls_fooof_params_meancentre_interact_{options["ntrials"]}_{drug}'
print('MPLS between drugs differences ....')
mpls_fooof_dict_drugdiff = run_plot_mpls(param_ketpsi2, detect_stats, 
                                         paths, param_labels_filtered, options,file_prefix=file_prefix,
                                         exclude_bad_fits=True, good_channels=common_good_ketpsi,
                                         only_retained=False, 
                                         common_retained=common_detect_ketpsi[options['cond_to_analyze']], 
                                         mc=2, do_plot=True,analyze_percentage=False, 
                                         analyze_power=True, analyze_aperiodic=True,
                                         analyze_all=True, drug='both')
# =============================================================================
# Correlation with the ASC measures
# =============================================================================

bpls_fooof_dict = run_plot_bpls(param_list, param_labels_filtered, options,
                  paths, do_predict=True, exclude_bad_fits=True, 
                  good_channels=common_good_ketpsi, only_retained=False,
                  common_retained=common_detect_ketpsi, only_frontal=True,
                  do_plot=True)



