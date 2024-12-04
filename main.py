# -*- coding: utf-8 -*-
"""
The analysis pipeline to analyze the power spectral density features of ketamine and 
psilocybin datasets obtained using fooof model.

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
from   collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib
# from   matplotlib import cm, colors, colorbar
# import matplotlib.patches as patches
import seaborn as sns
from matplotlib.colors import TwoSlopeNorm
# Import the FOOOF object
from fooof import FOOOF
from fooof import FOOOFGroup
from fooof.objs import fit_fooof_3d, combine_fooofs
from fooof.bands import Bands
# from fooof.analysis import get_band_peak_fg,get_band_peak_fm 
# Import plotting function for model parameters and components
# from fooof.plts.periodic import plot_peak_fits, plot_peak_params
# from fooof.plts.aperiodic import plot_aperiodic_params, plot_aperiodic_fits
from fooof.plts.annotate import plot_annotated_model
# MNE
import mne
from mne.viz import plot_topomap
# from mne.time_frequency import psd_welch
# Import Pyls
from pyls import behavioral_pls
from pyls import meancentered_pls
# from pyls import pls_regression
# import modeling objects
# from sklearn.linear_model import LinearRegression
# from sklearn.cross_decomposition import PLSRegression
# import scipy
# from scipy import stats
# import statsmodels
import os
from os.path import isfile, exists
from aux_functions import extract_params_fg,check_nans,\
                          plot_pls_analysis, plot_group_comparison, plot_example_psd,\
                              plot_violin_distributions, plot_behav_results
from scipy.stats import pearsonr
# from scipy.stats import ttest_ind
# from statsmodels.multivariate.manova import MANOVA
# =============================================================================
# Set global plotting attributes
# =============================================================================
plt.rcParams['axes.edgecolor'] = 'black'
plt.rcParams['axes.linewidth'] = 1.5  # Adjust the thickness to your liking
font = {'family': 'Times New Roman',
        'weight': 'normal'}  

# Apply font properties
plt.rc('font', **font)
sns.set(font='Times New Roman') 
# =============================================================================
# Define EEG data features and variables
# =============================================================================

ChanNames = ['Fp1', 'AF7', 'AF3', 'F1' ,'F3' ,'F5', 'F7' ,'FT7', 'FC5', 'FC3', 'FC1', 'C1',\
             'C3' ,'C5' ,'T7' ,'TP7' ,'CP5' ,'CP3' ,'CP1' ,'P1' ,'P3' ,'P5' ,'P7', 'P9',\
             'PO7', 'PO3' ,'O1', 'Iz', 'Oz', 'POz' ,'Pz', 'CPz' ,'Fpz' ,'Fp2', 'AF8',\
             'AF4' ,'AFz' ,'Fz' ,'F2' ,'F4' ,'F6', 'F8' ,'FT8' ,'FC6', 'FC4' ,'FC2',\
             'FCz' ,'Cz', 'C2', 'C4', 'C6' ,'T8', 'TP8', 'CP6', 'CP4', 'CP2', 'P2', 'P4',\
             'P6' ,'P8' ,'P10' ,'PO8' ,'PO4' ,'O2']

Nchan      = len(ChanNames)
freq_range = [0.5, 30]
sfreq = 256
nTrials = 'all'
bands = Bands({'delta' : [1, 4],
               'theta' : [4, 8],
               'alpha' : [8, 13],
               'beta'  : [13, 30]})
param_labels = ['Offset','Exponent','Delta_CF','Delta_PW','Delta_BW',\
                'Theta_CF','Theta_PW','Theta_BW','Alpha_CF',\
                'Alpha_PW','Alpha_BW','Beta_CF','Beta_PW','Beta_BW']
subs_psi    = ['3621', '4415', '4418','4419', '4420', '4421', '4426','4332', '4433', \
               '4460', '4476', '4502', '4515', '4518', '4591','4592']
subs_ket    = ['4431', '4446', '4447', '4458', '4482', '4487', '4488', '4548', '4494',\
               '4499', '4500', '4520', '4532', '4497', '4534', '4459', '4507','4422','4478']
Nket = len(subs_ket)
Npsi = len(subs_psi)

# =============================================================================
# Select Analysis
# =============================================================================
wght_th = 2.5 # Bootstap weights threshold

DoFOOOF = 0
DoPLSonParams = 0
DoPlot = 1
DoPLSonKetPsi = 0
DoPLSonBehavSpect = 1
DoPredict = 0

mc = 0

# =============================================================================
# Setup Paths
# =============================================================================
results_path = 'results/files/'
fig_path = 'results/figs/'
data_path = 'psd_data/'
if not exists (results_path):
    os.mkdir('results')
    os.mkdir(results_path)
if not exists (fig_path):
    os.mkdir(fig_path)
    
# =============================================================================
# Exatraction of a/periodic parameters using FOOOF model
# =============================================================================

param_list = {}

for drug in ['ket','psi']:
    print('Modeling '+drug+' PSD...')
    if drug == 'ket':
        subs = ['4431', '4446', '4447', '4458', '4482', '4487', '4488', '4548', '4494',\
               '4499', '4500', '4520', '4532', '4497', '4534', '4459', '4507','4422','4478']
        Nsubs = 19
        mode_list = ['_std_pla', '_std_ket', '_dev_pla', '_dev_ket']

    if drug == 'psi':
        print('')
        subs  = ['3621', '4415', '4418','4419', '4420', '4421', '4426','4332', '4433', \
                 '4460', '4476', '4502', '4515', '4518', '4591','4592']
        Nsubs = 16
        mode_list = ['_std_pla2', '_std_psi', '_dev_pla2', '_dev_psi']
    figname  = fig_path+'goodness_of_fit_for_'+str(nTrials)+'_'+drug
    filename = results_path+'fgs'+'_'+str(nTrials)+'_'+drug
    if (DoFOOOF):
        
        #---------Parameter extraction example---------#
        fm = FOOOF(min_peak_height=0.3, verbose=False,peak_width_limits=[1, 6],max_n_peaks=6)
        spec = sio.loadmat(data_path+'psd_std_pla2_roving_trial_series_csd_all.mat')['psd_std_pla2']
        X = spec[0,0,:]
        freqs = np.linspace(freq_range[0], freq_range[1], np.shape(X)[-1])
        fm.fit(freqs, X)
        fm.report(freqs,X,freq_range)
        if DoPlot:
            fm.plot(plot_peaks='shade', peak_kwargs={'color' : 'green'})
            plot_annotated_model(fm, annotate_aperiodic=False)
            ax = plt.gca()
        
        # ---------- Extracting FOOOF parameters for group of subjects ------------
        fg = FOOOFGroup(peak_width_limits=[1, 6], min_peak_height=0.3,#min_peak_height=0.15,
                                max_n_peaks=6, verbose=False)
        fgs_fit = {}
        for mode in mode_list:
            print('Running FOOOF for mode: '+mode)
            spectra = sio.loadmat(data_path+'psd'+mode+'_roving_trial_series_csd_'+str(nTrials)+'.mat')['psd'+mode]        
            freqs   = np.linspace(freq_range[0], freq_range[1], np.shape(spectra)[-1])
            fgs_fit[mode] = fit_fooof_3d(fg, freqs, spectra)
            fg_temp = combine_fooofs(fgs_fit[mode])
            if DoPlot:
                fg_temp.plot()#(save_fig=True, file_name=filename, file_path=filepath)
                plt.savefig((figname+mode+'.jpg'))
            F = fgs_fit[mode]
            with open(filename + mode+'_fitting.pkl', 'wb') as f:
                pickle.dump(F, f) 
    else:
        fgs_fit = {}
        for mode in mode_list:
            with open(filename + mode+'_fitting.pkl', 'rb') as f:
                fgs_fit[mode] = pickle.load(f)
            fg_temp = combine_fooofs(fgs_fit[mode])
            if DoPlot:
                fg_temp.plot()
    
    
    # =============================================================================
    # # Extract band-specific parameters
    # =============================================================================
        
    aperiod_param  = {}
    period_param   = {}
    fit_metric     = {} 
    param_temp = []
        
    for mode in mode_list:
        # Extract periodic and aperiodic parameters
        [aperiod_param[mode], period_param[mode], fit_metric[mode]] = extract_params_fg(fgs_fit[mode],bands,Nsubs)
 
        temp = np.concatenate([aperiod_param[mode],
                                period_param[mode]['delta'],
                                period_param[mode]['theta'],
                                period_param[mode]['alpha'],
                                period_param[mode]['beta']],axis=2)
        Nparam = np.size(temp,axis=2)
        param_temp.append(temp)
        print('group error mean for mode '+mode+' :' + str(np.mean(fit_metric[mode][:,:,0]))+'+/-'+str(np.std(fit_metric[mode][:,:,0])))
        print('group r^2 mean for mode '+mode+' :' + str(np.mean(fit_metric[mode][:,:,1]))+'+/-'+str(np.std(fit_metric[mode][:,:,1])))
        param_list[drug] = np.array(param_temp)
    
    # =========================================================================
    # PLS analysis on FOOOF model's parameters
    # =========================================================================
    
    filename = results_path+"mpls_fooof_params_meancentre_"+str(mc)+"_"+str(nTrials)+'_'+drug+".pkl"
    if (DoPLSonParams):
        Nparam = param_list[drug].shape[3]
        param_to_pls = np.reshape(param_list[drug],(4*Nsubs, Nchan*Nparam))
        mpls_fooof_params = meancentered_pls(param_to_pls,
                                             mean_centering=mc, n_perm=2000,
                                             n_boot=2000, groups= [Nsubs, Nsubs], n_cond=2, 
                                             n_proc='max',verbose=False)
    
        with open(filename, 'wb') as f:
            pickle.dump(mpls_fooof_params,f)
    else:
        with open(filename, 'rb') as f:
            mpls_fooof_params = pickle.load(f)
    
    Nparam = int(np.size(mpls_fooof_params.x_weights,axis=0)/Nchan)
    x_w = np.copy(np.reshape(mpls_fooof_params.bootres.x_weights_normed[:,0],(Nchan,Nparam)))
    if not(np.argwhere(np.isnan(x_w))) == []:
        print('nan Weights exist')
        check_nans(x_w, 'zero')
    
    x_w_df = pd.DataFrame(np.transpose(x_w), index = param_labels, columns = ChanNames)
    x_w_df[np.abs(x_w_df)< wght_th] = 0
    
    # =============== Plot PLS Results (bootstrap ratios,contrast,topomap) =================#
    if DoPlot:
        plot_pls_analysis(
            mpls_fooof_params, 
            fig_path, 
            nTrials, 
            drug, 
        )
        
# =============================================================================
# Ketamine-Psilocybin comparison
# =============================================================================
param_ket = param_list['ket']
param_psi = param_list['psi']

filename = results_path+'mpls_fooof_params_interact_meancentre_'+str(mc)+'_'+\
    str(nTrials)+'_'+drug+'.pkl'
if DoPLSonKetPsi:
    print('Analyzing ketamine-psilocybin-placebo comparison...')
    Nparam = param_list['ket'].shape[3]
  
    param_ketpsi = np.concatenate((param_ket[0,:], param_ket[2,:],
                                   param_ket[1,:], param_ket[3,:],
                                   param_psi[0,:], param_psi[2,:],
                                   param_psi[1,:], param_psi[3,:]))
    Nsubs = Nket+Npsi
    param_to_pls = np.reshape(param_ketpsi, (Nsubs*4, Nchan*Nparam))
    mpls_fooof_params_interact = meancentered_pls(param_to_pls,
                                                  mean_centering=1, n_perm=2000,
                                                  n_boot=2000, groups= [19,19,16,16],
                                                  n_cond=2, 
                                                  n_proc='max')
    with open(filename, 'wb') as f:
        pickle.dump(mpls_fooof_params_interact,f)
else:
    with open(filename, 'rb') as f:
            mpls_fooof_params_interact = pickle.load(f)

if DoPlot:
    plot_group_comparison(mpls_fooof_params_interact, fig_path)
    plot_example_psd(data_path, fig_path)
    plot_violin_distributions(param_list['ket'], param_labels, fig_path, mode_list,Nsubs)
      
# =============================================================================
# Correlations with the behavioural data
# =============================================================================
behav_path = 'PROVIDE ASC DATA PATH'
labels = pd.read_excel(behav_path,header=1,skiprows=0,nrows=1)

filename_single_std = results_path+"bpls_behav_spect_combined_std_"+str(nTrials)+\
    "_single_Behav_multChan_doPredict_"+str(DoPredict)+"_comb_behav.pkl"
filename_single_dev = "bpls_behav_spect_combined_dev_"+str(nTrials)+\
    "_single_Behav_multChan_doPredict_"+str(DoPredict)+"comb_pla.pkl"

front_ch =  ['Fp1', 'AF7', 'AF3', 'F1' ,'F3' ,'F5', 'F7',
              'Fpz' ,'Fp2', 'AF8','AF4' ,'AFz' ,'Fz' ,'F2' ,'F4' ,'F6', 'F8']

param_to_do = [0,1,8,9,10,11,12,13] # aperiodic, alpha, beta
behav_to_do = [11] # Impaired control and cognition

param_ket = np.reshape(param_ket, (4, Nket, Nchan, Nparam))
param_psi = np.reshape(param_psi, (4, Npsi, Nchan, Nparam))

# Read ASC data
behav_pla  = pd.read_excel(behav_path,names=labels,skiprows=1,nrows=35)
behav_drug = pd.read_excel(behav_path,names=labels,skiprows=36,nrows=35)
# replacing first and second variable to match the ordering of data #######
df = behav_pla
df.iloc[0], df.iloc[1] = df.iloc[1].copy(), df.iloc[0].copy()

df = behav_drug
df.iloc[0], df.iloc[1] = df.iloc[1].copy(), df.iloc[0].copy()
###########################################################################
N = 35

behav_vars = behav_pla.columns[2::]
                      
if DoPLSonBehavSpect:
    
    #============ PLS-behavioral on the spectro-behavioral data ==============#

    pla_Y = np.array(behav_pla [behav_vars[behav_to_do]])
    drug_Y = np.array(behav_drug[behav_vars[behav_to_do]])
    
    pla_Y_modif = pla_Y + np.random.multivariate_normal(np.zeros(pla_Y.shape[1]),
                                  0.000001*np.diagflat(pla_Y.std(axis=0)),N)
    drug_Y_modif = drug_Y + np.random.multivariate_normal(np.zeros(drug_Y.shape[1]),
                                  0.000001*np.diagflat(drug_Y.std(axis=0)),N)
    Nparam = len(param_to_do)
    Nbehav = len(behav_to_do)

    par_ket = param_ket[:,:,:,param_to_do]
    par_psi = param_psi[:,:,:,param_to_do]
    for i in [0,2]: # Std,Dev
    
        pla_X  = np.concatenate((par_ket[i,:,:,:],  par_psi[i,:,:,:]),axis=0)
        drug_X = np.concatenate((par_ket[i+1,:,:,:],par_psi[i+1,:,:,:]),axis=0)
        bpls_combined = defaultdict(dict)

        print('## Performing single behavioral-multiple spectral varibale B-PLS ...')

        front_idx = [ChanNames.index(key) for key in front_ch]
        Nc = len(front_idx)
        if DoPredict:
            X = np.squeeze(np.reshape(pla_X [:,front_idx,:],(N,Nc*len(param_to_do))))
        else:
            X = np.squeeze(np.reshape(drug_X[:,front_idx,:],(N,Nc*len(param_to_do))))
        X = X + np.random.multivariate_normal(np.zeros(X.shape[1]), # Add a small noise to avoid devision by zero
                            0.000001*np.diagflat(X.std(axis=0)),N)
        if np.any(X.std(axis=0) == 0):
            idx_zero = np.argwhere(X.std(axis=0) == 0)
            temp = np.squeeze(X[:,idx_zero])
            temp = temp + \
                np.random.multivariate_normal(np.zeros(X[:,idx_zero].shape[1]), # Add a small noise to avoid devision by zero
                        0.000001*np.diagflat(0.0000001+X[:,idx_zero].std(axis=0)),len(subs))
            X[:,idx_zero] = np.expand_dims(temp,axis=2)
        b = behav_to_do
        Y = np.expand_dims(drug_Y_modif[:,b],1) #- pla_Y_modif
        print('behavioral variable:'+behav_vars[behav_to_do[b]])
        bpls_temp = behavioral_pls(X,Y,n_perm=2000, n_boot=2000,n_proc='max', verbose=False)
        bpls_combined[behav_vars[behav_to_do[b]]] = bpls_temp

        print('+++ significant correlation--')
        print('++++ p_val='+str(bpls_temp.permres.pvals[0]))
        if i == 0:
            print('standard group done')
            bpls_single_std = bpls_combined
            with open(filename_single_std, 'wb') as f:
                pickle.dump(bpls_single_std,f)
        else:
            print('deviant group done')
            bpls_single_dev = bpls_combined
            with open(filename_single_dev, 'wb') as f:
                pickle.dump(bpls_single_dev,f)
else:
    with open(filename_single_std, 'rb') as f:
        bpls_single_std = pickle.load(f)
    with open(filename_single_dev, 'rb') as f:
        bpls_single_dev = pickle.load(f)
        
X_reg = np.squeeze(bpls_single_std['Impaired control and cognition'].x_scores[:,0])
Y_reg = np.squeeze(bpls_single_std['Impaired control and cognition'].y_scores[:,0]*\
                   bpls_single_std['Impaired control and cognition'].y_weights)

# Plotting --------------------------------------------------------------------
if DoPlot:
    plot_behav_results(X_reg, Y_reg, fig_path, nTrials, 
                       param_to_do, bpls_single_std, bpls_single_dev)

