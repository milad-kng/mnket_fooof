# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21 15:12:04 2025

@author: milad
"""
# =============================================================================
# Importage
# =============================================================================
from fooof import FOOOF
from fooof import FOOOFGroup
from fooof.objs import fit_fooof_3d, combine_fooofs
from fooof.bands import Bands

from fooof.plts.annotate import plot_annotated_model
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
import pickle
from os.path import isfile, exists
import pandas as pd
from collections import defaultdict
from aux_functions import extract_params_fg, check_nans, plot_impute_dist

from pyls import behavioral_pls
from pyls import meancentered_pls

from plot_helpers import plot_group_comparison, plot_heatmap, \
                         plot_topo_data, plot_behav_results, plot_psd_behavior_correlation
from scipy.stats import ttest_ind
from scipy.stats import ttest_rel
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression

# =============================================================================
# Run FOOOF
# =============================================================================
def run_fooof(drug, paths, options, do_plot=False):
    file_path = paths['results'] + '/files/'
    fig_path = paths['results'] + '/figs/'
    data_path = paths['data']
    print('Modeling '+drug+' PSD...')
    
    figname  = fig_path+'goodness_of_fit_for_'+str(options['ntrials'])+'_'+drug
    filename = file_path+'fgs'+'_'+str(options['ntrials'])+'_'+drug
    #---------Parameter extraction example---------#
    fm = FOOOF(min_peak_height=0.1, verbose=False,peak_width_limits=[1, 6],max_n_peaks=2)
    # spec = sio.loadmat(data_path+'psd_std_ket_roving_trial_series_csd_all.mat')['psd_std_ket']
    spec = sio.loadmat(data_path+'psd_std_pla2_roving_trial_series_csd_all.mat')['psd_std_pla2']
    X = spec[10,3,:]
    freqs = np.linspace(options['freq_range'][0], options['freq_range'][1], X.shape[-1])
    # if options['fit_domain'] == 'high':
    #     freq_mask = (freqs >= 20) & (freqs <= 40)
    # elif options['fit_domain'] == 'low':
    #     freq_mask = (freqs >= options['freq_range'][0]) & (freqs <= 20) 
    # else:
    #     # freq_mask = (freqs >= options['freq_range'][0]) & (freqs <= options['freq_range'][1])
    #     freq_mask = (freqs >= options['fit_range'][0]) & (freqs <= options['freq_range'][1])

   
    # fit_range = freqs[freq_mask]
    # X = X[freq_mask]
    fm.fit(freqs, X)
    fm.report(freqs,X,freq_range=options['fit_range'])
    if do_plot:
        fm.plot(plot_peaks='shade', peak_kwargs={'color' : 'green'})
        plot_annotated_model(fm, annotate_aperiodic=False)
        # ax = plt.gca()

        # Set y-axis limits
        plt.ylim([-1.5, 0.5])  # Adjust limits as needed
    
        # Save the figure
        plt.savefig('fooof_plot.png', dpi=600, bbox_inches='tight')  # or use .pdf/.svg
    
        plt.show()
    # ---------- Extracting FOOOF parameters for group of subjects ------------
    fg = FOOOFGroup(peak_width_limits=[1, 6], min_peak_height=0.1,#min_peak_height=0.15,
                            max_n_peaks=4, verbose=False)
    fgs_fit = {}
    for mode in options['mode_list']:
        print('Running FOOOF for mode: '+mode)
        spectra = sio.loadmat(data_path+'psd'+mode+'_roving_trial_series_csd_'+\
                              str(options['ntrials'])+'.mat')['psd'+mode]  
        # spectra = spectra[:, :, freq_mask]
        # freqs   = np.linspace(options['freq_range'][0], options['freq_range'][1]-5, np.shape(spectra)[-1])
        fgs_fit[mode] = fit_fooof_3d(fg, freqs, spectra, freq_range=options['fit_range'])
        fg_temp = combine_fooofs(fgs_fit[mode])
        if do_plot:
            fg_temp.plot()#(save_fig=True, file_name=filename, file_path=filepath)
            plt.savefig((figname+mode+'.jpg'))
        F = fgs_fit[mode]
        with open(filename + mode+'_fitting.pkl', 'wb') as f:
            pickle.dump(F, f) 
    return fgs_fit

def load_fooof(filename, options, do_plot=True):
    fgs_fit = {}
    for mode in options['mode_list']:
        with open(filename + mode+'_fitting.pkl', 'rb') as f:
            fgs_fit[mode] = pickle.load(f)
        fg_temp = combine_fooofs(fgs_fit[mode])
        if do_plot:
            fg_temp.plot()
    return fgs_fit
# =============================================================================
# Extracting parameters
# =============================================================================
def loop_param_extract(fgs_fit, options, paths, do_plot=True):
    aperiod_param  = {}
    period_param   = {}
    fit_metric     = {}
    detect_stats   = {}
    detectable_percent = {}
    param_list  = []
    
    for mode in options['mode_list']:
        # Extract periodic and aperiodic parameters
        
        [aperiod_param[mode], period_param[mode], fit_metric[mode], 
                     detect_stats[mode], detectable_percent[mode]] = \
                                    extract_params_fg(fgs_fit[mode], options)
        
        print('group error mean for mode '+mode+' :' + str(np.mean(fit_metric[mode][:,:,0]))+'+/-'+str(np.std(fit_metric[mode][:,:,0])))
        print('group r^2 mean for mode '+mode+' :' + str(np.mean(fit_metric[mode][:,:,1]))+'+/-'+str(np.std(fit_metric[mode][:,:,1])))
    
    # =========================================================================
    #     Replace NaN Values with the total_mean-2*std
    # =========================================================================
    if options['nan_policy'] == 'impute_normal' and\
        options['impute_param'] == 'group_mean':
            
        print('NaN values will be replaced by samples from a distribution with mean of mean-2*std of the group')

        band_labels = options['bands'].labels
        mean_per_subband = {}
        std_per_subband = {}
        param_prefix = ['_CF', '_PW', '_BW']
        for subband in band_labels:  # loop over subbands
            arrays = []
            for mode in options['mode_list']:
                arrays.append(period_param[mode][subband])  # shape (nsubs, nchans, nparams)
                    
            # concatenate across subjects dimension
            stacked = np.concatenate(arrays, axis=0)
            
            # mean over subs × chans, keep params
            mean_per_subband[subband] = np.nanmean(stacked, axis=(0,1))  # shape = (nparams,) 
            std_per_subband[subband]  = np.nanstd(stacked, axis=(0,1))  # shape = (nparams,) 
            if do_plot:
                for i_p in range(3):
                    if param_prefix[i_p] == '_CF' or param_prefix[i_p] == 'BW':
                        n_std = 0
                    elif param_prefix[i_p] == '_PW':
                        n_std = 2
                    plot_impute_dist(options['bands'], paths, 
                                     impute_mean=mean_per_subband[subband][i_p]-n_std*std_per_subband[subband][i_p],
                                     impute_std=std_per_subband[subband][i_p],
                                     label=f'{subband}{param_prefix[i_p]}',
                                     filename=f'{subband}{param_prefix[i_p]}')
            
                             
        for mode in options['mode_list']:
            for ind in range(options['nsubs']):
             
                for subband in band_labels:
                    
                        periodic_temp = period_param[mode][subband][ind,:,:]
                        # if subband == 'alpha' or subband == 'beta':
                        period_param[mode][subband][ind,:,:] = check_nans(periodic_temp,
                                                                    band_label=subband,
                                                                    options=options,
                                                                    impute_mean_pw=mean_per_subband[subband][1]-2*std_per_subband[subband][1] ,
                                                                    impute_std_pw=std_per_subband[subband][1],
                                                                    cf_mean=mean_per_subband[subband][0],
                                                                    cf_std=std_per_subband[subband][0],
                                                                    bw_mean=mean_per_subband[subband][2],
                                                                    bw_std=std_per_subband[subband][2])
           
    for mode in options['mode_list']:

        if options['fit_domain'] == 'high' :
            temp = np.concatenate([aperiod_param[mode],
                                    period_param[mode]['beta']], axis=2)
        elif options['fit_domain'] == 'broad':
            temp = np.concatenate([aperiod_param[mode],
                                    period_param[mode]['delta'],
                                    period_param[mode]['theta'],
                                    period_param[mode]['alpha'],
                                    period_param[mode]['beta']], axis=2)
        elif options['fit_domain'] == 'low':
            temp = np.concatenate([aperiod_param[mode],
                                    period_param[mode]['delta'],
                                    period_param[mode]['theta'],
                                    period_param[mode]['alpha'],
                                    period_param[mode]['beta']], axis=2)
        param_list.append(temp)

        print('group error mean for mode '+mode+' :' + str(np.mean(fit_metric[mode][:,:,0]))+'+/-'+str(np.std(fit_metric[mode][:,:,0])))
        print('group r^2 mean for mode '+mode+' :' + str(np.mean(fit_metric[mode][:,:,1]))+'+/-'+str(np.std(fit_metric[mode][:,:,1])))
    param_list = np.array(param_list)
    return param_list, detect_stats, detectable_percent, fit_metric
        
        
def run_plot_mpls(param_list, detect_stats, paths, param_labels, options, file_prefix=None,
                  exclude_bad_fits=True, good_channels=None,
                  only_retained=True, common_retained=None, mc=2, do_plot=True, 
                  analyze_percentage=True, analyze_power=True, analyze_aperiodic=True,
                  analyze_all=True, do_topo=False, detectable_percent=None,drug='both'):
    
    if drug == 'ket' or drug == 'psi':
        Nsubs = options['nsubs']
        reshaped_stats = []
        for mode in options['mode_list']:
            subbands = list(detect_stats[mode].keys())  # 4 subbands
            all_values = []

            for i, sb in enumerate(subbands):
                values = detect_stats[mode][sb]  # shape (19,)
                all_values.append(values.T)

            reshaped_stats.append(np.array(all_values ))
        reshaped_stats = np.array(reshaped_stats)
        reshaped_stats = np.transpose(reshaped_stats, (0, 2, 1)) 
    else:
        Nsubs = options['nket'] + options['npsi']
        reshaped_stats = None
    file_path = paths['results'] + '/files/'
    fig_path  = paths['results'] + '/figs/'
    
    bands = options['bands']
    chan_names = options['chan_names']
    if exclude_bad_fits:
        aperiodic_channel_indices = good_channels
        if only_retained:
            periodic_channel_indices = [ch for ch in common_retained if ch in good_channels]
            # periodic_channel_indices = [options['chan_names'][i] for i in common_good_indices]
        else:
            periodic_channel_indices = aperiodic_channel_indices
    else:
        aperiodic_channel_indices = range(len(chan_names))
        if only_retained:
            periodic_channel_indices = common_retained
        else:
            periodic_channel_indices = aperiodic_channel_indices

            
    
    if drug == 'ket' or drug == 'psi':
        analyses = {
            'percentage': {
                'overwrite': analyze_percentage,
                'param_to_analyze'   : ['alpha', 'beta'],
                'param_list': reshaped_stats,
                'n_cond' : 2
            },
            'periodic': {
                'overwrite': analyze_power,
                'param_to_analyze'   : param_labels[2::],#['Alpha_CF', 'Alpha_PW', 'Alpha_BW',
                                       # 'Beta_CF', 'Beta_PW', 'Beta_BW'],
                'param_list': param_list,
                'n_cond' : 2,
                'chan_idx' : periodic_channel_indices
            },
            'aperiodic': {
                'overwrite': analyze_aperiodic,
                'param_to_analyze'   : param_labels[0:2],
                'param_list': param_list,
                'n_cond' : 2,
                'chan_idx' : aperiodic_channel_indices
            },
            'all': {
                'overwrite': analyze_all,
                'param_to_analyze'   : param_labels,
                'param_list': param_list,
                'n_cond' : 2,
                'chan_idx' : periodic_channel_indices
            }
            }
        study = drug
    elif drug == 'both':
        if options['cond_to_analyze'] == 'both' or options['interact_mode'] == 'diff':
            n_cond = 2
            study = 'interact_diff'
        else:
            n_cond = 1
        if options['interact_mode'] == 'original':
            study = 'interact'

        analyses = {
            'periodic': {
                'overwrite': analyze_power,
                'param_to_analyze'   : param_labels[2::],#['Alpha_CF', 'Alpha_PW', 'Alpha_BW',
                                       # 'Beta_CF', 'Beta_PW', 'Beta_BW'],
                'param_list': param_list,
                'n_cond' : n_cond,
                'chan_idx' : periodic_channel_indices
            },
            'aperiodic': {
                'overwrite': analyze_aperiodic,
                'param_to_analyze'   : param_labels[0:2],
                'param_list': param_list,
                'n_cond' : n_cond,
                'chan_idx' : aperiodic_channel_indices
            },
            'all': {
                'overwrite': analyze_all,
                'param_to_analyze'   : param_labels,
                'param_list': param_list,
                'n_cond' : n_cond,
                'chan_idx' : periodic_channel_indices
            }
            }
        
    mpls_fooof_dict = {}
    for i_analysis, key in enumerate(analyses):
        
        filename  = f'{file_path}{file_prefix}_{key}.pkl'

        # Percentage
        if analyses[key]['overwrite'] and options['overwrite_pls_on_params'] or \
            not(isfile(filename)):
            print(f'PLS on {key}...')
            Nparam = len(analyses[key]['param_to_analyze'])
            if key == 'percentage' and drug != 'both':
                param_indices = [bands.labels.index(label) for label in analyses[key]['param_to_analyze']]
                param_to_pls  = np.take(analyses[key]['param_list'], 
                                        indices=param_indices, axis=-1)
                param_to_pls  = np.reshape(param_to_pls, (4*Nsubs, Nparam))
                groups = [Nsubs, Nsubs]
            else:
                if drug == 'ket' or drug == 'psi':
                    if key == 'periodic':
                        param_to_pls = analyses[key]['param_list'][:,:,:,2::]
                    if key == 'aperiodic':
                        param_to_pls = analyses[key]['param_list'][:,:,:,0:2]
                    if key == 'all':
                        param_to_pls = analyses[key]['param_list']
                    param_to_pls  = param_to_pls[:,:,analyses[key]['chan_idx'],:]
                    param_to_pls  = np.reshape(param_to_pls, 
                                               (4*Nsubs, 
                                                Nparam*len(analyses[key]['chan_idx'])))
                    groups = [Nsubs, Nsubs]
                elif drug == 'both':
                    if key == 'periodic':
                        param_to_pls = analyses[key]['param_list'][:,:,2::]
                    if key == 'aperiodic':
                        param_to_pls = analyses[key]['param_list'][:,:,0:2]
                    if key == 'all':
                        param_to_pls = analyses[key]['param_list']
                    
                    if options['interact_mode'] == 'diff':
                        groups = [options['nket'], options['npsi']]
                    if options['interact_mode'] == 'original':
                        groups = [options['nket'], options['nket'], 
                                  options['npsi'], options['npsi']]
                    param_to_pls  = param_to_pls[:,analyses[key]['chan_idx'],:]
                    param_to_pls  = np.reshape(param_to_pls, 
                                               (len(groups)*Nsubs, 
                                                Nparam*len(analyses[key]['chan_idx'])))
                    
            mpls_fooof = meancentered_pls(param_to_pls,
                                          mean_centering=mc, 
                                          n_perm=2000, n_boot=2000,
                                          groups= groups,
                                          n_cond=analyses[key]['n_cond'],
                                          n_proc='max',verbose=True)
            with open(filename, 'wb') as f:
                pickle.dump(mpls_fooof,f)
        else:
            with open(filename, 'rb') as f:
                mpls_fooof = pickle.load(f)
    
        mpls_fooof_dict[key] = mpls_fooof
        
        # Plotting
        if do_plot and key != 'percentage':
            labels = analyses[key]['param_to_analyze']
            channels_to_plot = [chan_names[i] for i in analyses[key]['chan_idx']]
            chan_idx_to_plot = analyses[key]['chan_idx']
            plot_group_comparison(mpls_fooof, fig_path, study=study, filename=key)
            
            plot_heatmap(mpls_fooof, labels,fig_path, channels_to_plot, 
                         chan_idx_to_plot, drug, filename=key)
    if drug != 'both':         
        if do_topo :
            mpls_fooof = mpls_fooof_dict['all']
            plot_topo_data(
                param_list, mpls_fooof, param_labels,fig_path, 
                options['ntrials'], Nsubs, detectable_percent, periodic_channel_indices, 
                label=drug, mode_list=options['mode_list'])
            
            mpls_fooof = mpls_fooof_dict['aperiodic']
            param_labels=['Offset', 'Exponent']
            plot_topo_data(
                param_list, mpls_fooof, param_labels,fig_path, 
                options['ntrials'], Nsubs, detectable_percent, aperiodic_channel_indices, 
                label=drug, mode_list=options['mode_list'])
        print(f'difference in the percentage of detectable channels across modes, p-value: {str(mpls_fooof_dict["percentage"].permres.pvals[0])}')
        print(f'# X_weights_normalized={str(mpls_fooof_dict["percentage"].bootres.x_weights_normed[:,0])}')
    
    return mpls_fooof_dict
    

def prep_ket_psi_param_diff(param_list, options):
    
    param_ket = param_list['ket']
    param_psi = param_list['psi']
    
    cond_to_analyze = options['cond_to_analyze']
    mode = options['interact_mode']
    if mode == 'diff':
        if cond_to_analyze == 'std':
            param_ketpsi = np.concatenate((
            param_ket[1, :] - param_ket[0, :], 
            param_psi[1, :] - param_psi[0, :]
        ))
    
        elif cond_to_analyze == 'dev':
            param_ketpsi = np.concatenate((
                param_ket[3, :] - param_ket[2, :], 
                param_psi[3, :] - param_psi[2, :]
            ))
        
        elif cond_to_analyze == 'both':
            param_ketpsi = np.concatenate((
                param_ket[1, :] - param_ket[0, :], 
                param_ket[3, :] - param_ket[2, :], 
                param_psi[1, :] - param_psi[0, :], 
                param_psi[3, :] - param_psi[2, :]
            ))
    
    if mode == 'original':
        param_ketpsi = np.concatenate((param_ket[0,:], param_ket[2,:],
                                        param_ket[1,:], param_ket[3,:],
                                        param_psi[0,:], param_psi[2,:],
                                        param_psi[1,:], param_psi[3,:]))
        # param_ketpsi = np.concatenate((
        #                                param_ket[1,:], param_ket[3,:],
        #                                param_psi[1,:], param_psi[3,:]))
    return param_ketpsi
    
def run_plot_bpls(param_list, param_labels_filtered, options,
                  paths, drug='both',do_predict=True, exclude_bad_fits=True, 
                  good_channels=None, only_retained=True,
                  common_retained=None, only_frontal=False, do_plot=True, do_anxiety=False):
    # =============================================================================
    # Correlations with the behavioural data
    # =============================================================================
    behav_drug, behav_pla = prep_behav_data(options)
    file_path = paths['results'] + '/files/'
    fig_path  = paths['results'] + '/figs/'
    chan_names = options['chan_names']

    frontal_channels =  ['Fp1', 'AF7', 'AF3', 'F1' ,'F3' ,'F5', 'F7',
                         'Fpz' ,'Fp2', 'AF8','AF4' ,'AFz' ,'Fz' ,'F2' ,
                         'F4' ,'F6', 'F8']
    frontal_indices = [chan_names.index(ch) for ch in frontal_channels if ch in chan_names]
    # param_to_do = [0,1,5,6,7,8,9,10] # aperiodic, alpha, beta
    # behav_to_do = [11] # Impaired control and cognition
    Nparam = len(param_labels_filtered)
    # Nbehav = len(behav_to_do)
    
    if drug == 'ket':
        behav_drug = behav_drug[0:options['nket']]
        behav_pla  = behav_pla [0:options['nket']]
        Nsubs = options['nket']
    elif drug == 'psi':
        behav_drug = behav_drug[options['nket']::]
        behav_pla  = behav_pla [options['nket']::]
        Nsubs = options['npsi']
    elif drug == 'both':
        Nsubs = options['nket'] + options['npsi']
    behav_vars = behav_pla.columns[2::]
    filename_single = f'{file_path}bpls_behav_spect_combined_{str(options["ntrials"])}_single_Behav_multChan_do_Predict_{str(do_predict)}_single_behav.pkl'        
    
    if do_anxiety:
        for drug in ['ket', 'psi']:
            if drug == 'ket':
                behav_drug_sub = behav_drug[0:options['nket']]
                behav_pla  = behav_pla [0:options['nket']]
                Nsubs = options['nket']
            elif drug == 'psi':
                behav_drug = behav_drug_sub[options['nket']::]
                behav_pla  = behav_pla [options['nket']::]
            Y = np.expand_dims(np.array(behav_drug_sub['Anxiety']),1)

            param_cond = np.reshape(param_list[drug], 
                                   (4, options[f'n{drug}'], 64, Nparam))
           
            # pla_X  = param_cond[0,:,:,:]
            X = np.reshape(param_cond[1,:,:,:],(options[f'n{drug}'], 64*8))
            axniety_bpls = behavioral_pls(X,Y, 
                            n_perm=2000, n_boot=2000,
                            n_proc='max', verbose=False)
            plot_psd_behavior_correlation(X, Y, range(options[f'n{drug}']), 
                                              xlabel='', 
                                              ylabel='',
                                              title='',
                                              fig_path=None, fname=None, nTrials=None)
    # behav_to_do = ['AV synesthia', 'Changed meaning of percept2',
    #                 'Elementry Imagry', 'Complex Imagery',
    #                 'Disembodiment', 'Experience of unity', 
    #                 'Impaired control and cognition']
    behav_to_do = [#'AV synesthia', 'Changed meaning of percept2',
                    # 'Complex Imagery',
                    'Disembodiment']
        #============ PLS-behavioral on the spectro-behavioral data ==============#
    for ib, behav in enumerate(behav_to_do):
        if options['overwrite_bpls']:
            print('## Performing single behavioral-multiple spectral varibale B-PLS...')

            # pla_Y  = np.array(behav_pla [behav_vars[behav_to_do]])
            drug_Y = np.array(behav_drug[behav])
            
            # pla_Y_modif = pla_Y + np.random.multivariate_normal(np.zeros(pla_Y.shape[1]),
            #                               0.000001*np.diagflat(pla_Y.std(axis=0)),Nsubs)
            drug_Y_modif = drug_Y #+ np.random.multivariate_normal(np.zeros(drug_Y.shape[1]),
                                   #       0.000001*np.diagflat(drug_Y.std(axis=0)),Nsubs)
    
            
            # behav_var = behav_vars[ib]#.tolist()
            channel_indices_dict = create_channel_list(chan_names, frontal_indices, good_channels, 
                                                    exclude_bad_fits=exclude_bad_fits, 
                                                    only_frontal=only_frontal,
                                                    only_retained=only_retained, 
                                                    common_retained=common_retained)
            bpls_single = {}
            i_pla = 0
            for stim in {'std', 'dev'}:
                channel_indices = channel_indices_dict[stim]
                Nchan  = len(channel_indices)
                if drug == 'both':
                    param_ket = np.reshape(param_list['ket'][:,:,channel_indices], 
                                           (4, options['nket'], Nchan, Nparam))
                    param_psi = np.reshape(param_list['psi'][:,:,channel_indices], 
                                           (4, options['npsi'], Nchan, Nparam))
                    pla_X  = np.concatenate((param_ket[i_pla,:,:,:], 
                                             param_psi[i_pla,:,:,:]),
                                            axis=0)
                    drug_X = np.concatenate((param_ket[i_pla+1,:,:,:],
                                             param_psi[i_pla+1,:,:,:]),
                                            axis=0)
                elif drug == 'ket' or drug == 'psi':
                    param_cond = np.reshape(param_list[:,:,channel_indices,:], 
                                           (4, options[f'n{drug}'], Nchan, Nparam))
                   
                    pla_X  = param_cond[i_pla,:,:,:]
                    drug_X = param_cond[i_pla+1,:,:,:]
                    
                # bpls_combined = defaultdict(dict)
    
    
                if do_predict:
                    # X = np.squeeze(np.reshape(pla_X [:,chan_idx,:],(N,Nc*Nparam)))
                    X = np.squeeze(np.reshape(pla_X ,(Nsubs,Nchan*Nparam)))
                else:
                    # X = np.squeeze(np.reshape(drug_X[:,chan_idx,:],(N,Nc*Nparam)))
                    X = np.squeeze(np.reshape(drug_X,(Nsubs,Nchan*Nparam)))
                X = X + np.random.multivariate_normal(np.zeros(X.shape[1]), # Add a small noise to avoid devision by zero
                                    0.000001*np.diagflat(X.std(axis=0)),Nsubs)
               
                Y = np.expand_dims(drug_Y_modif,1) #- pla_Y_modif
                print(f'behavioral variable: {behav} in {stim}')
                bpls_temp = behavioral_pls(X,Y, 
                                           n_perm=2000, n_boot=2000,
                                           n_proc='max', verbose=False)
    
                # print('+++ significant correlation--')
                print('++++ p_val='+str(bpls_temp.permres.pvals[0]))
                
                bpls_single[stim] = bpls_temp
                i_pla = i_pla + 2
            
            with open(filename_single, 'wb') as f:
                pickle.dump(bpls_single,f)
       
        else:
            with open(filename_single, 'rb') as f:
                bpls_single = pickle.load(f)
            channel_indices_dict = create_channel_list(chan_names, frontal_indices, good_channels, 
                                                    exclude_bad_fits=exclude_bad_fits, 
                                                    only_frontal=only_frontal,
                                                    only_retained=only_retained, 
                                                    common_retained=common_retained)
          
        if do_plot:
            for stim in ['std', 'dev']:
                channel_indices = channel_indices_dict[stim]
                Nchan  = len(channel_indices)

                if drug == 'both':
                    param_ket = np.reshape(param_list['ket'][:,:,channel_indices], 
                                           (4, options['nket'], Nchan, Nparam))
                    param_psi = np.reshape(param_list['psi'][:,:,channel_indices], 
                                           (4, options['npsi'], Nchan, Nparam))
                    pla_X  = np.concatenate((param_ket[i_pla,:,:,:], 
                                             param_psi[i_pla,:,:,:]),
                                            axis=0)
                    drug_X = np.concatenate((param_ket[i_pla+1,:,:,:],
                                             param_psi[i_pla+1,:,:,:]),
                                            axis=0)
                elif drug == 'ket' or drug == 'psi':
                    param_cond = np.reshape(param_list[:,:,channel_indices,:], 
                                           (4, options[f'n{drug}'], Nchan, Nparam))
                   
                    pla_X  = param_cond[0,:,:,:]
                    drug_X = param_cond[0+1,:,:,:]
                X_reg = np.squeeze(bpls_single[stim].x_scores[:,0])
                Y_reg = np.squeeze(bpls_single[stim].y_scores[:,0]*\
                                    bpls_single[stim].y_weights)
        
                # Plotting --------------------------------------------------------------------
                # if DoPlot:
                selected_channel_names = [chan_names[i] for i in sorted(channel_indices)]
            
                X_psd = pla_X[:, 3, 2]#.mean(axis=1)
                plot_behav_results(X_reg, Y_reg.ravel(), fig_path, options['ntrials'], 
                                   bpls_single[stim], 
                                   select_channels=selected_channel_names, 
                                   param_labels=param_labels_filtered, drug=drug,
                                   plot_orig_correlation=True, 
                                   X_psd=X_psd, Y_behav=Y_reg.ravel(), 
                                   behav_label=f'{behav} in {stim} for {drug}')
    # i_pla = 0
    # if do_plot:
    #     for stim in ['std', 'dev']:
    #         X_reg = np.squeeze(bpls_single[stim].x_scores[:,0])
    #         Y_reg = np.squeeze(bpls_single[stim].y_scores[:,0]*\
    #                             bpls_single[stim].y_weights)
    
    #         # Plotting --------------------------------------------------------------------
    #         # if DoPlot:
    #         selected_channel_names = [chan_names[i] for i in sorted(channel_indices)]
        
    #         X_psd = X_reg#pla_X[:, 9, 3]#.mean(axis=1)
    #         plot_behav_results(X_reg, Y_reg.ravel(), fig_path, options['ntrials'], 
    #                            bpls_single, 
    #                            select_channels=selected_channel_names, 
    #                            param_labels=param_labels_filtered, drug=drug,
    #                            plot_orig_correlation=True, X_psd=X_psd, Y_behav=Y_reg.ravel())
    return bpls_single
def prep_behav_data(options):
# =============================================================================
#     
# =============================================================================
    behav_path = 'D:/Science/Coding/mnket_fooof_paper/mnket_fooof/mnket_behavioural_data/Q_data.xlsx'#'PROVIDE ASC DATA PATH'
    labels = pd.read_excel(behav_path,header=1,skiprows=0,nrows=1)
    
    # Read ASC data
    behav_pla  = pd.read_excel(behav_path,names=labels,skiprows=1,nrows=35)
    behav_drug = pd.read_excel(behav_path,names=labels,skiprows=36,nrows=35)
    # replacing first and second variable to match the ordering of data #######
    df = behav_pla
    df.iloc[0], df.iloc[1] = df.iloc[1].copy(), df.iloc[0].copy()
    
    df = behav_drug
    df.iloc[0], df.iloc[1] = df.iloc[1].copy(), df.iloc[0].copy() 
    
    return behav_drug, behav_pla
# =============================================================================
# MMN Analysis
# =============================================================================
def mmn_vs_psd_analysis(options, erp_path ,mpls_fooof):
    ChanNames = options['chan_names']
    Nket = options['nket']
    Npsi = options['npsi']
    if options['do_mmn_analysis']:
        selected_channels = ['Cz', 'Fz', 'FCz']  # Channels to plot
        selected_idx = [ChanNames.index(ch) for ch in selected_channels]
        
        mmn_drug = {}
        mmn_pla  = {}
        
        for drug in ['ket', 'psi']:
            erp_high_drug = sio.loadmat(f'{erp_path}mmn_high_{drug}_roving.mat')[f'mmn_high_drug']
            erp_low_drug  = sio.loadmat(f'{erp_path}mmn_low_{drug}_roving.mat')[f'mmn_low_drug']
            mmn_d = erp_high_drug - erp_low_drug  # shape: subj × chan × time
            mmn_drug[drug] = mmn_d
        
            pla = 'pla1' if drug == 'ket' else 'pla2'
            erp_high_pla = sio.loadmat(f'{erp_path}mmn_high_{pla}_roving.mat')[f'mmn_high_pla']
            erp_low_pla  = sio.loadmat(f'{erp_path}mmn_low_{pla}_roving.mat')[f'mmn_low_pla']
            mmn_p = erp_high_pla - erp_low_pla
            mmn_pla[pla] = mmn_p
    
        # --- Average across subjects ---
        mmn_avg = {
            'Ketamine': -np.mean(mmn_drug['ket'], axis=0),  # shape: chan × time
            'Placebo1': np.mean(mmn_pla['pla1'], axis=0),
            'Psilocybin': np.mean(mmn_drug['psi'], axis=0),
            'Placebo2': np.mean(mmn_pla['pla2'], axis=0)
        }
        
        # --- Plotting ---
        time = np.arange(mmn_avg['Ketamine'].shape[1])  # X-axis (adjust if you have time in ms)
        
        for ch_idx in selected_idx:
            plt.figure(figsize=(10, 5))
            for label, color, ls in zip(mmn_avg.keys(), ['blue', 'gray', 'green', 'black'], ['-', '--', '-', '--']):
                plt.plot(time, mmn_avg[label][ch_idx], label=label, color=color, linestyle=ls)
            plt.title(f'MMN at {ChanNames[ch_idx]}')
            plt.xlabel('Time (samples or ms)')
            plt.ylabel('Amplitude (µV)')
            plt.ylim([-2,1.5])
            plt.axvline(0, color='k', linestyle=':')  # Optional: mark stimulus onset
            plt.legend()
            plt.grid(True)
            plt.tight_layout()
            plt.show()
        
        
        time_window_ms = [100, 200]
        time_vector = np.linspace(-100,400,mmn_pla['pla1'].shape[2])
        samples = np.searchsorted(time_vector, time_window_ms)
        
        
        placebo_mmn = mmn_pla['pla1'][:, :, samples[0]:samples[1]].mean(axis=2)
        ketamine_mmn = -mmn_drug['ket'][:, :, samples[0]:samples[1]].mean(axis=2)
     
        for iCh, ch_name in enumerate(selected_channels):
            t_stat, p_val = ttest_rel(placebo_mmn[:, iCh], ketamine_mmn[:, iCh])
            # Convert to one-tailed: test if placebo < ketamine (more negative)
            if t_stat > 0:
                p_val = 1 - p_val / 2
            else:
                p_val = p_val / 2
            print(f"\n{ch_name} MMN amplitude difference, p-value = {p_val:.4f}")
        
        # ket_mmn_diff = ketamine_mmn - placebo_mmn
        ket_mmn_diff = -ketamine_mmn - placebo_mmn
    
        placebo_mmn = mmn_pla['pla2'][:, :, samples[0]:samples[1]].mean(axis=2)
        psilocybin_mmn = mmn_drug['psi'][:, :, samples[0]:samples[1]].mean(axis=2)
        
        psi_mmn_diff = psilocybin_mmn - placebo_mmn
        
        Nsubs = Nket+Npsi
        mmn_diff_ketpsi = np.concatenate((ket_mmn_diff, 
                                          psi_mmn_diff), axis=0)
        mmn_to_pls = np.reshape(mmn_diff_ketpsi, (Nsubs, 64))
        mpls_mmn = meancentered_pls(mmn_to_pls,
                                            mean_centering=2, n_perm=2000,
                                            n_boot=2000, groups= [Nket,Npsi],
                                            n_cond=1, 
                                            n_proc='max')
        lr = LinearRegression()
        drug_labels = np.concatenate((np.repeat(1,Nket), np.repeat(0,Npsi)))
        le = LabelEncoder()
        y = le.fit_transform(drug_labels) 
        lr.fit(mpls_mmn.x_scores[:,0].reshape(-1,1), y)
        r2_mmn = lr.score(mpls_mmn.x_scores[:,0].reshape(-1,1), y)
        
        lr.fit(mpls_fooof.x_scores[:,0].reshape(-1,1), y)
        r2_fooof = lr.score(mpls_fooof.x_scores[:,0].reshape(-1,1), y)
        
        delta_r2 = r2_fooof - r2_mmn
        print(f"ΔR² = {delta_r2:.4f}")

def create_channel_list(chan_names, frontal_indices, good_channels, 
                        exclude_bad_fits=True, only_frontal=True,
                        only_retained=False, common_retained=None):
    channel_indices = {}
    for i, stim in enumerate(['std', 'dev']):
        if exclude_bad_fits:
            if only_retained:
                if only_frontal:
                    channel_indices[stim] = [ch for ch in common_retained[stim]
                                       if ch in good_channels and ch in frontal_indices]
                else:
                    channel_indices[stim] = channel_indices = [ch for ch in common_retained[stim] if ch in good_channels]

            else:
                if only_frontal:
                    channel_indices[stim] = [ch for ch in good_channels 
                                       if ch in frontal_indices]
                else:
                    channel_indices[stim] = good_channels
        else:
            if only_retained:
                if only_frontal:
                    channel_indices[stim] = [ch for ch in common_retained[stim]
                                       if ch in frontal_indices]
                else:
                    channel_indices[stim] = common_retained[stim]
            else:
                if only_frontal:
                    channel_indices[stim] = frontal_indices
                else:
                    channel_indices[stim] = range(len(chan_names))
    
    
    return channel_indices