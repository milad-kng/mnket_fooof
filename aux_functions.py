# -*- coding: utf-8 -*-
"""
Created on Fri Aug 22 13:36:30 2025

@author: milad
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 20:35:00 2024

@author: milad
"""
import numpy as np

import matplotlib.pyplot as plt

from fooof.bands import Bands
from fooof.analysis import get_band_peak_fg, get_band_peak_fm 

def get_analysis_options(drug):
    ChanNames = ['Fp1', 'AF7', 'AF3', 'F1' ,'F3' ,'F5', 'F7' ,'FT7', 'FC5', 'FC3', 'FC1', 'C1',
                 'C3' ,'C5' ,'T7' ,'TP7' ,'CP5' ,'CP3' ,'CP1' ,'P1' ,'P3' ,'P5' ,'P7', 'P9',
                 'PO7', 'PO3' ,'O1', 'Iz', 'Oz', 'POz' ,'Pz', 'CPz' ,'Fpz' ,'Fp2', 'AF8',
                 'AF4' ,'AFz' ,'Fz' ,'F2' ,'F4' ,'F6', 'F8' ,'FT8' ,'FC6', 'FC4' ,'FC2',
                 'FCz' ,'Cz', 'C2', 'C4', 'C6' ,'T8', 'TP8', 'CP6', 'CP4', 'CP2', 'P2', 'P4',
                 'P6' ,'P8' ,'P10' ,'PO8' ,'PO4' ,'O2']
    
    subs_psi = ['3621', '4415', '4418','4419', '4420', '4421', '4426','4332', '4433',
                '4460', '4476', '4502', '4515', '4518', '4591','4592']
    subs_ket = ['4431', '4446', '4447', '4458', '4482', '4487', '4488', '4548', '4494',
                '4499', '4500', '4520', '4532', '4497', '4534', '4459', '4507','4422','4478']
    if drug == 'ket':
        subs = subs_ket
        Nsubs = len(subs)
        mode_list = ['_std_pla', '_std_ket', '_dev_pla', '_dev_ket']
    
    if drug == 'psi':
        subs  = subs_psi
        Nsubs = len(subs)
        mode_list = ['_std_pla2', '_std_psi', '_dev_pla2', '_dev_psi']
    param_labels = ['Offset','Exponent','Delta_CF','Delta_PW','Delta_BW',
                    'Theta_CF','Theta_PW','Theta_BW','Alpha_CF',
                    'Alpha_PW','Alpha_BW','Beta_CF','Beta_PW','Beta_BW']
    params_to_remove=['Delta_CF', 'Delta_PW', 'Delta_BW',
                      'Theta_CF', 'Theta_PW', 'Theta_BW']
                       # 'Alpha_CF', 'Alpha_PW','Alpha_BW', ]
    #                   'Beta_BW', 'Beta_PW','Beta_CF']
    # params_to_remove=['Offset', 'Exponent']
    fit_domain = 'broad'
    # Define frequency bands using the Bands class
    if fit_domain == 'high':
        fit_range = [20, 40]
        bands = Bands({ 'beta'  : [20, 30]})
    elif fit_domain == 'low':
        fit_range = [1, 20]

        bands = Bands({'delta' : [1, 4],
                       'theta' : [4, 8],
                       'alpha' : [8, 13],
                       'beta'  : [13, 20]})
    elif fit_domain == 'broad':
        fit_range = [0.5, 30]
        bands = Bands({'delta' : [1, 4],
                       'theta' : [4, 8],
                       'alpha' : [8, 13],
                       'beta'  : [13,30]})
    cond_to_analyze = 'both'
    if cond_to_analyze == 'std':
        cond_select = np.arange(0, 2)
    if cond_to_analyze == 'dev':
        cond_select = np.arange(2, 4)
    if cond_to_analyze == 'both':
        cond_select = np.arange(0, 4)
    if cond_to_analyze == 'both':
        n_cond = 2
    else:
        n_cond = 1
    
    options = {
    'chan_names': ChanNames,
    'nchan': len(ChanNames),
    'freq_range': [0.5, 30],
    'fit_range' : fit_range,
    'sfreq': 256,
    'ntrials': 'all',
    'bands': bands,
    'param_labels': param_labels,
    'subs_psi': subs_psi,
    'subs_ket': subs_ket,
    'npsi': len(subs_psi),
    'nket': len(subs_ket),
    'nsubs': Nsubs,
    'mode_list': mode_list,
    'wght_th': 2.5,
    'nan_policy' : 'impute_normal',
    'impute_param' : 'group_mean',
    'interact_mode':'diff',
    'overwrite_fooof': 1,
    'overwrite_hyper_param_tune': 0,
    'overwrite_pls_on_params': 1,
    'overwrite_mmn_analysis': 0,
    'overwrite_bpls':1,
    'params_to_remove' : params_to_remove,
    'fit_domain': fit_domain,  # or 'low', 'high'
    'mpls_mode': 'separate',  # or 'together'
    'mc': 0,
    'cond_to_analyze': cond_to_analyze
}

    
    return options


def extract_params_fg(fgs, options):
    bands = options['bands']
    N = options['nsubs']
    nan_policy = options['nan_policy']
    
    n_channels = 64
    aper_param = np.zeros((N,n_channels,2))
    labels = bands.labels
    
    per_param = {}
    fit_metric = np.zeros((N,n_channels,4))# 1.r^2, 2.error, 3. n_peaks
    percent_non_nan_chans = {}
    detectable_count = {band: np.zeros(n_channels) for band in bands.labels}
    # Init
    for i in range(bands.n_bands):
        per_param[bands.labels[i]]  = np.zeros((N,n_channels,3))
        percent_non_nan_chans[bands.labels[i]] = np.zeros((N))
    # Main loop    
    for ind, fg in enumerate(fgs):
        aper_param[ind,:,:] = fg.get_params('aperiodic_params')# offset, exponent
        fit_metric[ind,:,0] = fg.get_params('error')
        fit_metric[ind,:,1] = fg.get_params('r_squared')
        fit_metric[ind,:,2] = fg.n_peaks_
        
        for band_label in labels:
            periodic_temp = get_band_peak_fg(fg, bands[band_label])
            valid_channels = np.isfinite(periodic_temp)

            percent_non_nan_chans[band_label][ind] = 100 * np.sum(valid_channels[:,0]) / n_channels
            
            detectable_count[band_label] += valid_channels[:,0].astype(int)
            if nan_policy == 'zero':
                per_param[band_label][ind,:,:] = check_nans(periodic_temp,
                                                  band_label,
                                                  options)
            elif nan_policy == 'impute_normal':
                if  options['impute_param'] == 'fixed':
                    per_param[band_label][ind,:,:] = check_nans(periodic_temp,
                                                      band_label,
                                                      options, impute_mean_pw=0.08,
                                                      impute_std_pw=0.02)
                if  options['impute_param'] == 'group_mean':
                    per_param[band_label][ind,:,:] = periodic_temp
    # for band_label in per_param.keys():
    #     values = per_param[band_label]
    #     valid_channels = np.isfinite(values[:,:,1])
    #     mean_non_nan   = np.nanmean(values[valid_channels,:],axis=0)   
    #     std_non_nan    = np.nanstd (values[valid_channels,:],axis=0)   
    #     for ind in range(N):
    #         if any(std_non_nan<0):
    #             std_non_nan[std_non_nan<0] = 0
    #         per_param[band_label][ind,:,:] = check_nans(values[ind,:,:],
    #                                         band_label=band_label,
    #                                         options=options, 
    #                                         impute_mean_pw=mean_non_nan[1]-2*std_non_nan[1] ,
    #                                         impute_std_pw=std_non_nan[1],
    #                                         cf_mean=mean_non_nan[0],
    #                                         cf_std=std_non_nan[0],
    #                                         bw_mean=mean_non_nan[2],
    #                                         bw_std=std_non_nan[2])
                           
    detectable_percent = {
        band: 100 * detectable_count[band] / N
        for band in labels}
    # plot_impute_dist( impute_mean=mean_non_nan-2*std_non_nan, impute_std=0.02, threshold=0.3)

    return aper_param, per_param, fit_metric, percent_non_nan_chans, detectable_percent
# def nan_stats():
    
from scipy.stats import truncnorm


def check_nans(data, band_label, options, impute_mean_pw=0.08, impute_std_pw=0.02, 
               threshold=0.3, min_peak_width=1.0, max_peak_width=6.0, 
               bw_mean=3.5, bw_std=0.5, cf_mean=None, cf_std=1):
    """
    Replace NaN rows in (n_channels, 3) array based on the specified policy.
    
    Assumes that if one value in a row is NaN, all are NaN.
    
    Parameters
    ----------
    data : np.ndarray
        2D array of shape (n_channels, 3): [center freq, power, bandwidth].
    nan_policy : str
        One of ['zero', 'mean', 'impute', 'uniform'].
    impute_mean : float
        Mean of the truncated normal distribution for power imputation.
    impute_std : float
        Standard deviation of the truncated normal distribution for power imputation.
    threshold : float
        Upper bound for truncated normal or uniform sampling.
        
    Returns
    -------
    data : np.ndarray
        Array with NaNs replaced according to the policy.
    stat : float
        Proportion of NaN rows.
    """
    
    mean_band_features = {'delta': [2.5,  3],
                          'theta': [6,    4],
                          'alpha': [10.5, 5],
                          'beta' : [21.5, 17]}
    nan_policy = options['nan_policy']
    bands = options['bands']
    # Identify rows with any NaNs (assumes full-row NaNs)
    nan_rows = np.isnan(data).any(axis=1)

    if nan_policy == 'zero':
        data[nan_rows] = [0, 0, 0]

    # elif nan_policy == 'mean':
    #     mean_freq = np.nanmean(data[:, 0])
    #     mean_power = np.nanmean(data[:, 1])
    #     mean_bw = np.nanmean(data[:, 2])
    #     data[nan_rows] = [mean_freq, mean_power, mean_bw]

    elif nan_policy == 'impute_normal':
        a = (0 - impute_mean_pw) / impute_std_pw   # lower bound (at 0)
        b = np.inf                                 # no upper bound
        imputed_power = truncnorm.rvs(
            a, b, loc=impute_mean_pw, scale=impute_std_pw, size=np.sum(nan_rows))
        
        
        if options['impute_param'] == 'fixed':
            bw_mean = (min_peak_width + max_peak_width) / 2
            bw_std = 0.5  
            cf_mean = np.nanmean(data[:, 0])

        a_bw = (min_peak_width - bw_mean) / bw_std
        b_bw = (max_peak_width - bw_mean) / bw_std
        imputed_bw = truncnorm.rvs(a_bw, b_bw, loc=bw_mean, scale=bw_std, 
                                   size=np.sum(nan_rows))
        
        a_cf = (1+bands.bands[band_label][0] - cf_mean) / cf_std
        b_cf = (bands.bands[band_label][1]-1 - cf_mean) / cf_std
        imputed_cf = truncnorm.rvs(a_cf, b_cf, loc=cf_mean, scale=cf_std, 
                                   size=np.sum(nan_rows))
        # mean_bw = np.nanmean(data[:, 2])
        # if np.isnan(mean_freq):
        #     mean_freq = mean_band_features[band_label][0]
        # if np.isnan(mean_bw):
        #     mean_bw = mean_band_features[band_label][1]
        # imputed_rows = np.column_stack([
        #     np.full(np.sum(nan_rows), mean_freq),
        #     imputed_power,
        #     imputed_bw])
        imputed_rows = np.column_stack([
            imputed_cf,
            imputed_power,
            imputed_bw])
        data[nan_rows] = imputed_rows

    else:
        raise ValueError("Nan policy not understood. Choose from ['zero', 'mean', 'impute'].")

    return data


def plot_impute_dist(bands, paths, impute_mean=0.25, impute_std=0.05, threshold=[-np.inf, np.inf], 
                     min_peak_width=1, max_peak_width=6,label=None, filename=None):

    fig_path = paths['results'] + '/figs/'

    if label[-3::] == '_PW':
        a = (0 - impute_mean) / impute_std   # lower bound (at 0)
        b = np.inf                                 # no upper bound
        samples = truncnorm.rvs(
            a, b, loc=impute_mean, scale=impute_std, size=10000)
        threshold = [0, b]
    if label[-3::] == '_BW':
        a_bw = (min_peak_width - impute_mean) / impute_std
        b_bw = (max_peak_width - impute_mean) / impute_std
        samples = truncnorm.rvs(a_bw, b_bw, loc=impute_mean, scale=impute_std, 
                                   size=10000)
        threshold = [min_peak_width, max_peak_width]
    if label[-3::] == '_CF':
        a_cf = (1+bands.bands[label[:-3]][0] - impute_mean) / impute_std
        b_cf = (bands.bands[label[:-3]][1]-1 - impute_mean) / impute_std
        samples = truncnorm.rvs(a_cf, b_cf, loc=impute_mean, scale=impute_std, 
                                   size=10000)
        threshold = [1+bands.bands[label[:-3]][0], bands.bands[label[:-3]][1]-1]
    # Plot the histogram
    plt.figure(figsize=(8, 5))
    plt.hist(samples, bins=50, density=True, alpha=0.7, edgecolor='k')
    plt.axvline(impute_mean, color='red', linestyle='--', label=f'Mean ({impute_mean:.2f})')
    plt.axvline(threshold[0], color='blue', linestyle='--', label=f'Lower Threshold ({str(threshold[0])})')
    plt.axvline(threshold[1], color='green', linestyle='--', label=f'Upper Threshold ({str(threshold[1])})')

    plt.title(f'Truncated Normal Distribution for Imputation in {label}')
    # plt.xlabel('Power Value')
    # plt.ylabel('Density')
    # plt.xticks(fontsize=20)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f'{fig_path}{filename}.png', dpi=300)
    plt.show()

# def check_nans(data, nan_policy):
#     """Check an array for nan values, and replace, based on policy."""

#     # Find where there are nan values in the data
#     nan_inds = np.where(np.isnan(data))
#     stat = len(nan_inds[0])/3/len(data[:, 1])
#     # Apply desired nan policy to data
#     if nan_policy == 'zero':
#         data[nan_inds] = 0
#     elif nan_policy == 'mean':
#         data[nan_inds] = np.nanmean(data)
#     # elif nan_policy == 'impute':
        
#     else:
#         raise ValueError('Nan policy not understood.')

#     return data, stat
# def adjacent_values(vals, q1, q3):
#     upper_adjacent_value = q3 + (q3 - q1) * 1.5
#     upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

#     lower_adjacent_value = q1 - (q3 - q1) * 1.5
#     lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
#     return lower_adjacent_value, upper_adjacent_value

