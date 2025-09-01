# -*- coding: utf-8 -*-
"""
Created on Fri Aug 22 13:34:06 2025

@author: milad
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

def qc_plot_detectable(stats, detectable_percent, options, paths, drug,
                       colors={ 'delta': '#9467bd',
                                'theta': '#1f77b4',
                                'alpha': '#ff7f0e',
                                'beta' : '#2ca02c'}, cond_to_analyze='both',
                       do_plot=True):
    fig_path = paths['results'] + '/figs/'
    bands = options['bands']
    fig, axs = plt.subplots(len(options['mode_list']), len(bands.labels),
                            figsize=(14, 10), sharey=True, squeeze=False)
    df, df_retention = plot_retention_summary(stats, detectable_percent,
                                              axs, options)
    # df['Retention'] = df['Retention'].round(2)
    
    if do_plot:
        plot_retention_boxplots(df, options, colors, output_file=f"{fig_path}{drug}_boxplot_retention.png")
        plot_group_barplot(df_retention, colors, output_file=f"{fig_path}{drug}_barplot_retention.png")
    
    return df, df_retention
def qc_filter(param_list, detectable_percent, fit_metric, options):
    
    retained_channels, common_retained = \
            compute_common_retained_channels(detectable_percent,
                                             options, threshold=50)
    
    if options['fit_domain'] == 'high':
        param_labels = ['Offset','Exponent','Beta_CF','Beta_PW','Beta_BW']
        param_indices_to_remove = []
    else:
        param_labels = options['param_labels']
        param_indices_to_remove = [param_labels.index(param) for param in options['params_to_remove']]#[2,3,4,5,6,7]
    param_list = np.delete(param_list, param_indices_to_remove, axis=3)
    param_labels_filtered = [item for i, item in enumerate(param_labels) if i not in param_indices_to_remove]
    
    # =========================================================================
    #     Remove bad channels for all conditions
    # =========================================================================
    good_channels = get_good_channels_by_r2(fit_metric, options)
    good_channel_indices = [options['chan_names'].index(ch) for ch in good_channels]
    # if exclude_bad_fits:
    #     print('excluding bad model fits...')
    #     param_list = param_list[:,:,good_channel_indices,:]

    # else:
    #     print('keeping all channels ...')
    #     good_channels = chan_names
    #     good_channel_indices = range(1, len(chan_names))
    print(f"Good channels in ALL conditions: {len(good_channels)}")

    return param_list, common_retained, param_labels_filtered, good_channel_indices

def plot_retention_summary(stats, detectable_percent, axs, options):
    data = []
    retention_summary = []
    
    bands = options['bands']
    for i, mode in enumerate(options['mode_list']):
        for j, band in enumerate(bands.labels):
            values = np.array(stats[mode][band])
            if np.any(values > 100):
                print(f"⚠️ {mode} - {band} has values > 100: {values[values > 100]}")
            for subj_val in values:
                data.append({'Subband': band, 'Condition': mode, 'Retention': subj_val})

            n_subjects = len(values)
            n_above_50 = np.sum(values > 50)
            percentage = 100 * n_above_50 / n_subjects
            retention_summary.append({'Condition': mode, 'Subband': band, 'PercentOver50': percentage})

            ax = axs[i, j]
            ax.bar(range(64), detectable_percent[mode][band])
            ax.set_ylim(0, 120)
            if i == len(options['mode_list']) - 1:
                ax.set_xlabel("Channels")
            if j == 0:
                ax.set_ylabel(f"{mode}")
            if i == 0:
                ax.set_title(f"{band}")
            ax.axhline(50, color='red', linestyle='--', linewidth=1)

    return pd.DataFrame(data), pd.DataFrame(retention_summary)

def plot_retention_boxplots(df, options, color_map, output_file=None):
    sns.set(style="whitegrid")
    fontsize = 17
    fig, axes = plt.subplots(1, 4, figsize=(15, 5), sharey=True)
    for ax, band in zip(axes, options['bands'].labels):
        sns.boxplot(data=df[df['Subband'] == band],
                    x='Condition',
                    y='Retention',
                    ax=ax,
                    order=options['mode_list'],
                    color=color_map[band])
        ax.set_ylim(0, 110)

        ax.set_title(band, fontsize=fontsize)
        ax.set_ylabel('Proportion of channels retained (%)' if band == 'delta' else '', fontsize=fontsize)
        ax.set_xlabel('')
        ax.tick_params(axis='x', rotation=45, labelsize=fontsize)
    fig.tight_layout()
    if output_file:
        fig.savefig(output_file, dpi=500)
    else:
        plt.show()

def plot_group_barplot(df_retention, color_map, output_file=None):
    plt.figure(figsize=(8, 6))
    fontsize = 16
    sns.barplot(data=df_retention,
                x='Condition',
                y='PercentOver50',
                hue='Subband',
                palette=color_map)
    plt.ylabel('Participants with >50% channels retained (%)', fontsize=fontsize)
    plt.ylim(0, 100)
    plt.xticks(rotation=45,fontsize=fontsize)
    plt.legend(title='')
    plt.tight_layout()
    if output_file:
        plt.savefig(output_file, dpi=500)
    else:
        plt.show()

def get_good_channels_by_r2(fit_metric, options, r2_threshold=0.9):
    """
    Return a dictionary of channels with R^2 above threshold for each condition and band.
    """
    good_channel_modes = []
    chan_names = np.array(options['chan_names'])
    for mode in options['mode_list']:
        r2_vals   = fit_metric[mode][:, :, 1]  
        median_r2 = np.median(r2_vals, axis=0)
    
        # Get list of good channels for this mode
        good_channels_mode_indices = median_r2 >= r2_threshold
        good_channel_modes.append(chan_names[good_channels_mode_indices])
    
    # Intersect across all modes, preserving chan_names order
    good_channels = [ch for ch in chan_names if all(ch in chans for chans in good_channel_modes)]
    
    print(f"Good channels in ALL conditions: {len(good_channels)}")
    
    return good_channels
def compute_common_retained_channels(detectable_percent, options, threshold=50):
    """
    Compute common retained channels based on alpha/beta detectability across modes,
    without using sets or helper functions.

    Returns:
        retained_channels: dict[mode] -> list of retained channel indices
        common_retained: dict -> dict with keys 'std', 'dev', and 'both'
                         each containing sorted lists of retained channels
    """
    retained_channels = {}
    common_retained = {}

    retained_channels = {}

    # Step 1: compute retained indices for each mode
    for mode in options['mode_list']:
        alpha_detect = detectable_percent[mode]['alpha']
        beta_detect  = detectable_percent[mode]['beta']
        
        alpha_idx = [i for i in range(64) if alpha_detect[i] >= threshold]
        beta_idx  = [i for i in range(64) if beta_detect[i] >= threshold]
        
        if len(alpha_idx) >= 15 and len(beta_idx) >= 15:
            # both valid → AND
            retained_idx = [i for i in alpha_idx if i in beta_idx]
        elif len(alpha_idx) >= 15:
            # only alpha valid
            retained_idx = alpha_idx
        elif len(beta_idx) >= 15:
            # only beta valid
            retained_idx = beta_idx
        else:
            # neither valid
            retained_idx = []
        
        retained_channels[mode] = retained_idx


    # Step 2: compute common channels across first two modes (std)
    std_modes = options['mode_list'][:2]
    std_common = []
    for ch in range(64):
        if all(ch in retained_channels[mode] for mode in std_modes):
            std_common.append(ch)

    # Step 3: compute common channels across last two modes (dev)
    dev_modes = options['mode_list'][-2:]
    dev_common = []
    for ch in range(64):
        if all(ch in retained_channels[mode] for mode in dev_modes):
            dev_common.append(ch)

    # Step 4: compute channels common to both std and dev
    both_common = [ch for ch in std_common if ch in dev_common]

    # Step 5: store result
    common_retained = {
        'std': sorted(std_common),
        'dev': sorted(dev_common),
        'both': sorted(both_common)
    }

    return retained_channels, common_retained

def filter_common_channels_between_drugs(param_list, good_channels,
                                         common_retained, options):
        
    # Step 1: Channels that are good in both drugs
    common_good_channels = [ch for ch in good_channels['ket']  if ch in good_channels['psi']]
    
    # Step 2: Channels that are detected in both drugs for the given condition
    conds = ['std', 'dev', 'both']
    common_detect_chans = {}
    for cond_to_analyze in conds:
        list1 = common_retained['ket'][cond_to_analyze]
        list2 = common_retained['psi'][cond_to_analyze]
        common_detect_chans[cond_to_analyze] = [ch for ch in list2 if ch in list1]
        # list1 +[ch for ch in list2 if ch not in list1]
    
    # Step 3: Channels that are both good and commonly detected
    # common_final_channels = [ch for ch in common_detect_chans if ch in common_good_channels]
    
    # Step 4: Index and slice param_list
    # final_param_list = {}
    # final_chan_names = np.array(common_final_channels)  # same for both drugs
    
    # for drug in ['ket', 'psi']:
    #     retained_list = list(common_retained[drug][options['cond_to_analyze']])
    #     chan_indices = [retained_list.index(ch) for ch in common_final_channels]
    #     final_param_list[drug] = param_list[drug][:,:,chan_indices,:]
    # print(f"Shared good channels across ketamine and psilocybin: {len(common_final_channels)}")

    return common_detect_chans, common_good_channels







