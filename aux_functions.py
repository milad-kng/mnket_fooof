# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 20:35:00 2024

@author: milad
"""
import numpy as np
import scipy.io as sio
import pickle
# from   os.path import isfile
# import os
import pandas as pd
from   collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib
from   matplotlib import cm, colors, colorbar
import matplotlib.patches as patches
import seaborn as sns
from matplotlib.colors import TwoSlopeNorm
# Import the FOOOF object
from fooof import FOOOF
from fooof import FOOOFGroup
from fooof.objs import fit_fooof_3d, combine_fooofs
from fooof.bands import Bands
from fooof.analysis import get_band_peak_fg,get_band_peak_fm 
# Import plotting function for model parameters and components
from fooof.plts.periodic import plot_peak_fits, plot_peak_params
from fooof.plts.aperiodic import plot_aperiodic_params, plot_aperiodic_fits
from fooof.plts.annotate import plot_annotated_model
# MNE
import mne
from mne.viz import plot_topomap
from mne.time_frequency import psd_welch
# Import Pyls
from pyls import behavioral_pls
from pyls import meancentered_pls
# from pyls import pls_regression
# import modeling objects
from sklearn.linear_model import LinearRegression
# from sklearn.cross_decomposition import PLSRegression
# import scipy
from scipy import stats
# import statsmodels
# class options():
#     def __init__(self,):

def extract_params_fg(fgs,bands,N):
    aper_param = np.zeros((N,64,2))
    labels = bands.labels
    per_param = {}
    fit_metric = np.zeros((N,64,2))
    for i in range(bands.n_bands):
        per_param[bands.labels[i]]  = np.zeros((N,64,3))
    for ind, fg in enumerate(fgs):
        aper_param[ind,:,:] = fg.get_params('aperiodic_params')# offset, exponent
        fit_metric[ind,:,0] = fg.get_params('error')
        fit_metric[ind,:,1] = fg.get_params('r_squared')
        for band_label in labels:
            per_param[band_label][ind,:,:] = get_band_peak_fg(fg, bands[band_label])
            per_param[band_label][ind,:,:] = check_nans(per_param[band_label][ind,:,:],
                                                        nan_policy='zero')
    return aper_param, per_param, fit_metric

def check_nans(data, nan_policy):
    """Check an array for nan values, and replace, based on policy."""

    # Find where there are nan values in the data
    nan_inds = np.where(np.isnan(data))

    # Apply desired nan policy to data
    if nan_policy == 'zero':
        data[nan_inds] = 0
    elif nan_policy == 'mean':
        data[nan_inds] = np.nanmean(data)
    else:
        raise ValueError('Nan policy not understood.')

    return data
# def adjacent_values(vals, q1, q3):
#     upper_adjacent_value = q3 + (q3 - q1) * 1.5
#     upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

#     lower_adjacent_value = q1 - (q3 - q1) * 1.5
#     lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
#     return lower_adjacent_value, upper_adjacent_value

import matplotlib.transforms

def plot_pls_analysis(
    mpls_fooof_params, 
    fig_path, 
    nTrials, 
    drug, 
):
    # Variables
    mc = mpls_fooof_params.inputs.mean_centering 
    wght_th = 2.5
    ChanNames = ['Fp1', 'AF7', 'AF3', 'F1' ,'F3' ,'F5', 'F7' ,'FT7', 'FC5', 'FC3', 'FC1', 'C1',\
                 'C3' ,'C5' ,'T7' ,'TP7' ,'CP5' ,'CP3' ,'CP1' ,'P1' ,'P3' ,'P5' ,'P7', 'P9',\
                 'PO7', 'PO3' ,'O1', 'Iz', 'Oz', 'POz' ,'Pz', 'CPz' ,'Fpz' ,'Fp2', 'AF8',\
                 'AF4' ,'AFz' ,'Fz' ,'F2' ,'F4' ,'F6', 'F8' ,'FT8' ,'FC6', 'FC4' ,'FC2',\
                 'FCz' ,'Cz', 'C2', 'C4', 'C6' ,'T8', 'TP8', 'CP6', 'CP4', 'CP2', 'P2', 'P4',\
                 'P6' ,'P8' ,'P10' ,'PO8' ,'PO4' ,'O2']

    Nchan      = len(ChanNames)
    sfreq = 256
    param_labels = ['Offset','Exponent','Delta_CF','Delta_PW','Delta_BW',\
                    'Theta_CF','Theta_PW','Theta_BW','Alpha_CF',\
                    'Alpha_PW','Alpha_BW','Beta_CF','Beta_PW','Beta_BW']
    # Compute weights
    Nparam = int(np.size(mpls_fooof_params.x_weights, axis=0) / Nchan)
    x_w = np.copy(np.reshape(mpls_fooof_params.bootres.x_weights_normed[:, 0], (Nchan, Nparam)))

    # Check for NaN weights
    if not (np.argwhere(np.isnan(x_w))).size == 0:
        print("NaN Weights exist")
        check_nans(x_w, "zero")

    # Convert weights to DataFrame
    x_w_df = pd.DataFrame(np.transpose(x_w), index=param_labels, columns=ChanNames)
    x_w_df[np.abs(x_w_df) < wght_th] = 0

    # Plot heatmap of weights
    fig, ax = plt.subplots(1, 1, figsize=(30, 10))
    sns.heatmap(x_w_df, cmap="bwr", vmax=5, vmin=-5, linewidths=1, linecolor="black", ax=ax)

    # Customize heatmap appearance
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=25)

    plt.xticks(ticks=range(len(ChanNames)), labels=ChanNames, fontsize=18, rotation="vertical")
    plt.yticks(ticks=range(len(param_labels)), labels=param_labels, fontsize=25)

    plt.title("PLS on combined (Chan x Parameters) Features - X_weights_normed are shown")
    plt.xlabel("Channels", {"fontsize": 25})
    plt.ylabel("Parameters", {"fontsize": 25})

    # Adjust tick label positions
    dx, dy = 5 / 30, 0
    offset = matplotlib.transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
    for label in ax.xaxis.get_majorticklabels():
        label.set_transform(label.get_transform() + offset)

    dx, dy = 0, -8 / 30
    offset = matplotlib.transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
    for label in ax.yaxis.get_majorticklabels():
        label.set_transform(label.get_transform() + offset)

    plt.savefig(
        f"{fig_path}/heatmap_fooof_meancentre_{mc}_{nTrials}_{drug}.png",
        bbox_inches="tight",
        dpi=600,
    )

    # Plot contrasts
    idx_sig = np.squeeze(np.where(mpls_fooof_params.permres.pvals < 0.05))
    mpls_fooof = {
        "LV1": mpls_fooof_params.bootres.contrast[:, idx_sig],
        "UL": mpls_fooof_params.bootres.contrast_ci[:, idx_sig, 1],
        "LL": mpls_fooof_params.bootres.contrast_ci[:, idx_sig, 0],
        "Condition": np.array(["std", "dev", "std", "dev"]),
        "Group": np.array(["pla", "pla", "ket", "ket"] if drug == "ket" else ["pla", "pla", "psi", "psi"]),
    }
    mpls_fooof = pd.DataFrame(mpls_fooof)
    
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    sns.barplot(
        data=mpls_fooof,
        x="Condition" if mc else "Group",
        y="LV1",
        hue="Group" if mc else "Condition",
        palette=["#2986cc", "#D43F3A"],
        ax=ax,
    )
    ax.set_xlabel("Condition" if mc else "Group", fontsize=20)
    ax.errorbar(
        x=np.sort([r.get_x() for r in ax.patches]) + (np.diff(np.sort([r.get_x() for r in ax.patches])).min() / 2),
        y=mpls_fooof["LV1"],
        fmt="none",
        yerr=mpls_fooof["UL"].values - mpls_fooof["LV1"].values,
        uplims=True,
        ecolor="black",
    )
    
    plt.savefig(
        f"{fig_path}/LV_{mc}_{nTrials}_{drug}.png",
        bbox_inches="tight",
        dpi=400,
    )
    
    # Plot topomaps
    cnorm = TwoSlopeNorm(vmin=x_w.min(), vcenter=0, vmax=x_w.max())
    for i, param in enumerate(param_labels):
        X = x_w[:, i]
        montage = mne.channels.make_standard_montage("biosemi64")
        info = mne.create_info(montage.ch_names, sfreq, "eeg").set_montage(montage)
        fig, ax1 = plt.subplots(ncols=1)
        im, cn = plot_topomap(
            X, pos=info, axes=ax1, show=False, cmap="RdBu_r", ch_type="eeg", contours=0, size=5, cnorm=cnorm
        )
        ax1.set_title(f"BSRs for {param}", {"fontsize": 14})
        plt.savefig(
            f"{fig_path}/BSRs_{param}_{mc}_{nTrials}_{drug}.png",
            dpi=500,
        )



from scipy.io import loadmat

def plot_group_comparison(mpls_fooof_params_interact, fig_path, fontsize=20):
    """
    Plots barplots with error bars for group comparison.
    """
    # Prepare data
    idx_sig = np.squeeze(np.where(mpls_fooof_params_interact.permres.pvals < 0.05))
    mpls_fooof = {
        'LV1': mpls_fooof_params_interact.bootres.contrast[:, idx_sig],
        'UL': mpls_fooof_params_interact.bootres.contrast_ci[:, idx_sig, 1],
        'LL': mpls_fooof_params_interact.bootres.contrast_ci[:, idx_sig, 0],
        'Condition': np.array(['std', 'dev', 'std', 'dev', 'std', 'dev', 'std', 'dev']),
        'Group': np.array(['pla1', 'pla1', 'ket', 'ket', 'pla2', 'pla2', 'psi', 'psi']),
    }
    mpls_fooof = pd.DataFrame(mpls_fooof)

    # Barplot
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    sns.barplot(data=mpls_fooof, x='Group', y="LV1", hue='Condition', palette=["#2986cc", "#D43F3A"])
    
    # Add error bars
    x = [r.get_x() for r in ax.patches]
    nx = np.sort(x)
    abs_err = mpls_fooof['UL'].values - mpls_fooof['LV1'].values
    ax.errorbar(x=nx + (np.diff(nx).min() / 2), y=mpls_fooof['LV1'], fmt='none', yerr=abs_err, lolims=True, ecolor='black')
    abs_err = -mpls_fooof['LL'].values + mpls_fooof['LV1'].values
    ax.errorbar(x=nx + (np.diff(nx).min() / 2), y=mpls_fooof['LV1'], fmt='none', yerr=abs_err, uplims=True, ecolor='black')

    # Customize plot
    ax.grid(False)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    ax.set_ylabel("Latent variable 1", fontsize=fontsize)
    ax.set_xlabel("Group", fontsize=fontsize)
    plt.setp(ax.get_legend().get_texts(), fontsize='12')
    plt.setp(ax.get_legend().get_title(), fontsize='12')
    ax.set_facecolor('none')

    # Save figure
    plt.savefig(f"{fig_path}/LV_3groups.png", bbox_inches='tight', dpi=600)
    plt.show()
def plot_example_psd(data_path, fig_path, fontsize=30):
    """
    Plots an example PSD fit with FOOOF.
    """
    freq_range = [0.5, 30]
    fm = FOOOF(min_peak_height=0.3, verbose=False, peak_width_limits=[1, 8], max_n_peaks=6)
    spec = loadmat(data_path + 'psd_std_psi_roving_trial_series_csd_all.mat')['psd_std_psi']
    X = spec[3, 11, :]
    freqs = np.linspace(freq_range[0], freq_range[1], np.shape(X)[-1])
    fm.fit(freqs, X)
    
    # Plot FOOOF model
    plot_annotated_model(fm, annotate_aperiodic=False)
    ax = plt.gca()
    ax.grid(False)
    ax.set_xlabel('Frequency (Hz)', fontsize=fontsize)
    ax.set_ylabel('log(Power)', fontsize=fontsize)
    ax.set_facecolor('none')
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.savefig(f"{fig_path}/example_PSD_FOOOF_chan_11.jpg", bbox_inches='tight', dpi=600)
    plt.close()

def plot_violin_distributions(param_list, param_labels, fig_path, mode_list, Nsubs, fontsize=20):
    """
    Plots violin distributions for specified parameters and channels.
    """
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))
    c = 0
    for i in [0, 1, 9, 12]:  # Example indices for channels
        ch_idx = {0: 28, 1: 28, 9: 30, 12: 11}[i]  # Channel mapping
        paras_df = np.array(param_list[:, :, ch_idx, i]).tolist()
        
        ax = axes[np.unravel_index(c, (2, 2))]
        parts = ax.violinplot(paras_df, showmedians=False, showmeans=False, showextrema=False)
        for pc in parts['bodies']:
            pc.set_facecolor('#D43F3A')
            pc.set_edgecolor('black')
            pc.set_alpha(1)
        
        quartile1, medians, quartile3 = np.percentile(paras_df, [25, 50, 75], axis=1)
        whiskers = np.array([
            np.clip(sorted_array, q1, q3)
            for sorted_array, q1, q3 in zip(paras_df, quartile1, quartile3)
        ])
        whiskers_min, whiskers_max = whiskers[:, 0], whiskers[:, 1]
        
        inds = np.arange(1, len(medians) + 1)
        ax.scatter(inds, medians, marker='o', color='white', s=30, zorder=3)
        ax.vlines(inds, quartile1, quartile3, color='k', linestyle='-', lw=5)
        ax.vlines(inds, whiskers_min, whiskers_max, color='k', linestyle='-', lw=1)
        ax.set_title(param_labels[i], fontsize=fontsize)
        ax.set_xticks(ticks=[1, 2, 3, 4])
        ax.set_xticklabels(labels=mode_list)
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)
        plt.xticks(fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        ax.grid(False)
        c += 1
    
    plt.tight_layout()
    plt.savefig(f"{fig_path}/distribution_visualization.png", dpi=600)
    plt.show()

from scipy.stats import pearsonr


def plot_behav_results(X_reg, Y_reg, fig_path, nTrials, 
                       param_to_do, bpls_single_std, bpls_single_dev, wght_th=0.5, 
                       DoPredict=True, DoPlot=True, sfreq=256, Nchan=64, ChanNames=None):
    front_ch =  ['Fp1', 'AF7', 'AF3', 'F1' ,'F3' ,'F5', 'F7',
                  'Fpz' ,'Fp2', 'AF8','AF4' ,'AFz' ,'Fz' ,'F2' ,'F4' ,'F6', 'F8']
    param_labels = ['Offset','Exponent','Delta_CF','Delta_PW','Delta_BW',\
                    'Theta_CF','Theta_PW','Theta_BW','Alpha_CF',\
                    'Alpha_PW','Alpha_BW','Beta_CF','Beta_PW','Beta_BW']
    subs_psi    = ['3621', '4415', '4418','4419', '4420', '4421', '4426','4332', '4433', \
                   '4460', '4476', '4502', '4515', '4518', '4591','4592']
    subs_ket    = ['4431', '4446', '4447', '4458', '4482', '4487', '4488', '4548', '4494',\
                   '4499', '4500', '4520', '4532', '4497', '4534', '4459', '4507','4422','4478']
    # Compute correlation
    r, p = pearsonr(X_reg, Y_reg)
    Y_reg = Y_reg[:, np.newaxis]
    X_reg = X_reg[:, np.newaxis]
    data = np.concatenate((X_reg, Y_reg), axis=1)
    data_df = pd.DataFrame(data, index=np.concatenate((subs_ket, subs_psi)), 
                           columns=['x_score', 'y_score'])
    
    # Scatter plot with regression line
    sns.lmplot(x='x_score', y='y_score', data=data_df, ci=95)
    ax = plt.gca()
    ax.set_title(f'ASC correlation with spectral params: r = {r:.2f}, p = {p:.4f}')
    ax.set_facecolor('None')
    plt.xticks(fontsize=20)
    plt.xlabel('x_scores', fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylabel('y_scores', fontsize=20)
    plt.savefig(fig_path + f'ASC_reg_chan_impContCog_{nTrials}_.png', bbox_inches='tight', dpi=500)
    plt.show()

    # Compute PLS weights
    bpls_weights_std = [bpls_single_std[key].y_weights * bpls_single_std[key].bootres.x_weights_normed
                        for key in bpls_single_std if bpls_single_std[key].permres.pvals < 0.05]
    bpls_weights_dev = [bpls_single_dev[key].y_weights * bpls_single_dev[key].bootres.x_weights_normed
                        for key in bpls_single_dev if bpls_single_dev[key].permres.pvals < 0.05]
    
    for weights in [bpls_weights_std, bpls_weights_dev]:
        if weights:
            x_w = np.reshape(weights, (len(front_ch), len(param_to_do))).copy()
            x_w_df = pd.DataFrame(np.transpose(x_w), index=[param_labels[i] for i in param_to_do],
                                  columns=front_ch)
            x_w_df[np.abs(x_w_df) < wght_th] = 0
            
            # Plot heatmap
            fig, ax = plt.subplots(1, 1, figsize=(30, 10))
            sns.heatmap(x_w_df, cmap='bwr', vmax=5, vmin=-5, linewidths=1, linecolor='black')
            cbar = ax.collections[0].colorbar
            cbar.ax.tick_params(labelsize=25)
            plt.xticks(ticks=list(range(len(x_w_df.columns))), labels=x_w_df.columns, fontsize=25)
            plt.yticks(fontsize=28, rotation='horizontal')
            plt.xlabel('Channels', fontsize=25)
            plt.ylabel('Parameters', fontsize=25)
            plt.savefig(fig_path + f'/heatmap_bpls_ASC_{nTrials}_{DoPredict}_comb_behav.png', dpi=600)
            plt.show()
            
            if DoPlot:
                x_w_expand = []
                i = 0
                for c in range(Nchan):
                    ch = ChanNames[c] if ChanNames else f'Chan{c}'
                    if ch in front_ch:
                        x_w_expand.append(x_w[i, :])
                        i += 1
                    else:
                        x_w_expand.append(np.zeros((len(param_to_do))))
                x_w_expand = np.reshape(np.stack(x_w_expand, axis=0), (Nchan, len(param_to_do)))
                
                cnorm = TwoSlopeNorm(vmin=x_w.min(), vcenter=0, vmax=max(0.1, x_w.max()))
                montage = mne.channels.make_standard_montage('biosemi64')
                info = mne.create_info(montage.ch_names, sfreq, 'eeg').set_montage(montage)

                for i, param in enumerate(param_to_do):
                    X = x_w_expand[:, i].copy()
                    X[np.abs(X) < wght_th] = 0
                    fig, ax1 = plt.subplots(ncols=1)
                    im, _ = plot_topomap(X, pos=info, axes=ax1, show=False, cmap="RdBu_r", ch_type='eeg', 
                                         contours=0, size=5, names=info.ch_names, sensors=False, cnorm=cnorm)
                    ax1.set_title(f'ASC BSRs for {param_labels[param]}', fontsize=14)
                    
                    # Add colorbar
                    cbar_ax = fig.add_axes([0.85, 0.1, 0.04, 0.9])
                    clb = fig.colorbar(im, cax=cbar_ax)
                    clb.ax.tick_params(labelsize=25)
                    plt.savefig(fig_path + f'/BSRs_ASC_{param_labels[param]}_{nTrials}_{DoPredict}_comb_behav.png', dpi=900)
                    plt.show()


