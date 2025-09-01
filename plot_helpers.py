# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21 16:56:49 2025

@author: milad
"""

import matplotlib.transforms
import numpy as np
# from   os.path import isfile
# import os

import pandas as pd
from scipy.io import loadmat
import matplotlib.pyplot as plt
import matplotlib
from   matplotlib import cm, colors, colorbar
from matplotlib.colors import TwoSlopeNorm

import seaborn as sns
# MNE
import mne
from mne.viz import plot_topomap
#FOOOOF
from fooof import FOOOF

from fooof.plts.annotate import plot_annotated_model


def plot_group_comparison(mpls_fooof_params, fig_path, 
                          study = 'interact', fontsize=20,
                          cond=[], filename=None):
    """
    Plots barplots with error bars for group comparison.
    """
    # Prepare data
    all_idx_sig = np.squeeze(np.where(mpls_fooof_params.permres.pvals < 0.06))
    all_idx_sig = np.atleast_1d(all_idx_sig)
    for idx_sig in all_idx_sig:
        mpls_fooof = {
            f'LV{str(idx_sig+1)}': mpls_fooof_params.bootres.contrast[:, idx_sig],
            'UL': mpls_fooof_params.bootres.contrast_ci[:, idx_sig, 1],
            'LL': mpls_fooof_params.bootres.contrast_ci[:, idx_sig, 0]}
        if study == 'interact':
            mpls_fooof.update({
                'Condition': np.array(['std', 'dev', 'std', 'dev', 'std', 'dev', 'std', 'dev']),
                'Group': np.array(['pla1', 'pla1', 'ket', 'ket', 'pla2', 'pla2', 'psi', 'psi']),
            })
        if study == 'interact_diff':
            mpls_fooof.update({
                'Condition': np.array(['std', 'dev', 'std', 'dev']),
                'Group': np.array(['diff_ket', 'diff_ket', 'diff_psi', 'diff_psi']),
            })
        if study == 'interact_std':
            mpls_fooof.update({
                'Condition': np.array(['std', 'std']),
                'Group': np.array(['diff_ket', 'diff_psi']),
            })
        if study == 'ket' or  study == 'ket_both':
            mpls_fooof.update({
                'Group': np.array(['pla1', 'ket', 'pla1', 'ket']),
                'Condition': np.array(['std', 'std', 'dev', 'dev']),
            })
        if study == 'ket_std':
            mpls_fooof.update({
                'Group': np.array(['pla1', 'ket']),
                'Condition': np.array(['std', 'std']),
            })
        if study == 'ket_dev':
            mpls_fooof.update({
                'Group': np.array(['pla1', 'ket']),
                'Condition': np.array(['dev', 'dev']),
            })
        if study == 'psi':
            mpls_fooof.update({
                
                'Group': np.array(['pla2', 'psi', 'pla2', 'psi']),
                'Condition': np.array(['std', 'std', 'dev', 'dev']),
                
            })
        if study == 'psi_std':
            mpls_fooof.update({
                
                'Group': np.array(['pla2', 'psi']),
                'Condition': np.array(['std', 'std']),
                
            })
        if study == 'psi_dev':
            mpls_fooof.update({
                
                'Group': np.array(['pla2', 'psi']),
                'Condition': np.array(['dev', 'dev']),
                
            })

        mpls_fooof = pd.DataFrame(mpls_fooof)
    
        # Barplot
        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        if study == 'interact' or study == 'interact_diff':
            sns.barplot(data=mpls_fooof, x='Group', y=f'LV{str(idx_sig+1)}', hue='Condition', 
                        palette=["#2986cc", "#D43F3A"])
        else:
            sns.barplot(data=mpls_fooof, x='Condition', y=f'LV{str(idx_sig+1)}', hue='Group', 
                        palette=["#2986cc", "#D43F3A"])

        
        # Add error bars
        
        x = [r.get_x() for r in ax.patches]
        nx = np.sort(x)
        abs_err = mpls_fooof['UL'].values - mpls_fooof[f'LV{str(idx_sig+1)}'].values
        ax.errorbar(x=nx + (np.diff(nx).min() / 2), y=mpls_fooof[f'LV{str(idx_sig+1)}'], fmt='none', yerr=abs_err, lolims=True, ecolor='black')
        abs_err = -mpls_fooof['LL'].values + mpls_fooof[f'LV{str(idx_sig+1)}'].values
        ax.errorbar(x=nx + (np.diff(nx).min() / 2), y=mpls_fooof[f'LV{str(idx_sig+1)}'], fmt='none', yerr=abs_err, uplims=True, ecolor='black')
    
        # Customize plot
        p_val = float(mpls_fooof_params.permres.pvals[idx_sig])
        ax.set_title(f'Permutation P-value: {p_val:.4f}')
        ax.grid(False)
        plt.xticks(fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        ax.set_ylabel(f'Latent variable {str(idx_sig+1)}', fontsize=fontsize)
        ax.set_xlabel("Group", fontsize=fontsize)
        plt.setp(ax.get_legend().get_texts(), fontsize='12')
        plt.setp(ax.get_legend().get_title(), fontsize='12')
        ax.set_facecolor('none')
    
        # Save figure
        plt.savefig(f"{fig_path}/LV_{study}_idxSig_{str(idx_sig)}_{filename}.png", bbox_inches='tight', dpi=600)
    plt.show()
def plot_example_psd(data_path, fig_path, fontsize=30):
    """
    Plots an example PSD fit with FOOOF.
    """
    freq_range = [0.5, 30]
    fm = FOOOF(min_peak_height=0.3, verbose=False, peak_width_limits=[1, 6], max_n_peaks=6)
    spec = loadmat(data_path + 'psd_std_psi_roving_trial_series_csd_all.mat')['psd_std_psi']
    # spec = loadmat(data_path + 'psd_std_pla2_roving_trial_series_csd_all.mat')['psd_std_pla2']

    X = spec[3, 11, :]
    # X = spec[15, 2, :]
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
    # plt.close()

def plot_violin_distributions(param_list, param_labels, fig_path, mode_list, Nsubs, fontsize=20):
    """
    Plots violin distributions for specified parameters and channels.
    """
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))
    c = 0
    for i in [0, 1, 3, 6]:  # Example indices for channels
        # ch_idx = {0: 28, 1: 28, 9: 30, 12: 11}[i]  # Channel mapping
        ch_idx = {0: 28, 1: 28, 3: 30, 6: 11}[i]  # Channel mapping

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

def plot_heatmap(
    mpls_fooof, 
    param_labels,
    fig_path, 
    channels, channel_indices, drug='' ,filename=None
):
    # Variables
    mc = mpls_fooof.inputs.mean_centering 
    wght_th = 2.5
    NChan   = 64
    if drug == 'ket':
        Nsubs = 19
    else:
        if drug == 'psi':
            Nsubs = 16
        else:
            Nsubs = 35
    ChanNames = ['Fp1', 'AF7', 'AF3', 'F1' ,'F3' ,'F5', 'F7' ,'FT7', 'FC5', 'FC3', 'FC1', 'C1',\
                  'C3' ,'C5' ,'T7' ,'TP7' ,'CP5' ,'CP3' ,'CP1' ,'P1' ,'P3' ,'P5' ,'P7', 'P9',\
                  'PO7', 'PO3' ,'O1', 'Iz', 'Oz', 'POz' ,'Pz', 'CPz' ,'Fpz' ,'Fp2', 'AF8',\
                  'AF4' ,'AFz' ,'Fz' ,'F2' ,'F4' ,'F6', 'F8' ,'FT8' ,'FC6', 'FC4' ,'FC2',\
                  'FCz' ,'Cz', 'C2', 'C4', 'C6' ,'T8', 'TP8', 'CP6', 'CP4', 'CP2', 'P2', 'P4',\
                  'P6' ,'P8' ,'P10' ,'PO8' ,'PO4' ,'O2']

    NGoodChan      = len(channels)
    # param_labels = ['Offset','Exponent','Delta_CF','Delta_PW','Delta_BW',\
    #                 'Theta_CF','Theta_PW','Theta_BW','Alpha_CF',\
    #                 'Alpha_PW','Alpha_BW','Beta_CF','Beta_PW','Beta_BW']
    # Compute weights
    Nparam = int(np.size(mpls_fooof.x_weights, axis=0) / NGoodChan)
    idx_sig = np.squeeze(np.where(mpls_fooof.permres.pvals < 0.07))
    idx_sig = np.atleast_1d(idx_sig)
    for idx in idx_sig:
        x_w = np.copy(np.reshape(mpls_fooof.bootres.x_weights_normed[:, idx],
                                 (NGoodChan, Nparam)))
    
        # Check for NaN weights
        # if not (np.argwhere(np.isnan(x_w))).size == 0:
        #     print("NaN Weights exist")
        #     check_nans(x_w, "zero")
    
        # Convert weights to DataFrame
        x_w_df = pd.DataFrame(np.transpose(x_w), index=param_labels, 
                              columns=channels)
        x_w_df[np.abs(x_w_df) < wght_th] = 0
    
        # Plot heatmap of weights
        fig, ax = plt.subplots(1, 1, figsize=(NGoodChan*0.6, len(param_labels)*0.6))
        sns.heatmap(x_w_df, cmap="bwr", vmax=5, vmin=-5, linewidths=1, 
                    linecolor="black", ax=ax)
    
        # Customize heatmap appearance
        cbar = ax.collections[0].colorbar
        cbar.ax.tick_params(labelsize=25)
    
        plt.xticks(ticks=range(len(channels)), labels=channels, fontsize=18, 
                   rotation="vertical")
        plt.yticks(ticks=range(len(param_labels)), labels=param_labels, 
                   fontsize=25, rotation='horizontal')
    
        plt.title(f'PLS on combined (Chan x Parameters) Features - X_weights_normed are shown for {drug}: LV{str(idx)}',
                  {"fontsize": 25})
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
            f"{fig_path}/heatmap_fooof_meancentre_{mc}_{drug}_{filename}.png",
            bbox_inches="tight",
            dpi=600,
        )
        
def plot_topo_data(param_list,
    mpls_fooof, param_labels,fig_path, 
    nTrials,Nsubs, detectable_percent, retained_channels_idx, label='', 
    mode_list=['mode1', 'mode2','mode3','mode4'], threshold=50): 
    
    
    # Plot topomaps
    # cnorm = TwoSlopeNorm(vmin=x_w.min(), vcenter=0, vmax=x_w.max())
    wght_th = 2.5
    NChan   = 64
    sfreq = 256
    ChanNames = ['Fp1', 'AF7', 'AF3', 'F1' ,'F3' ,'F5', 'F7' ,'FT7', 'FC5', 'FC3', 'FC1', 'C1',\
                  'C3' ,'C5' ,'T7' ,'TP7' ,'CP5' ,'CP3' ,'CP1' ,'P1' ,'P3' ,'P5' ,'P7', 'P9',\
                  'PO7', 'PO3' ,'O1', 'Iz', 'Oz', 'POz' ,'Pz', 'CPz' ,'Fpz' ,'Fp2', 'AF8',\
                  'AF4' ,'AFz' ,'Fz' ,'F2' ,'F4' ,'F6', 'F8' ,'FT8' ,'FC6', 'FC4' ,'FC2',\
                  'FCz' ,'Cz', 'C2', 'C4', 'C6' ,'T8', 'TP8', 'CP6', 'CP4', 'CP2', 'P2', 'P4',\
                  'P6' ,'P8' ,'P10' ,'PO8' ,'PO4' ,'O2']
    NRetChan      = len(retained_channels_idx)
    
    Nparam = len(param_labels)
    retained_channels = np.array(ChanNames)[retained_channels_idx]
    # param_list = np.reshape(mpls_fooof.inputs.X,(len(mode_list), Nsubs, NGoodChan,
    #                                              len(param_labels)))
    idx_sig = np.squeeze(np.where(mpls_fooof.permres.pvals < 0.05))
    if idx_sig == []:
        idx_sig = [0]
    idx_sig = np.atleast_1d(idx_sig)
    n_cons = len(mode_list)
    for idx in idx_sig:
        x_w = np.copy(np.reshape(mpls_fooof.bootres.x_weights_normed[:, idx],
                                 (NRetChan, Nparam)))
        x_w_df = pd.DataFrame(np.transpose(x_w), index=param_labels, 
                              columns=retained_channels)
        x_w_df[np.abs(x_w_df) < wght_th] = 0
        fig, axes = plt.subplots(Nparam, n_cons, figsize=(4*n_cons, 3*Nparam))  # flexible sizing
        if Nparam == 1 and n_cons == 1:
            axes = np.array([[axes]])
        elif Nparam == 1 or n_cons == 1:
            axes = axes.reshape(Nparam, n_cons)
        
        montage = mne.channels.make_standard_montage("biosemi64")
        info = mne.create_info(montage.ch_names, sfreq, "eeg").set_montage(montage)
        
        for i, param_name in enumerate(param_labels):
            x_w_topo = np.zeros(NChan)
            sig_channels_mask = np.zeros(NChan)
            sig_channels = np.abs(x_w[:, i]) > 2.5
            vmin = param_list[:, :, :, i].mean(axis=1).min()
            vmax = param_list[:, :, :, i].mean(axis=1).max()
            cnorm = TwoSlopeNorm(vmin=vmin, vcenter=0.5*(vmin+vmax), vmax=vmax)
        
            for j, con in enumerate(mode_list):
                if param_name.startswith(('delta', 'theta', 'alpha')):
                    subband_name = param_name[:5]
                    detect_mask = detectable_percent[con][subband_name] < threshold
                elif param_name.startswith('beta'):
                    subband_name = param_name[:4]
                    detect_mask = detectable_percent[con][subband_name] < threshold
                elif param_name in ['offset', 'exponent']:
                    detect_mask = np.ones(NChan)
                else:
                    detect_mask = np.zeros(NChan)
        
                sig_channels_mask[retained_channels_idx] = sig_channels
        
                # Fill topography
                x_w_topo = np.mean(np.squeeze(param_list[j, :, :, i]), axis=0)
        
                ax = axes[i, j]
                mask_params = dict(markersize=7, markerfacecolor="y")
        
                im, _ = plot_topomap(
                    x_w_topo, pos=info, axes=ax, show=False, cmap="plasma",
                    ch_type="eeg", size=10, cnorm=cnorm,
                    mask=sig_channels_mask, mask_params=mask_params,
                    sphere="auto"
                )
                ax.set_title(f"{param_name} | {con[1::]}", fontsize=20)
        
        fig.suptitle(f"Topomaps ( Trials={nTrials}, Study={label})", fontsize=14)
        plt.tight_layout(rect=[0, 0, 1, 0.95])  # leave space for suptitle
        plt.savefig(f"{fig_path}/Topomaps_{nTrials}_{label}.png", dpi=600)
        plt.show()
        # for i, param_name in enumerate(param_labels):
        #     x_w_topo = np.zeros(NChan)
        #     sig_channels_mask = np.zeros(NChan)
        #     sig_channels = np.abs(x_w[:,i]) > 2.5
        #     vmin = param_list[:,:,:,i].mean(axis=1).min()
        #     vmax = param_list[:,:,:,i].mean(axis=1).max()
        #     cnorm = TwoSlopeNorm(vmin=vmin, 
        #                          vcenter=0.5*(vmin+vmax), 
        #                          vmax=vmax)
            
        #     for i_con, con in enumerate(mode_list):
        #         if param_name[0:5] == 'delta' or param_name[0:5] == 'theta' or \
        #             param_name[0:5] == 'alpha':
        #             subband_name = param_name[0:5]
        #             # mask detectable
        #             detect_mask = detectable_percent[con][subband_name]<threshold
        #         else:
        #             if param_name[0:4] == 'beta':
        #                 subband_name = param_name[0:4]
        #                 # mask detectable
        #                 detect_mask = detectable_percent[con][subband_name]<threshold
        #             if param_name == 'offset' or param_name == 'exponent':
        #                 subband_name == param_name
        #                 # mask detectable
        #                 detect_mask = np.ones(NChan)
        #         sig_channels_mask[retained_channels_idx] = sig_channels 
        #         # x_w_topo[channel_indices] = np.mean(np.squeeze(param_list[i_con,:,:,i]),
        #         #                                     axis=0)
        #         x_w_topo = np.mean(np.squeeze(param_list[i_con,:,:,i]),
        #                                             axis=0)
        #         montage = mne.channels.make_standard_montage("biosemi64")
        #         info = mne.create_info(montage.ch_names, sfreq, "eeg").set_montage(montage)
        #         # layout = mne.find_layout(info)
        #         # pos = layout.pos[:, :2] 
        #         fig, ax1 = plt.subplots(ncols=1)
        #         # im, cn = plot_topomap(
        #         #     x_w_topo, pos=info, axes=ax1, show=False, cmap="RdBu_r", ch_type="eeg", contours=0, size=5, cnorm=cnorm
        #         # )
        #         mask_params = dict(markersize=7, markerfacecolor="y")
        #         im, cn = plot_topomap(
        #             x_w_topo, pos=info, axes=ax1, show=False, cmap="plasma", 
        #             ch_type="eeg", size=5, cnorm=cnorm, 
        #             mask=sig_channels_mask, mask_params = mask_params,
        #             # names=info.ch_names,
        #             sphere="auto"
        #         )
                
          
        #         # for ch_name in np.array(ChanNames)[sig_channels]:
        #         #     if ch_name in info.ch_names:
        #         #         idx = info.ch_names.index(ch_name)
        #         #         x, y = pos[idx]  # Correct indexing
        #         #         ax1.text(x, y, '*', fontsize=14, ha='center', va='center', color='black')
        #         # ax_x_start = 0.85
        #         # ax_x_width = 0.04
        #         # ax_y_start = 0.1
        #         # ax_y_height = 0.8
        #         # cbar_ax = fig.add_axes([ax_x_start, ax_y_start, ax_x_width, ax_y_height])
        #         # clb = fig.colorbar(im, cax=cbar_ax)
        #         ax1.set_title(f'{param_name} in mode {con} ', fontsize=14)
        #         # ax1.set_title('', fontsize=14)
        #         # clb.ax.set_title('',fontsize=16) # title on top of colorbar
        #         plt.savefig(
        #             f"{fig_path}/BSRs_{param_name}_{nTrials}_{label}_{con}.png",
        #             dpi=500)
        #         plt.close(fig)


def plot_behav_results(X_reg, Y_reg, fig_path, nTrials, 
                       bpls_single, wght_th=2.5, plot_orig_correlation = False,
                       X_psd=None, Y_behav=None,DoPredict=True, DoPlot=False,
                       sfreq=256, Nchan=64, ChanNames=None,
                       select_channels=None, param_labels=None, drug='both',
                       behav_label=None):
    Nparam = len(param_labels)
    # front_ch =  select_channels#['Fp1', 'AF7', 'AF3', 'F1' ,'F3' ,'F5', 'F7',
                 # 'Fpz' ,'Fp2', 'AF8','AF4' ,'AFz' ,'Fz' ,'F2' ,'F4' ,'F6', 'F8']
    # param_labels = ['Offset','Exponent',
    #                 #'Delta_CF','Delta_PW','Delta_BW',\
    #                 #'Theta_CF','Theta_PW','Theta_BW',
    #                 'Alpha_CF',\
    #                 'Alpha_PW','Alpha_BW','Beta_CF','Beta_PW','Beta_BW']
    subs_psi    = ['3621', '4415', '4418','4419', '4420', '4421', '4426','4332', '4433', \
                   '4460', '4476', '4502', '4515', '4518', '4591','4592']
    subs_ket    = ['4431', '4446', '4447', '4458', '4482', '4487', '4488', '4548', '4494',\
                   '4499', '4500', '4520', '4532', '4497', '4534', '4459', '4507','4422','4478']
    
    if plot_orig_correlation:
        # Compute and plot correlation
        if drug == 'both':
            subs_list = np.concatenate((subs_ket, subs_psi))
        elif drug == 'ket':
            subs_list = subs_ket
        elif drug == 'psi':
            subs_list = subs_psi
        plot_psd_behavior_correlation(X_psd, Y_behav, subs_list, 
                                          # xlabel='Exponent at F5', 
                                          xlabel='ALPHA CF at F1', 

                                          ylabel=f'{behav_label}',
                                          title='ASC correlation with spectral params',
                                          fig_path=fig_path, fname=None, nTrials=nTrials)

    # Compute and plot correlation between scores
    plot_psd_behavior_correlation(X_reg, Y_reg, subs_list, 
                                      xlabel='PSD scores', 
                                      ylabel=f'{behav_label}',
                                      title='ASC correlation with spectral params',
                                      fig_path=fig_path, fname=None, nTrials=nTrials)
 
    # Compute PLS weights
    
    if bpls_single.permres.pvals < 1:
        print(f"p-value: {bpls_single.permres.pvals[0]:.4f}")
        weights = bpls_single.y_weights * bpls_single.bootres.x_weights_normed
                             
        # bpls_weights_dev = [bpls_single_dev[key].y_weights * bpls_single_dev[key].bootres.x_weights_normed
        #                     for key in bpls_single_dev if bpls_single_dev[key].permres.pvals < 0.05]
        # bpls_weights_std = [bpls_single_std[key].bootres.x_weights_normed
        #                     for key in bpls_single_std if bpls_single_std[key].permres.pvals < 0.05]
        # bpls_weights_dev = [bpls_single_dev[key].bootres.x_weights_normed
        #                     for key in bpls_single_dev if bpls_single_dev[key].permres.pvals < 0.05]
        
        x_w = np.reshape(weights, (len(select_channels), Nparam)).copy()
        x_w_df = pd.DataFrame(np.transpose(x_w), index=param_labels,
                              columns=select_channels)
        x_w_df[np.abs(x_w_df) < wght_th] = 0
        
        # Plot heatmap
        fig, ax = plt.subplots(1, 1, figsize=(30, 10))
        sns.heatmap(x_w_df, cmap='bwr', vmax=5, vmin=-5, linewidths=1, 
                    linecolor='black')
        cbar = ax.collections[0].colorbar
        cbar.ax.tick_params(labelsize=25)
        plt.xticks(ticks=list(range(len(x_w_df.columns))), 
                   labels=x_w_df.columns, fontsize=25,
                   rotation="vertical")
        plt.yticks(fontsize=28, rotation='horizontal')
        plt.xlabel('Channels', fontsize=25)
        plt.ylabel('Parameters', fontsize=25)
        plt.title(f"{behav_label}, permute p-val:{bpls_single.permres.pvals[0]:.4f}", fontsize=25)
        # Adjust tick label positions
        dx, dy = 20 / 30, 0
        offset = matplotlib.transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
        for label in ax.xaxis.get_majorticklabels():
            label.set_transform(label.get_transform() + offset)
        plt.savefig(fig_path + f'/heatmap_bpls_ASC_{nTrials}_{DoPredict}_comb_behav.png', dpi=600)
        plt.show()

def plot_psd_behavior_correlation(X, Y, subs_list, 
                                  xlabel='', 
                                  ylabel='',
                                  title='',
                                  fig_path=None, fname=None, nTrials=None, 
                                  alpha=0.05, random_state=0, n_boot=2000):
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy.stats import pearsonr, spearmanr
    from numpy import arctanh, tanh
    from math import sqrt
    from scipy.stats import norm
    """
    Computes Pearson correlation between PSD values and behavioral scores,
    and plots a scatterplot with regression line.
    
    Parameters:
    - X: np.array of shape (n_subjects,)
    - Y: np.array of shape (n_subjects,)
    - subs_ket, subs_psi: subject ID arrays to use for index
    - xlabel, ylabel: axis labels
    - title_prefix: text for the title before showing r and p
    - fig_path: optional path to save the figure
    - fname: optional filename to save the figure
    - nTrials: optional trial count to append to filename
    """
    # Ensure correct shape
    X = X[:, np.newaxis]
    Y = Y[:, np.newaxis]
    n = len(X)
    r, p = pearsonr(X.ravel(), Y.ravel())
    # Fisher r-to-z CI (analytic)
    zr = arctanh(r)
    se = 1 / sqrt(n - 3)
    z_crit = norm.ppf(1 - alpha/2)
    lo, hi = zr - z_crit * se, zr + z_crit * se
    r_lo, r_hi = tanh([lo, hi])
   
    # Bootstrap CI
    rng = np.random.default_rng(random_state)
    r_boot = []
    rho_boot = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        xb, yb = X[idx].ravel(), Y[idx].ravel()
        r_boot.append(pearsonr(xb, yb)[0])
        rho_boot.append(spearmanr(xb, yb)[0])
    r_boot = np.array(r_boot)
    rho_boot = np.array(rho_boot)
    r_boot_lo, r_boot_hi = np.percentile(r_boot, [100*alpha/2, 100*(1-alpha/2)])
    rho_boot_lo, rho_boot_hi = np.percentile(rho_boot, [100*alpha/2, 100*(1-alpha/2)])
   
    # --------------------
    # Spearman correlation
    # --------------------
    rho, p_rho = spearmanr(X, Y)
   
    # Combine data into DataFrame
    data = np.column_stack((X, Y))
    data_df = pd.DataFrame(data, index=subs_list, 
                           columns=['x_score', 'y_score'])
   
    # Plot
    sns.lmplot(x='x_score', y='y_score', data=data_df, ci=95)
    ax = plt.gca()
    ax.set_title(
        f'{title}\n'
        f'Pearson r = {r:.2f}, 95% CI (Fisher) [{r_lo:.2f}, {r_hi:.2f}], '
        f'95% CI (Boot) [{r_boot_lo:.2f}, {r_boot_hi:.2f}], p = {p:.4f}\n'
        f'Spearman ρ = {rho:.2f}, 95% CI (Boot) [{rho_boot_lo:.2f}, {rho_boot_hi:.2f}], '
        f'p = {p_rho:.4f}',
        fontsize=11
    )
    ax.set_facecolor('None')
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlabel(xlabel, fontsize=20)
    plt.ylabel(ylabel, fontsize=20)
   
    # Save or show
    if fig_path and fname:
        if nTrials is not None:
            fname = f'{fname}_{nTrials}.png'
        plt.savefig(fig_path + fname, bbox_inches='tight', dpi=500)
    plt.show()
   
    return {
        "pearson_r": r, "pearson_p": p, "pearson_CI_fisher": (r_lo, r_hi),
        "pearson_CI_boot": (r_boot_lo, r_boot_hi),
        "spearman_rho": rho, "spearman_p": p_rho,
        "spearman_CI_boot": (rho_boot_lo, rho_boot_hi)
    }

def plot_violin_distributions(param_list, fig_path, param_labels, options):
    param_ket = param_list['ket']
    param_psi = param_list['psi']
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10,10))
    fontsize = 20
    c = 0
    for i in [0,1,3,6]:
        if i == 0 or i == 1:
            ch_idx = 26#28
        if i == 3:
            ch_idx = 26#30
        if i == 6:
            ch_idx = 47#11
        data = {'Group':   ['Standard']*(2*(19+16)) +\
                           ['Deviant' ]*(2*(19+16)),
                'Condition':  ['Placebo 1']*19 + ['Ketamine']*19 + ['Placebo 2']*16 + ['Psilocybin']*16 +\
                              ['Placebo 1']*19 + ['Ketamine']*19 + ['Placebo 2']*16 + ['Psilocybin']*16,
                'Values': np.concatenate([param_ket[0,:,ch_idx,i],
                                          param_ket[1,:,ch_idx,i],
                                          param_psi[0,:,ch_idx,i],
                                          param_psi[1,:,ch_idx,i],
                                          param_ket[2,:,ch_idx,i],
                                          param_ket[3,:,ch_idx,i],
                                          param_psi[2,:,ch_idx,i],
                                          param_psi[3,:,ch_idx,i]])}
         
        paras_df = pd.DataFrame(data)
        ax = axes[np.unravel_index(c,(2,2))]
         # sns.stripplot(x='Group', y='Values',  hue='Condition',data=paras_df, ax=ax,
         #               jitter=True, color='black', size=3)
        sns.violinplot(x='Group', y='Values', hue='Condition',data=paras_df,ax=ax,
                        palette=["#2986cc","#D43F3A",'#00bdff','#41b011'])
        
        sns.stripplot(x='Group', y='Values', hue='Condition', data=paras_df, 
                      color='black',size=5, jitter=True, ax=ax, dodge=True)
         # mode_list = ['_std_pla', '_std_ket', 'std_psi','_dev_pla', '_dev_ket','dev_psi']
         # ax.set_xticks(ticks=[1,2,3,4,5,6])
         # ax.set_xticklabels(labels=mode_list)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.legend_.remove()
    
        ax.set_title(param_labels[i]+' (channel: '+options['chan_names'][ch_idx]+')',
                     fontsize=fontsize)
        ax.tick_params(axis='x', labelsize=fontsize)  
        ax.tick_params(axis='y', labelsize=fontsize) 
        ax.grid(False)
        c = c+1
        ax.set_facecolor('none')
        # ax.patch.set_edgecolor('black')  
        # ax.patch.set_linewidth(0.5)  
     # Add spacing between subplots
    plt.tight_layout()
    
     # Show the plots
    # plt.savefig(fig_path+'others/distribution_visulaization2.png',dpi=600)
    # plt.show()
