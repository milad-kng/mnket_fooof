function [] = mnket_tf_group(options)
close all

subjs = options.subjects.all;

psd_std_pla = zeros(length(subjs),64,128);
psd_dev_pla = psd_std_pla;
psd_std_ket = psd_std_pla;
psd_dev_ket = psd_std_pla;

if strcmp(options.analysis,'MNKET')
    conds = {'placebo', 'ketamine'};
else
    conds = {'placebo', 'psilocybin'};
end

if options.tf.spectrum.overwrite 
    for cond = conds
        options.condition = char(cond);
        for subj_num = 1:length(subjs)
            fprintf('\nSubject#%d',subj_num)
            idCell = subjs(subj_num);
            id = char(idCell);
            [details, paths] = mn_subjects(id, options);
            load(details.PSD,"pwelch_psd")
            if (strcmp(cond,'placebo'))
                psd_std_pla(subj_num,:,:) = pwelch_psd{1}';
                psd_dev_pla(subj_num,:,:) = pwelch_psd{2}';
            else
                psd_std_ket(subj_num,:,:) = pwelch_psd{1}';
                psd_dev_ket(subj_num,:,:) = pwelch_psd{2}';
            end
        end
    end
    if strcmp(options.analysis,'MNKET')
        save(paths.psd_std_pla,"psd_std_pla")
        save(paths.psd_std_ket,"psd_std_ket")
        save(paths.psd_dev_pla,"psd_dev_pla")
        save(paths.psd_dev_ket,"psd_dev_ket")
        save(paths.subj_list,"subjs")
    else 
        if strcmp(options.analysis,'MNPSI')
            psd_std_pla2 = psd_std_pla;
            psd_dev_pla2 = psd_dev_pla;
            psd_std_psi  = psd_std_ket;
            psd_dev_psi  = psd_dev_ket;
            save(paths.psd_std_pla2,"psd_std_pla2")
            save(paths.psd_std_psi,"psd_std_psi")
            save(paths.psd_dev_pla2,"psd_dev_pla2")
            save(paths.psd_dev_psi,"psd_dev_psi")
            save(paths.subj_list,"subjs")
        end
    end
else
    [~, paths] = mn_subjects('0000', options);
    if strcmp(options.analysis,'MNKET')
        load(paths.psd_std_pla,"psd_std_pla")
        load(paths.psd_std_ket,"psd_std_ket")
        load(paths.psd_dev_pla,"psd_dev_pla")
        load(paths.psd_dev_ket,"psd_dev_ket")
        load(paths.subj_list,"subjs")
    else 
        if strcmp(options.analysis,'MNPSI')

            load(paths.psd_std_pla2,"psd_std_pla2")
            load(paths.psd_std_psi,"psd_std_psi")
            load(paths.psd_dev_pla2,"psd_dev_pla2")
            load(paths.psd_dev_psi,"psd_dev_psi")
            load(paths.subj_list,"subjs")
        end
    end
end

% N = size(all_pwr_dev_ket,1);   
if strcmp(options.analysis,'MNKET')

    for chan = 1:length(options.erp.channels)
        fprintf('\nchannel#%d',chan)
        freqs = linspace(options.preproc.highpassfreq,options.preproc.lowpassfreq,...
        size(psd_std_pla,3));
    
        h1 = figure;
        avg_tf_plot_per_chan(freqs,squeeze(psd_std_pla(:,chan,:)),...
            squeeze(psd_std_ket(:,chan,:)))
        title(strcat('placebo vs. ketamine average power spectra of STDs for channel',...
            options.erp.channels(chan)))
        legend('std-pla','std-pla-CI','std-pla-CI','std-ket','std-ket-CI','std-ket-CI')
        
        h2 = figure;
    
        avg_tf_plot_per_chan(freqs,squeeze(psd_dev_pla(:,chan,:)),...
            squeeze(psd_dev_ket(:,chan,:)))
        legend('dev-pla','dev-pla-CI','dev-pla-CI','dev-ket','dev-ket-CI','dev-ket-CI')
        title(strcat('placebo vs. ketamine average power spectra of DEVs for channel ',...
            options.erp.channels(chan)))
        if strcmp(options.erp.channels(chan),'Cz')
            % savefig(h1,paths.psd_std_fig)
            % savefig(h2,paths.psd_dev_fig)
        end
    end
end
% Combining PWPE signals
% mmn_pla_eps3 = zeros(length(subjs),64,129);
% mmn_ket_eps3 = mmn_pla_eps3;
mmn_high_pla = zeros(length(subjs),64,129);
mmn_low_pla  = mmn_high_pla;
mmn_high_drug = mmn_high_pla;
mmn_low_drug  = mmn_high_pla;
dest_loc = 'D:\Science\Coding\mnket_fooof_paper\mnket_fooof\erp_data\';
if options.tf.collectMMN
    for cond = conds
        options.condition = char(cond);
        for subj_num = 1:length(subjs)
            fprintf('\nSubject#%d',subj_num)
            idCell = subjs(subj_num);
            id = char(idCell);
            [details, ~] = mn_subjects(id, options);
            D = spm_eeg_load(details.erpfile);
            if (strcmp(cond,'placebo'))
                mmn_low_pla(subj_num,:,:)  = D(1:64,:,2);
                mmn_high_pla(subj_num,:,:) = D(1:64,:,3);
            else
                mmn_low_drug(subj_num,:,:)  = D(1:64,:,2);
                mmn_high_drug(subj_num,:,:) = D(1:64,:,3);
            end
        end
    end
    if strcmp(options.analysis,'MNKET')
        save(strcat(dest_loc,'mmn_low_pla1_' ,options.erp.type,'.mat'),"mmn_low_pla")
        save(strcat(dest_loc,'mmn_low_ket_' ,options.erp.type,'.mat'),"mmn_low_drug")
        save(strcat(dest_loc,'mmn_high_pla1_',options.erp.type,'.mat'),"mmn_high_pla")
        save(strcat(dest_loc,'mmn_high_ket_',options.erp.type,'.mat'),"mmn_high_drug")
    else
        if strcmp(options.analysis,'MNPSI')
            save(strcat(dest_loc,'mmn_low_pla2_' ,options.erp.type,'.mat'),"mmn_low_pla")
            save(strcat(dest_loc,'mmn_low_psi_' ,options.erp.type,'.mat'),"mmn_low_drug")
            save(strcat(dest_loc,'mmn_high_pla2_',options.erp.type,'.mat'),"mmn_high_pla")
            save(strcat(dest_loc,'mmn_high_psi_',options.erp.type,'.mat'),"mmn_high_drug")
        end
    end
end
end
% -------- Plot aux function --------- %
function [] = avg_tf_plot_per_chan(f,pwr_pla,pwr_ket)
f = log10(f);
pwr_pla = log10(pwr_pla);
pwr_ket = log10(pwr_ket);

N = size(pwr_pla,1);   
% base_diff = 
mean_pla = mean(pwr_pla,1);
mean_ket = mean(pwr_ket,1);

sem_pwr_pla = std(pwr_pla,0,1)/sqrt(N);
sem_pwr_ket = std(pwr_ket,0,1)/sqrt(N);

CI95 = tinv([0.025 0.975], N-1);          % Calculate 95% Probability Intervals Of t-Distribution

CI95_pla = bsxfun(@times, sem_pwr_pla, CI95(:));
CI95_ket = bsxfun(@times, sem_pwr_ket, CI95(:));

plot(f, mean_pla,'color','b')                      
hold on
plot(f, CI95_pla+mean_pla,':','Color','b')

plot(f, mean_ket,'color','r')                      
hold on
plot(f, CI95_ket+mean_ket,':','color','r')

hold off

ylabel('log-power')
xlabel('freq (Hz)')
grid
end

% Matrook
% for cond = {'placebo', 'ketamine'}
%     options.condition = char(cond);
%     for subj_num = 1:length(subjs)
%         fprintf('\nSubject#%d',subj_num)
%         idCell = subjs(subj_num);
% %         for chan = 1:64%length(options.erp.channels)
%             id = char(idCell);
%             [details, ~] = mnket_subjects(id, options);
%             load(details.CSD,"pwelch_csd")
% %             load(details.fullStd_pwr,"fft_welch_std_channels")
% %             load(details.fullDev_pwr,"fft_welch_dev_channels")
% %             fft_welch_std_channels = 10*log10(fft_welch_std_channels);
% %             fft_welch_dev_channels = 10*log10(fft_welch_dev_channels);
% %             if chan == 3
% %                 load(details.subbandsDev,"eeg_subbands_dev")
% %                 load(details.subbandsStd,"eeg_subbands_std")
% %                 eeg_subbands_std = eeg_subbands_std - eeg_subbands_std(:,1);
% %                 eeg_subbands_dev = eeg_subbands_dev - eeg_subbands_dev(:,1);
% %             end
%             if (strcmp(cond,'placebo'))
%                 all_pwr_std_pla(subj_num,chan,:) = pwelch_csd{1}(:,chan,chan)';
%                 all_pwr_dev_pla(subj_num,chan,:) = pwelch_csd{2}(:,chan,chan)';
% %                 all_pwr_dev_pla(subj_num,chan,:) = fft_welch_dev_channels(:,chan)';
% %                 all_pwr_dev_pla(subj_num,chan,:) = all_pwr_dev_pla(subj_num,chan,:) - ...
% %                     all_pwr_dev_pla(subj_num,chan,1);
% %                 if chan == 3
% %                     all_subbands_pla_dev(subj_num,:,:) = mean(eeg_subbands_dev.^2,2);
% %                     all_subbands_pla_std(subj_num,:,:) = mean(eeg_subbands_std.^2,2);
% %                 end
%             else
%                 all_pwr_std_ket(subj_num,chan,:) = pwelch_csd{1}(:,chan,chan)';
%                 all_pwr_dev_ket(subj_num,chan,:) = pwelch_csd{2}(:,chan,chan)';
% %                 all_pwr_std_ket(subj_num,chan,:) = fft_welch_std_channels(:,chan)';
% %                 all_pwr_std_ket(subj_num,chan,:) = all_pwr_std_ket(subj_num,chan,:) - ...
% %                     all_pwr_std_ket(subj_num,chan,1);
% %                 all_pwr_dev_ket(subj_num,chan,:) = fft_welch_dev_channels(:,chan)';
% %                 all_pwr_dev_ket(subj_num,chan,:) = all_pwr_dev_ket(subj_num,chan,:) - ...
% %                     all_pwr_dev_ket(subj_num,chan,1);
% %                 if chan == 3
% %                     all_subbands_ket_dev(subj_num,:,:) = mean(eeg_subbands_dev.^2,2);
% %                     all_subbands_ket_std(subj_num,:,:) = mean(eeg_subbands_std.^2,2);
% %                 end
%             end
% %         end
%     end
% end
