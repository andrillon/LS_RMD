%%
clear all;
close all;

%%
run ../LS_RMD_localdef.m
addpath((path_fieldtrip));
ft_defaults;
rmpath(genpath([path_fieldtrip filesep 'external' filesep 'signal' ]))
% 
addpath(genpath(path_LSCPtools));
addpath(genpath(path_IRASA))
files=dir([preproc_path filesep 'ICAcleaned_etrial_ft_*.mat']);


%%
res_mat=[];
redo=0; complete=0;

absThr=250;
nFc=0;
n_rlockerror=0; nE=0; rlockerrors=[];

for nF=1:length(files)
    file_name = files(nF).name;
    folder_name = files(nF).folder;
    SubID=file_name(1:end-4);
    seps=findstr(SubID,'ft_');
    SubID=SubID(seps(1)+3:end);
    tic;
    fprintf('... working on %s (%g/%g)\n',SubID,nF,length(files))
    
%     if redo==1 || exist([preproc_path filesep SubID '_TF_perTrial_varwin_ICAcleaned.mat'])==0 %Variable window
    if redo==1 || exist([preproc_path SubID '_TF_perTrial_ICAcleaned_RespLocked.mat'])==0 %Fixed window

        subj_pow=[];

        
        load([folder_name filesep file_name]);
        behav_table=readtable([folder_name filesep 'behav_' file_name(findstr(file_name,'ft_')+3:end-4) '.csv']);
        if size(behav_table,1)~=length(data.trial)
            error('different trial number in EEG and behav data')
        end
        
        % fix channel names
        mylabels=data.label;
        for nCh=1:length(mylabels)
            findspace=findstr(mylabels{nCh},' ');
            if isempty(findspace)
                newlabels{nCh}=mylabels{nCh};
            else
                if ismember(mylabels{nCh}(1),{'1','2','3','4','5','6','7','8','9'})
                    newlabels{nCh}=mylabels{nCh}(findspace+1:end);
                else
                    newlabels{nCh}=mylabels{nCh}(1:findspace-1);
                end
            end
        end
        
        fprintf('%2.0f-%2.0f\n',0,0)

        % remove NaN values
        cfg=[];
        cfg.trials=find(~isnan(behav_table.RT));
        data=ft_redefinetrial(cfg,data);
        
        % re-locking
        cfg=[];
        cfg.offset =-round(behav_table.RT(~isnan(behav_table.RT))*data.fsample);
        data_rlock=ft_redefinetrial(cfg,data);
        
        % crop trials
        cfg=[];
        cfg.toilim =[-1 1];
        data_rlock=ft_redefinetrial(cfg,data_rlock);
        
        
        % estimate the ERP
        tlck = ft_timelockanalysis([],data);
        figure; subplot(1,2,1); plot(tlck.time, tlck.avg);% legend(tlck.label);
        rlck = ft_timelockanalysis([],data_rlock);
        hold on; subplot(1,2,2); plot(rlck.time, rlck.avg);% legend(tlck.label);
        
        % subtract the ERP from the data
        data_minus_erp = data;
        data_rlock_minus_erp = data_rlock;
        for k = 1:numel(data.trial)
            data_minus_erp.trial{k} = data.trial{k} - tlck.avg(:,1:length(data.trial{k}));
            try
                data_rlock_minus_erp.trial{k} = data_rlock.trial{k} - rlck.avg(:,1:length(data_rlock.trial{k}));
            catch
                warning(['!!! Matrix dimensions do not agree - ' SubID])
                n_rlockerror=n_rlockerror+1;
                nE=nE+1;
                rlockerrors{nE}=SubID;
                continue
            end
        end
        
        % get the TF decomposition
        cfg              = [];
        cfg.output       = 'pow';
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.foi          = 2:0.5:30;                         % analysis 2 to 30 Hz in steps of 2 Hz
        cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
%         cfg.t_ftimwin    = 2./cfg.foi;   % 2 cycles per time window
        cfg.keeptrials   = 'yes';
        cfg.toi          = -1.0:0.05:3.5;                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
        TFRhann_induced = ft_freqanalysis(cfg, data_minus_erp);
        
        cfg.toi          = -1.0:0.05:1;                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
        TFRhann_rlock_induced = ft_freqanalysis(cfg, data_rlock_minus_erp);
        
        TFRhann_induced.powspctrm_log=log10(TFRhann_induced.powspctrm);
        TFRhann_rlock_induced.powspctrm_log=log10(TFRhann_rlock_induced.powspctrm);
        TFRhann_rlock_induced.powspctrm_bsl=TFRhann_rlock_induced.powspctrm_log-repmat(nanmean(TFRhann_induced.powspctrm_log(:,:,:,TFRhann_induced.time<0),4),[1 1 1 length(TFRhann_rlock_induced.time)]);
%         temp_powspctrm_bsl=squeeze(nanmean((TFRhann_induced.powspctrm_bsl(:,:,:,TFRhann_induced.time>-0.5 & TFRhann_induced.time<1.5)),1));
%         if size(temp_powspctrm_bsl,3)>39
%             temp_powspctrm_bsl=temp_powspctrm_bsl(:,:,2:end-1);
%             warning('More than expected time points - removing first and last');
%         end
        
        
        
        %         save([preproc_path filesep SubID '_TF_perTrial_varwin_ICAcleaned.mat'],'TFRhann','newlabels'); %Variable window
        save([preproc_path filesep SubID '_TF_perTrial_ICAcleaned_RespLocked_ERPremoved.mat'],'TFRhann_rlock_induced','newlabels'); %Fixed window
        save([preproc_path filesep SubID '_TF_perTrial_ICAcleaned_ERPremoved.mat'],'TFRhann_induced','newlabels'); %Fixed window
        
        %%%% evoked + induced
       cfg.toi          = -1.0:0.05:3.5;                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
        TFRhann = ft_freqanalysis(cfg, data);
        
        cfg.toi          = -1.0:0.05:1;                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
        TFRhann_rlock = ft_freqanalysis(cfg, data_rlock);
        
        TFRhann.powspctrm_log=log10(TFRhann.powspctrm);
        TFRhann_rlock.powspctrm_log=log10(TFRhann_rlock.powspctrm);
        TFRhann_rlock.powspctrm_bsl=TFRhann_rlock.powspctrm_log-repmat(nanmean(TFRhann.powspctrm_log(:,:,:,TFRhann.time<0),4),[1 1 1 length(TFRhann_rlock.time)]);
%         temp_powspctrm_bsl=squeeze(nanmean((TFRhann_induced.powspctrm_bsl(:,:,:,TFRhann_induced.time>-0.5 & TFRhann_induced.time<1.5)),1));

        
        %         save([preproc_path filesep SubID '_TF_perTrial_varwin_ICAcleaned.mat'],'TFRhann','newlabels'); %Variable window
        save([preproc_path filesep SubID '_TF_perTrial_ICAcleaned_RespLocked.mat'],'TFRhann_rlock','newlabels'); %Fixed window
        
    else
%         load([preproc_path filesep SubID '_TF_perTrial_varwin_ICAcleaned.mat']); %Variable window
%         load([preproc_path filesep SubID '_TF_perTrial_ICAcleaned_RespLocked.mat']); %Response locked, evoked
%         load([preproc_path filesep SubID '_TF_perTrial_ICAcleaned_RespLocked_ERPremoved.mat']);  %Induced, response locked
        load([preproc_path filesep SubID '_TF_perTrial_ICAcleaned_ERPremoved.mat']);  %Induced, stimulus locked


    end
%   

    nFc=nFc+1;
    
%     all_TFRhann(nFc,:,:,:)=squeeze(nanmean(TFRhann_rlock.powspctrm_bsl,1)); %Response locked, evoked
%     all_TFRhann(nFc,:,:,:)=squeeze(nanmean(TFRhann_rlock_induced.powspctrm_bsl,1)); %Induced, response locked
    all_TFRhann(nFc,:,:,:)=squeeze(nanmean(TFRhann_induced.powspctrm_bsl,1)); %Induced, stimulus locked


%     all_TFRhann(nFc,:,:,:)=squeeze(nanmean((TFRhann.powspctrm_bsl(:,:,:,TFRhann.time>-0.5 & TFRhann.time<1.5)),1));
%     TFtimes=TFRhann_induced.time(TFRhann_induced.time>-0.5 & TFRhann_induced.time<1.5);

    if length(SubID)==4 && SubID(1)=='A' % OLD (MONASH) - UP & DOWN 90%COH
        thisgroup=2;
        thisagegroup=1;
    elseif length(SubID)==7 % YOUNG (MONASH) - DOWN 50%COH
        thisgroup=1;
        thisagegroup=0;
    elseif length(SubID)==11 % YOUNG (TRINITY) - DOWN 50%COH
        thisgroup=3;
        thisagegroup=0;
    elseif length(SubID)==5 && SubID(3)=='8' % YOUNG - UP & DOWN 90%COH
        thisgroup=4;
        thisagegroup=0;
    elseif length(SubID)==5 && SubID(3)=='9'% OLD - UP & DOWN 90%COH
        thisgroup=5;
        thisagegroup=1;
    else
        thisgroup=NaN;
        thisagegroup=NaN;
    end
    all_group(nFc)=thisgroup; % STORE HERE THE INFO ABOUT THE GROUP
    all_agegroup(nFc)=thisagegroup; % STORE HERE THE INFO ABOUT THE AGE GROUP  
end

%%
%All 90% participants
% TFtimes=TFRhann_rlock.time; %Response locked, evoked
% TFtimes=TFRhann_rlock_induced.time; %Response locked, induced
TFtimes=TFRhann_induced.time; %Stimulus locked, induced

% faxis=TFRhann_rlock.freq; %Response locked, evoked
% faxis=TFRhann_rlock_induced.freq; %Response locked, induced
faxis=TFRhann_induced.freq; %Stimulus locked, induced


myLabels={'Fz','Cz','Pz','Oz'};
f1=figure;
maxAbsValues=[];
for nCh=1:length(myLabels)
    subplot(2,2,nCh);
    simpleTFplot(squeeze(nanmean(all_TFRhann((all_group==4|all_group==5),match_str(newlabels,myLabels{nCh}),:,:),1)),faxis,TFtimes,0,0);
    format_fig;
    title(myLabels{nCh})
    colorbar;
    maxAbsValues=[maxAbsValues max(max(abs(squeeze(nanmean(all_TFRhann(:,match_str(newlabels,myLabels{nCh}),:,:),1)))))];
end
for nCh=1:length(myLabels)
    subplot(2,2,nCh);
caxis([-1 1]*max(maxAbsValues))
end
%All younger participants
% figure;
% for nCh=1:length(myLabels)
%     subplot(2,2,nCh);
%     simpleTFplot(squeeze(nanmean(all_TFRhann((all_group==1|all_group==3|all_group==4),match_str(newlabels,myLabels{nCh}),:,:),1)),faxis,TFtimes,0,0);
%     format_fig;
%     title(myLabels{nCh})
% colorbar;
% end
% 
% %Younger 50% coherence
% maxAbsValues=[];
% 
% f1=figure;
% for nCh=1:length(myLabels)
%     subplot(2,2,nCh);
%     simpleTFplot(squeeze(nanmean(all_TFRhann((all_group==1|all_group==3),match_str(newlabels,myLabels{nCh}),:,:),1)),faxis,TFtimes,0,0);
%     format_fig;
%     title(myLabels{nCh})
% colorbar;
%     maxAbsValues=[maxAbsValues max(max(abs(squeeze(nanmean(all_TFRhann(:,match_str(newlabels,myLabels{nCh}),:,:),1)))))];
% end
% 
%Younger 90% coherence
f2=figure;
for nCh=1:length(myLabels)
    subplot(2,2,nCh);
    simpleTFplot(squeeze(nanmean(all_TFRhann((all_group==4),match_str(newlabels,myLabels{nCh}),:,:),1)),faxis,TFtimes,0,0);
    format_fig;
    title(myLabels{nCh})
colorbar;
    maxAbsValues=[maxAbsValues max(max(abs(squeeze(nanmean(all_TFRhann(:,match_str(newlabels,myLabels{nCh}),:,:),1)))))];
end

%Older 90% coherence
f3=figure;
for nCh=1:length(myLabels)
    subplot(2,2,nCh);
    simpleTFplot(squeeze(nanmean(all_TFRhann((all_group==5),match_str(newlabels,myLabels{nCh}),:,:),1)),faxis,TFtimes,0,0);
    format_fig;
    title(myLabels{nCh})
colorbar;
    maxAbsValues=[maxAbsValues max(max(abs(squeeze(nanmean(all_TFRhann(:,match_str(newlabels,myLabels{nCh}),:,:),1)))))];
end

for nCh=1:length(myLabels)
    figure(f1);
    subplot(2,2,nCh);
caxis([-1 1]*max(maxAbsValues))

figure(f2);
    subplot(2,2,nCh);
caxis([-1 1]*max(maxAbsValues))

figure(f3);
    subplot(2,2,nCh);
caxis([-1 1]*max(maxAbsValues))
end
%% Difference TF
% Young 90% - Old 90%
f1=figure;
for nCh=1:length(myLabels)
    subplot(2,2,nCh);
    TF_groupA=squeeze(nanmean(all_TFRhann((all_group==4),match_str(newlabels,myLabels{nCh}),:,:),1));
    TF_groupB=squeeze(nanmean(all_TFRhann((all_group==5),match_str(newlabels,myLabels{nCh}),:,:),1));
    simpleTFplot(TF_groupA-TF_groupB,faxis,TFtimes,0,0);
    format_fig;
    title(myLabels{nCh})
colorbar;
end

%Young 90% - Young 50%
% for nCh=1:length(myLabels)
%     subplot(2,2,nCh);
%     TF_groupA=squeeze(nanmean(all_TFRhann((all_group==4),match_str(newlabels,myLabels{nCh}),:,:),1));
%     TF_groupB=squeeze(nanmean(all_TFRhann((all_group==1|all_group==3),match_str(newlabels,myLabels{nCh}),:,:),1));
%     simpleTFplot(TF_groupA-TF_groupB,faxis,TFtimes,0,0);
%     format_fig;
%     title(myLabels{nCh})
% colorbar;
% end
% 
% %Young 50% - Old 90%
% for nCh=1:length(myLabels)
%     subplot(2,2,nCh);
%     TF_groupA=squeeze(nanmean(all_TFRhann((all_group==1|all_group==3),match_str(newlabels,myLabels{nCh}),:,:),1));
%     TF_groupB=squeeze(nanmean(all_TFRhann((all_group==5),match_str(newlabels,myLabels{nCh}),:,:),1));
%     simpleTFplot(TF_groupA-TF_groupB,faxis,TFtimes,0,0);
%     format_fig;
%     title(myLabels{nCh})
% colorbar;
% end

 %% Topographies 25Hz tag
cfg = [];
cfg.layout = 'acticap-64ch-standard2.mat';
cfg.channel=newlabels;
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);
correspCh=[];
for nCh=1:length(layout.label)-2
    correspCh(nCh)=match_str(newlabels,layout.label(nCh));
end

maxAbsValues=[];

f4=figure;
cmap=cbrewer('seq','YlOrRd',64); % select a sequential colorscale from yellow to red (64)
cmap(cmap<0)=0;
% faxis=TFRhann.freq
% faxis=faxis(faxis>=min(TFRhann.freq) & faxis<=max(TFRhann.freq));

% subplot(1,3,1); format_fig;
% % temp_topo=squeeze(nanmean(nanmean(nanmean(all_TFRhann(all_group==1|all_group==3,correspCh,faxis>2 & faxis<4,TFtimes>.200 & TFtimes<.500),1),3),4)); %Delta
% % temp_topo=squeeze(nanmean(nanmean(nanmean(all_TFRhann(all_group==1|all_group==3,correspCh,faxis>4 & faxis<8,TFtimes>.200 & TFtimes<.600),1),3),4)); %Theta
% % temp_topo=squeeze(nanmean(nanmean(nanmean(all_TFRhann(all_group==1|all_group==3,correspCh,faxis>8 & faxis<11,TFtimes>.500 & TFtimes<1.20),1),3),4)); %Alpha
% % temp_topo=squeeze(nanmean(nanmean(nanmean(all_TFRhann(all_group==1|all_group==3,correspCh,faxis>12 & faxis<16,TFtimes>.200 & TFtimes<1.0),1),3),4)); %Mu
% temp_topo=squeeze(nanmean(nanmean(nanmean(all_TFRhann(all_group==1|all_group==3,correspCh,faxis>16 & faxis<29,TFtimes>.200 & TFtimes<1.0),1),3),4)); %Beta
% simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
% colormap(cmap);
% title('Younger 50% Coherence','FontSize',10)
% maxAbsValues=[maxAbsValues max(max(abs(temp_topo)))];

subplot(1,2,1);
% temp_topo=squeeze(nanmean(nanmean(nanmean(all_TFRhann(all_group==4,correspCh,faxis>2 & faxis<4,TFtimes>.200 & TFtimes<.500),1),3),4)); %Delta
% temp_topo=squeeze(nanmean(nanmean(nanmean(all_TFRhann(all_group==4,correspCh,faxis>4 & faxis<8,TFtimes>.200 & TFtimes<.600),1),3),4)); %Theta
% temp_topo=squeeze(nanmean(nanmean(nanmean(all_TFRhann(all_group==4,correspCh,faxis>8 & faxis<11,TFtimes>.500 & TFtimes<1.20),1),3),4)); %Alpha
% temp_topo=squeeze(nanmean(nanmean(nanmean(all_TFRhann(all_group==4,correspCh,faxis>12 & faxis<16,TFtimes>.200 & TFtimes<1.0),1),3),4)); %Mu
temp_topo=squeeze(nanmean(nanmean(nanmean(all_TFRhann(all_group==4,correspCh,faxis>16 & faxis<29,TFtimes>.200 & TFtimes<1.0),1),3),4)); %Beta
simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
title('Younger 90% Coherence','FontSize',10)
maxAbsValues=[maxAbsValues max(max(abs(temp_topo)))];
 
subplot(1,2,2);
% temp_topo=squeeze(nanmean(nanmean(nanmean(all_TFRhann(all_group==5,correspCh,faxis>2 & faxis<4,TFtimes>.200 & TFtimes<.500),1),3),4)); %Delta
% temp_topo=squeeze(nanmean(nanmean(nanmean(all_TFRhann(all_group==5,correspCh,faxis>4 & faxis<8,TFtimes>.200 & TFtimes<.600),1),3),4)); %Theta
% temp_topo=squeeze(nanmean(nanmean(nanmean(all_TFRhann(all_group==5,correspCh,faxis>8 & faxis<11,TFtimes>.500 & TFtimes<1.20),1),3),4)); %Alpha
% temp_topo=squeeze(nanmean(nanmean(nanmean(all_TFRhann(all_group==5,correspCh,faxis>12 & faxis<16,TFtimes>.200 & TFtimes<1.0),1),3),4)); %Mu
temp_topo=squeeze(nanmean(nanmean(nanmean(all_TFRhann(all_group==5,correspCh,faxis>16 & faxis<29,TFtimes>.200 & TFtimes<1.0),1),3),4)); %Beta
simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
title('Older 90% Coherence','FontSize',10)
maxAbsValues=[maxAbsValues max(max(abs(temp_topo)))];
hb=colorbar('Position',[0.9195    0.3373    0.0143    0.2881]);

figure(f4);
caxis([-1 1]*max(maxAbsValues));