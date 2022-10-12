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
redo=1; complete=0;

absThr=250;
nFc=0;
for nF=1:length(files)
    file_name = files(nF).name;
    folder_name = files(nF).folder;
    SubID=file_name(1:end-4);
    seps=findstr(SubID,'ft_');
    SubID=SubID(seps(1)+3:end);
    tic;
    fprintf('... working on %s (%g/%g)\n',SubID,nF,length(files))
    
%     if redo==1 || exist([preproc_path filesep SubID '_TF_perTrial_varwin_ICAcleaned.mat'])==0 %Variable window
    if redo==1 || exist([preproc_path SubID '_TF_perTrial_ICAcleaned.mat'])==0 %Fixed window

        subj_pow=[];

        
        load([folder_name filesep file_name]);
        
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

        % crop trials
        cfg=[];
        cfg.toilim =[-0.5 3.5];
        data=ft_redefinetrial(cfg,data);
        
        % estimate the ERP
        tlck = ft_timelockanalysis([],data);
        figure; plot(tlck.time, tlck.avg); legend(tlck.label);
        
        % subtract the ERP from the data
        data_minus_erp = data;
        for k = 1:numel(data.trial)
            data_minus_erp.trial{k} = data.trial{k} - tlck.avg(:,1:length(data.trial{k}));
        end
        
        % get the TF decomposition
        cfg              = [];
        cfg.output       = 'pow';
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.foi          = 2:0.5:30;                         % analysis 2 to 30 Hz in steps of 2 Hz
        cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
%         cfg.t_ftimwin    = 2./cfg.foi;   % 2 cycles per time window
        cfg.toi          = -1.0:0.05:3.5;                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
        cfg.keeptrials   = 'yes';
        TFRhann_induced = ft_freqanalysis(cfg, data_minus_erp);

%         save([preproc_path filesep SubID '_TF_perTrial_varwin_ICAcleaned.mat'],'TFRhann','newlabels'); %Variable window
        save([preproc_path filesep SubID '_TF_perTrial_ICAcleaned_ERPremoved.mat'],'TFRhann_induced','newlabels'); %Fixed window
        
    else
%         load([preproc_path filesep SubID '_TF_perTrial_varwin_ICAcleaned.mat']); %Variable window
        load([preproc_path filesep SubID '_TF_perTrial_ICAcleaned_ERPremoved.mat']); %Fixed window
    end
%   

    nFc=nFc+1;
    TFRhann_induced.powspctrm_log=log10(TFRhann_induced.powspctrm);
    TFRhann_induced.powspctrm_bsl=TFRhann_induced.powspctrm_log-repmat(nanmean(TFRhann_induced.powspctrm_log(:,:,:,TFRhann_induced.time<0),4),[1 1 1 length(TFRhann_induced.time)]);
    temp_powspctrm_bsl=squeeze(nanmean((TFRhann_induced.powspctrm_bsl(:,:,:,TFRhann_induced.time>-0.5 & TFRhann_induced.time<1.5)),1));
    if size(temp_powspctrm_bsl,3)>39
        temp_powspctrm_bsl=temp_powspctrm_bsl(:,:,2:end-1);
        warning('More than expected time points - removing first and last');
    end
    all_TFRhann_induced(nFc,:,:,:)=temp_powspctrm_bsl;
%     all_TFRhann(nFc,:,:,:)=squeeze(nanmean((TFRhann.powspctrm_bsl(:,:,:,TFRhann.time>-0.5 & TFRhann.time<1.5)),1));
    TFtimes=TFRhann_induced.time(TFRhann_induced.time>-0.5 & TFRhann_induced.time<1.5);
% if exist('subj_osci')~=0
%         all_osci(nFc,:,:,:)=subj_osci;
%     end
%     if exist('subj_mixed')~=0
%         all_mixed(nFc,:,:,:)=subj_mixed;
%     end
%     if exist('subj_frac')~=0
%         all_frac(nFc,:,:,:)=subj_frac;
%     end
% 
%     

% Response locking

behavfiles=dir([data_path 'Behav' filesep file_name(1:end-4) '.mat']); %DP 28/07 adding Megan older adult + younger adult data
file_name = this_file(1).name;
file_folder = this_file(1).folder;
FileID=file_name(1:end-4);
fprintf('... working on subject %s (%g/%g) - file %s (%g/%g)\n',SubID,nF,length(folders),FileID,k,length(files))

if exist('behavfiles')~=0 && ~isempty(behavfiles)
    behav_data=load([behavfiles.folder filesep behavfiles.name]);
else
    behav_data=[];
    warning('cannot find the behavioural data!!');
end
cfg = [];
cfg.offset = XXX
TFRhann_induced_rlck = ft_redefinetrial(cfg,data_minus_erp)

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
%     %     save([powerspec_path filesep 'cohorts_FFT_perBlock_byElec_ICAcleaned.mat'],'YoungM_pow','YoungM_SNR','YoungM_SNRtag','YoungT_pow','YoungT_SNR','YoungT_SNRtag','YoungHN_pow','YoungHN_SNR','YoungHN_SNRtag',...
%     %          'OldA_pow','OldA_SNR','OldA_SNRtag','OldHN_pow','OldHN_SNR','OldHN_SNRtag');
    
end
% save([preproc_path filesep 'all_FFT_perBlock_byElec_ICAcleaned_v2.mat'],'all_pow','all_group','all_agegroup');

% %%
% %All participants
% myLabels={'Fz','Cz','Pz','Oz'};
% figure;
% maxAbsValues=[];
% for nCh=1:length(myLabels)
%     subplot(2,2,nCh);
%     simpleTFplot(squeeze(nanmean(all_TFRhann(:,match_str(newlabels,myLabels{nCh}),:,:),1)),TFRhann.freq,TFtimes,0,0);
%     format_fig;
%     title(myLabels{nCh})
%     colorbar;
%     maxAbsValues=[maxAbsValues max(max(abs(squeeze(nanmean(all_TFRhann(:,match_str(newlabels,myLabels{nCh}),:,:),1)))))];
% end
% for nCh=1:length(myLabels)
%     subplot(2,2,nCh);
% caxis([-1 1]*max(maxAbsValues))
% end
% %All younger participants
% figure;
% for nCh=1:length(myLabels)
%     subplot(2,2,nCh);
%     simpleTFplot(squeeze(nanmean(all_TFRhann((all_group==1|all_group==3|all_group==4),match_str(newlabels,myLabels{nCh}),:,:),1)),TFRhann.freq,TFtimes,0,0);
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
%     simpleTFplot(squeeze(nanmean(all_TFRhann((all_group==1|all_group==3),match_str(newlabels,myLabels{nCh}),:,:),1)),TFRhann.freq,TFtimes,0,0);
%     format_fig;
%     title(myLabels{nCh})
% colorbar;
%     maxAbsValues=[maxAbsValues max(max(abs(squeeze(nanmean(all_TFRhann(:,match_str(newlabels,myLabels{nCh}),:,:),1)))))];
% end
% 
% %Younger 90% coherence
% f2=figure;
% for nCh=1:length(myLabels)
%     subplot(2,2,nCh);
%     simpleTFplot(squeeze(nanmean(all_TFRhann((all_group==4),match_str(newlabels,myLabels{nCh}),:,:),1)),TFRhann.freq,TFtimes,0,0);
%     format_fig;
%     title(myLabels{nCh})
% colorbar;
%     maxAbsValues=[maxAbsValues max(max(abs(squeeze(nanmean(all_TFRhann(:,match_str(newlabels,myLabels{nCh}),:,:),1)))))];
% end
% 
% %Older 90% coherence
% f3=figure;
% for nCh=1:length(myLabels)
%     subplot(2,2,nCh);
%     simpleTFplot(squeeze(nanmean(all_TFRhann((all_group==5),match_str(newlabels,myLabels{nCh}),:,:),1)),TFRhann.freq,TFtimes,0,0);
%     format_fig;
%     title(myLabels{nCh})
% colorbar;
%     maxAbsValues=[maxAbsValues max(max(abs(squeeze(nanmean(all_TFRhann(:,match_str(newlabels,myLabels{nCh}),:,:),1)))))];
% end
% 
% for nCh=1:length(myLabels)
%     figure(f1);
%     subplot(2,2,nCh);
% caxis([-1 1]*max(maxAbsValues))
% 
% figure(f2);
%     subplot(2,2,nCh);
% caxis([-1 1]*max(maxAbsValues))
% 
% figure(f3);
%     subplot(2,2,nCh);
% caxis([-1 1]*max(maxAbsValues))
% end
% %% Difference TF
% % Young 90% - Old 90%
% f2=figure;
% for nCh=1:length(myLabels)
%     subplot(2,2,nCh);
%     TF_groupA=squeeze(nanmean(all_TFRhann((all_group==4),match_str(newlabels,myLabels{nCh}),:,:),1));
%     TF_groupB=squeeze(nanmean(all_TFRhann((all_group==5),match_str(newlabels,myLabels{nCh}),:,:),1));
%     simpleTFplot(TF_groupA-TF_groupB,TFRhann.freq,TFtimes,0,0);
%     format_fig;
%     title(myLabels{nCh})
% colorbar;
% end
% 
% %Young 90% - Young 50%
% for nCh=1:length(myLabels)
%     subplot(2,2,nCh);
%     TF_groupA=squeeze(nanmean(all_TFRhann((all_group==4),match_str(newlabels,myLabels{nCh}),:,:),1));
%     TF_groupB=squeeze(nanmean(all_TFRhann((all_group==1|all_group==3),match_str(newlabels,myLabels{nCh}),:,:),1));
%     simpleTFplot(TF_groupA-TF_groupB,TFRhann.freq,TFtimes,0,0);
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
%     simpleTFplot(TF_groupA-TF_groupB,TFRhann.freq,TFtimes,0,0);
%     format_fig;
%     title(myLabels{nCh})
% colorbar;
% end
% 
%  %% Topographies 25Hz tag
% cfg = [];
% cfg.layout = 'biosemi64.lay';
% cfg.channel=newlabels;
% cfg.center      = 'yes';
% layout=ft_prepare_layout(cfg);
% correspCh=[];
% for nCh=1:length(layout.label)-2
%     correspCh(nCh)=match_str(newlabels,layout.label(nCh));
% end
% 
% maxAbsValues=[];
% 
% f4=figure;
% cmap=cbrewer('seq','YlOrRd',64); % select a sequential colorscale from yellow to red (64)
% cmap(cmap<0)=0;
% faxis=TFRhann.freq
% faxis=faxis(faxis>=min(TFRhann.freq) & faxis<=max(TFRhann.freq));
% 
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
% 
% subplot(1,3,2);
% % temp_topo=squeeze(nanmean(nanmean(nanmean(all_TFRhann(all_group==4,correspCh,faxis>2 & faxis<4,TFtimes>.200 & TFtimes<.500),1),3),4)); %Delta
% % temp_topo=squeeze(nanmean(nanmean(nanmean(all_TFRhann(all_group==4,correspCh,faxis>4 & faxis<8,TFtimes>.200 & TFtimes<.600),1),3),4)); %Theta
% % temp_topo=squeeze(nanmean(nanmean(nanmean(all_TFRhann(all_group==4,correspCh,faxis>8 & faxis<11,TFtimes>.500 & TFtimes<1.20),1),3),4)); %Alpha
% % temp_topo=squeeze(nanmean(nanmean(nanmean(all_TFRhann(all_group==4,correspCh,faxis>12 & faxis<16,TFtimes>.200 & TFtimes<1.0),1),3),4)); %Mu
% temp_topo=squeeze(nanmean(nanmean(nanmean(all_TFRhann(all_group==4,correspCh,faxis>16 & faxis<29,TFtimes>.200 & TFtimes<1.0),1),3),4)); %Beta
% simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
% title('Younger 90% Coherence','FontSize',10)
% maxAbsValues=[maxAbsValues max(max(abs(temp_topo)))];
%  
% subplot(1,3,3);
% % temp_topo=squeeze(nanmean(nanmean(nanmean(all_TFRhann(all_group==5,correspCh,faxis>2 & faxis<4,TFtimes>.200 & TFtimes<.500),1),3),4)); %Delta
% % temp_topo=squeeze(nanmean(nanmean(nanmean(all_TFRhann(all_group==5,correspCh,faxis>4 & faxis<8,TFtimes>.200 & TFtimes<.600),1),3),4)); %Theta
% % temp_topo=squeeze(nanmean(nanmean(nanmean(all_TFRhann(all_group==5,correspCh,faxis>8 & faxis<11,TFtimes>.500 & TFtimes<1.20),1),3),4)); %Alpha
% % temp_topo=squeeze(nanmean(nanmean(nanmean(all_TFRhann(all_group==5,correspCh,faxis>12 & faxis<16,TFtimes>.200 & TFtimes<1.0),1),3),4)); %Mu
% temp_topo=squeeze(nanmean(nanmean(nanmean(all_TFRhann(all_group==5,correspCh,faxis>16 & faxis<29,TFtimes>.200 & TFtimes<1.0),1),3),4)); %Beta
% simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
% title('Older 90% Coherence','FontSize',10)
% maxAbsValues=[maxAbsValues max(max(abs(temp_topo)))];
% hb=colorbar('Position',[0.9195    0.3373    0.0143    0.2881]);
% 
% figure(f4);
% caxis([-1 1]*max(maxAbsValues));