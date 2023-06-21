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
nFc=0; nFy = 0; nFo = 0;
data_young=[];
data_old=[];
count_y=0;
count_o=0;
for nF=1:length(files)
    file_name = files(nF).name;
    folder_name = files(nF).folder;
    SubID=file_name(1:end-4);
    seps=findstr(SubID,'ft_');
    SubID=SubID(seps(1)+3:end);
    tic; nFc=nFc+1;
    fprintf('... working on %s (%g/%g)\n',SubID,nF,length(files))
    
    % Set groups
    if length(SubID)==5 && SubID(3)=='8' % YOUNG - UP & DOWN 90%COH
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
    
    load([preproc_path SubID '_TF_perTrial_ICAcleaned_RespLocked.mat']); %Fixed window
    if thisagegroup==0
        nFy=nFy+1;
        young_TFRhann(nFy,:,:,:)=squeeze(nanmean(TFRhann_rlock.powspctrm_bsl,1));
    elseif thisagegroup==1
        nFo=nFo+1;
        old_TFRhann(nFo,:,:,:)=squeeze(nanmean(TFRhann_rlock.powspctrm_bsl,1));
    else
        continue;
    end
    behav_table=readtable([folder_name filesep 'behav_' file_name(findstr(file_name,'ft_')+3:end-4) '.csv']);
    meanRT(nFc)=nanmean(behav_table.RT);
    meanAcc(nFc)=(nnz(behav_table.RT>0))/(length(behav_table.RT));
    TFtimes=TFRhann_rlock.time;
    faxis=TFRhann_rlock.freq;
end

%%
cfg = [];
cfg.layout = 'acticap-64ch-standard2.mat';
cfg.channel=TFRhann_rlock.label;
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);

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

% So I need to avg across the time windows of -200-100ms
% And over the beta / theta freqs

%% Cluster Permutation
clusteralpha=0.05;
montecarloalpha=0.05;
nperm=1000;
tailFlag=0;
layout=[];

dat=cat(1,tf_groupA,tf_groupB); %tf_groupA: subj * freq * time
design=[ones(size(tf_groupA,1),1) ; 2*ones(size(tf_groupB,1),1)];


cfg_neighb=[];
cfg_neighb.method = 'template';
cfg_neighb.layout=layout;
cfg_neighb.channel=layout.label;
neighbours = ft_prepare_neighbours(cfg_neighb);

[SW_clus]=get_clusterperm_lme_lsneurom(SW_est,clusteralpha,montecarloalpha,nperm,neighbours);

%% Uncorrected t/p maps

figure; zvalim=6;
temp_topo=[]; temp_pV=[];
for nCh=1:length(newlabels)
    [~,P,~,stats]=ttest2(squeeze(nanmean(nanmean(old_TFRhann(:,nCh,faxis>4 & faxis<7,TFtimes>=-.20 & TFtimes<=0.10),3),4)), ...
        squeeze(nanmean(nanmean(young_TFRhann(:,nCh,faxis>4 & faxis<7,TFtimes>=-.20 & TFtimes<=0.10),3),4)));
    temp_topo(nCh)=stats.tstat;
    temp_pV(nCh)=P;
end
temp_topo=temp_topo(correspCh);
temp_pV=temp_pV(correspCh);
temp_topo(temp_pV>0.05)=0;
temp_topo(match_str(layout.label,{'TP7','TP8'}))=NaN;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
title(['t test - Theta Power (TF) * Group'])
colormap(parula);
colorbar; 
caxis([-1 1]*zvalim);

figure; zvalim=6;
temp_topo=[]; temp_pV=[];
for nCh=1:length(newlabels)
    [~,P,~,stats]=ttest2(squeeze(nanmean(nanmean(old_TFRhann(:,nCh,faxis>12 & faxis<29,TFtimes>=-.20 & TFtimes<=0.10),3),4)), ...
        squeeze(nanmean(nanmean(young_TFRhann(:,nCh,faxis>12 & faxis<29,TFtimes>=-.20 & TFtimes<=0.10),3),4)));
    temp_topo(nCh)=stats.tstat;
    temp_pV(nCh)=P;
end
temp_topo=temp_topo(correspCh);
temp_pV=temp_pV(correspCh);
temp_topo(temp_pV>0.05)=0;
temp_topo(match_str(layout.label,{'TP7','TP8'}))=NaN;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
title(['t test - Beta Power (TF) * Group'])
colormap(parula);
colorbar; 
caxis([-1 1]*zvalim);