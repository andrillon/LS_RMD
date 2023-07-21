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

%% Trimming Window
young_TFRhann=young_TFRhann(:,:,:,TFtimes>=-.70 & TFtimes<=0.50);
old_TFRhann=old_TFRhann(:,:,:,TFtimes>=-.70 & TFtimes<=0.50);
TFwindow=TFtimes(TFtimes>=-.70 & TFtimes<=0.50);
% faxis=TFRhann.freq;

%% Trim to electrode of interest & mean
channels_to_plot={'Fz'};
young_TFRhann=squeeze(young_TFRhann(:,match_str(newlabels,channels_to_plot{1}),:,:));
old_TFRhann=squeeze(old_TFRhann(:,match_str(newlabels,channels_to_plot{1}),:,:));

%% Cluster Permutation
% clusteralpha=0.05;
% montecarloalpha=0.05;
% nperm=1000; nFreq=length(faxis); nTime=length(TFwindow);
% tailFlag=0;
% channels_to_plot={'Fz'};
% 
% dat=cat(1,young_TFRhann,old_TFRhann); %tf_groupA: subj * freq * time
% design=[ones(size(young_TFRhann,1),1) ; 2*ones(size(old_TFRhann,1),1)];
% 
% cfg                     = [];
% cfg.neighbours          = [];
% cfg.dimord              = 'rpt_freq_times';
% cfg.freq                = faxis;
% cfg.time                = TFwindow;
% cfg.channel             = {'eeg'};
% cfg.dim                 = [nFreq nTime];
% cfg.resampling          = 'permutation';
% 
% cfg.numrandomization    = nperm;                   % number of randomizations, can be 'all'
% cfg.correctm            = 'cluster';              % apply multiple-comparison correction, 'no', 'max', cluster', 'bonferoni', 'holms', 'fdr' (default = 'no')
% cfg.alpha               = montecarloalpha;                   % critical value for rejecting the null-hypothesis per tail (default = 0.05)
% cfg.clusteralpha        = clusteralpha;                   % critical value for rejecting the null-hypothesis per tail (default = 0.05)
% cfg.tail                = tailFlag;                      % -1, 1 or 0 (default = 0)
% cfg.correcttail         = 'no';                % correct p-values or alpha-values when doing a two-sided test, 'alpha','prob' or 'no' (default = 'no')
% cfg.ivar                = 1;                      % number or list with indices, independent variable(s)
% cfg.uvar                = [];                      % number or list with indices, unit variable(s)
% cfg.wvar                = [];                     % number or list with indices, within-cell variable(s)
% cfg.cvar                = [];                     % number or list with indices, control variable(s)
% cfg.feedback            = 'textbar';              % 'gui', 'text', 'textbar' or 'no' (default = 'text')
% cfg.randomseed          = 'yes';                  % 'yes', 'no' or a number (default = 'yes')
% 
% cfg.clusterstatistic    = 'maxsum';               % no, max, maxsize, maxsum, wcm init: maxsum
% cfg.clusterthreshold    = 'nonparametric_common'; % parametric, nonparametric_individual, nonparametric_common
% cfg.clusteralpha        = clusteralpha;                   %
% cfg.clustercritval      = [] ;                    %
% cfg.clustertail         = cfg.tail;               %
% 
% cfg.statistic           = 'indepsamplesT';
% cfg.avgoverchan         = 'no';
% [stat, cfg]             = ft_statistics_montecarlo(cfg, dat, design);
% --------------------------------------------------------------------------------------------------
%% Mana Added
% Load the structure that I sent to you
load([preproc_path filesep 'TFR_struct.mat']);

clusteralpha = 0.01;
montecarloalpha = 0.05;
nperm = 1000;
design = [ones(1,size(old_TFRhann,1)) , 2*ones(1,size(young_TFRhann,1))];

TFR_young = TFR_struct;
pow_young  = zeros([size(young_TFRhann,1), 1, size(young_TFRhann,2), size(young_TFRhann,3)]);
pow_young(:,1,:,:) = young_TFRhann;
TFR_young.powspctrm = pow_young;

TFR_old = TFR_struct;
pow_old  = zeros([size(old_TFRhann,1), 1, size(old_TFRhann,2), size(old_TFRhann,3)]);
pow_old(:,1,:,:) = old_TFRhann;
TFR_old.powspctrm = pow_old;

cfg = [];
cfg.channel          =  'all';
cfg.latency          = 'all';
cfg.avgovertime      = 'no';
cfg.avgoverfreq      ='no';
cfg.avgoverchan      = 'no';
cfg.data = 'powspctrm';
cfg.neighbours = [];
cfg.clusteralpha     = clusteralpha;
cfg.montecarloalpha = montecarloalpha;
cfg.statistic        = 'ft_statfun_indepsamplesT';
cfg.numrandomization = nperm;
cfg.correctm         =   'cluster';
cfg.method           = 'montecarlo';
cfg.design           = design;
cfg.ivar                = 1;                      
cfg.uvar                = [];                     
cfg.wvar                = [];                    
cfg.cvar                = [];                                     
stat = ft_freqstatistics(cfg, TFR_old, TFR_young);

%---------------------------------------------------------------------------------------------------


%% plotting result
stat.newmask=zeros(size(stat.mask));
stat.newmask(stat.negclusterslabelmat~=0)=1;
stat.newmask(stat.posclusterslabelmat~=0)=1;
TF_diff=squeeze(mean(TFR_old.powspctrm,1)-mean(TFR_young.powspctrm,1));
% mytval(~stat.mask)=0;
h=simpleTFplot(TF_diff,faxis,TFwindow,0,0); %colormap(cmap); %ylim([freqs(1) 8])
caxis([-1 1]*max(max(abs(TF_diff))));
newmask = squeeze(stat.newmask);
newmask(newmask(:) == 0) = .5;
alpha(h, newmask); hold on;
if ~isempty(newmask==1)
    contour(TFwindow, faxis,newmask,.5,'Color','k','LineWidth',3)
end

%%% pos clusters
pos_clus_labels=unique(stat.posclusterslabelmat);
pos_clus_labels(pos_clus_labels==0)=[];
for nclus=1:length(pos_clus_labels)
    temp_mat=repmat(stat.posclusterslabelmat,[size(old_TFRhann,1) 1 1]);
    temp_mat(temp_mat~=pos_clus_labels(nclus))=0;
    temp_mat(temp_mat==pos_clus_labels(nclus))=1;
    temp_pow=old_TFRhann;
    temp_pow(temp_mat==0)=NaN;
    TF_Old_PosClus{nclus}=squeeze(nanmean(nanmean(temp_pow,2),3));
    temp_mat=repmat(stat.posclusterslabelmat,[size(young_TFRhann,1) 1 1]);
    temp_mat(temp_mat~=pos_clus_labels(nclus))=0;
    temp_mat(temp_mat==pos_clus_labels(nclus))=1;
    temp_pow=young_TFRhann;
    temp_pow(temp_mat==0)=NaN;
    TF_Young_PosClus{nclus}=squeeze(nanmean(nanmean(temp_pow,2),3));
end

neg_clus_labels=unique(stat.negclusterslabelmat);
neg_clus_labels(neg_clus_labels==0)=[];
for nclus=1:length(neg_clus_labels)
    temp_mat=repmat(stat.negclusterslabelmat,[size(old_TFRhann,1) 1 1]);
    temp_mat(temp_mat~=neg_clus_labels(nclus))=0;
    temp_mat(temp_mat==neg_clus_labels(nclus))=1;
    temp_pow=old_TFRhann;
    temp_pow(temp_mat==0)=NaN;
    TF_Old_NegClus{nclus}=squeeze(nanmean(nanmean(temp_pow,2),3));
    temp_mat=repmat(stat.negclusterslabelmat,[size(young_TFRhann,1) 1 1]);
    temp_mat(temp_mat~=neg_clus_labels(nclus))=0;
    temp_mat(temp_mat==neg_clus_labels(nclus))=1;
    temp_pow=young_TFRhann;
    temp_pow(temp_mat==0)=NaN;
    TF_Young_NegClus{nclus}=squeeze(nanmean(nanmean(temp_pow,2),3));
end

% trying multicolour contours
 cmap=cbrewer('seq','oranges',10);
 cmap_negclus={cmap(1:3,:)}; %set colors
cmap_posclus={cmap(end-2:end,:)};

figure;
h1=simpleTFplot(TF_diff,faxis,TFwindow,0,0); %recreate TF diff plot
caxis([-1 1]*max(max(abs(TF_diff))));
% cmap_negclus={[1 0 0],[0 1 0],[0 0 1],[0 1 1]}; %set colors
% cmap_posclus={[1 0 1],[1 1 0]};

stat.newmask2=zeros(size(stat.mask));
stat.newmask2(stat.negclusterslabelmat~=0 | stat.posclusterslabelmat~=0)=1;
newmask = squeeze(stat.newmask2);
newmask(newmask(:) == 0) = .5;
alpha(h1, newmask); hold on;

for nclus=1:length(neg_clus_labels)
    stat.newmask=zeros(size(stat.mask));
    stat.newmask(stat.negclusterslabelmat==nclus)=nclus;
    newmask = squeeze(stat.newmask);
    newmask(newmask(:) == 0) = .5;
    if ~isempty(newmask==nclus)
        contour(TFwindow, 1:length(faxis),newmask,1,'Color',cmap_negclus{1}(nclus,:),'LineWidth',3)
    end
end

for nclus=1:length(pos_clus_labels)
    stat.newmask=zeros(size(stat.mask));
    stat.newmask(stat.posclusterslabelmat==nclus)=nclus;
    newmask = squeeze(stat.newmask);
    newmask(newmask(:) == 0) = .5;
    if ~isempty(newmask==nclus)
        contour(TFwindow, 1:length(faxis),newmask,1,'Color',cmap_posclus{1}(nclus,:),'LineWidth',3)
    end
end

%% Uncorrected T/p Maps
%{
[H,P,CI,STATS] = ttest2(old_TFRhann,young_TFRhann,'dim',1); %stats
t=squeeze(STATS.tstat);
t(P>.05)=0;

figure;
h1=simpleTFplot(t,faxis,TFwindow,0,0); colorbar;
%}

%% Behaviour Correlation

clearvars('TF_All_PosClus','pos_clus_stats','TF_All_NegClus','neg_clus_stats');

for nC=1:length(pos_clus_labels)
    TF_All_PosClus=[TF_Old_PosClus{nC}' TF_Young_PosClus{nC}'];
    [rho,pV]=corr(TF_All_PosClus', meanRT','type','spearman','rows','pairwise');
    pos_clus_stats(nC,1)=rho;
    pos_clus_stats(nC,2)=pV;
end

for nC=1:length(neg_clus_labels)
    TF_All_NegClus=[TF_Old_NegClus{nC}' TF_Young_NegClus{nC}'];
    [rho,pV]=corr(TF_All_NegClus', meanRT','type','spearman','rows','pairwise');
    neg_clus_stats(nC,1)=rho;
    neg_clus_stats(nC,2)=pV;
end

%% CPP Stats Correlation

clearvars('TF_All_PosClus','pos_clus_CPPons_stats','TF_All_NegClus','neg_clus_CPPons_stats','pos_clus_CPPrslope_stats','neg_clus_CPPrslope_stats');
load([preproc_path 'CPP_stats']);

for nC=1:length(pos_clus_labels)
    TF_All_PosClus=[TF_Old_PosClus{nC}' TF_Young_PosClus{nC}'];
    [rho,pV]=corr(TF_All_PosClus', CPP_onsets','type','spearman','rows','pairwise');
    pos_clus_CPPons_stats(nC,1)=rho;
    pos_clus_CPPons_stats(nC,2)=pV;
end

for nC=1:length(neg_clus_labels)
    TF_All_NegClus=[TF_Old_NegClus{nC}' TF_Young_NegClus{nC}'];
    [rho,pV]=corr(TF_All_NegClus', CPP_onsets','type','spearman','rows','pairwise');
    neg_clus_CPPons_stats(nC,1)=rho;
    neg_clus_CPPons_stats(nC,2)=pV;
end

for nC=1:length(pos_clus_labels)
    TF_All_PosClus=[TF_Old_PosClus{nC}' TF_Young_PosClus{nC}'];
    [rho,pV]=corr(TF_All_PosClus', CPPr_slopes','type','spearman','rows','pairwise');
    pos_clus_CPPrslope_stats(nC,1)=rho;
    pos_clus_CPPrslope_stats(nC,2)=pV;
end

for nC=1:length(neg_clus_labels)
    TF_All_NegClus=[TF_Old_NegClus{nC}' TF_Young_NegClus{nC}'];
    [rho,pV]=corr(TF_All_NegClus', CPPr_slopes','type','spearman','rows','pairwise');
    neg_clus_CPPrslope_stats(nC,1)=rho;
    neg_clus_CPPrslope_stats(nC,2)=pV;
end