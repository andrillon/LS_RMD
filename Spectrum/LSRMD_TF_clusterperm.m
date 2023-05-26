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

%% Cluster Permutation
clusteralpha=0.05;
montecarloalpha=0.05;
nperm=1000; nFreq=length(faxis); nTime=length(TFwindow);
tailFlag=0;
channels_to_plot={'Fz'};

dat=cat(1,squeeze(young_TFRhann(:,match_str(newlabels,channels_to_plot{1}),:,:)),squeeze(old_TFRhann(:,match_str(newlabels,channels_to_plot{1}),:,:))); %tf_groupA: subj * freq * time
design=[ones(size(young_TFRhann,1),1) ; 2*ones(size(old_TFRhann,1),1)];

cfg                     = [];
cfg.neighbours          = [];
cfg.dimord              = 'rpt_freq_times';
cfg.freq                = faxis;
cfg.time                = TFwindow;
cfg.channel             = {'eeg'};
cfg.dim                 = [1 nFreq nTime];

cfg.numrandomization    = nperm;                   % number of randomizations, can be 'all'
cfg.correctm            = 'cluster';              % apply multiple-comparison correction, 'no', 'max', cluster', 'bonferoni', 'holms', 'fdr' (default = 'no')
cfg.alpha               = montecarloalpha;                   % critical value for rejecting the null-hypothesis per tail (default = 0.05)
cfg.clusteralpha        = clusteralpha;                   % critical value for rejecting the null-hypothesis per tail (default = 0.05)
cfg.tail                = tailFlag;                      % -1, 1 or 0 (default = 0)
cfg.correcttail         = 'no';                % correct p-values or alpha-values when doing a two-sided test, 'alpha','prob' or 'no' (default = 'no')
cfg.ivar                = 1;                      % number or list with indices, independent variable(s)
cfg.uvar                = 2;                      % number or list with indices, unit variable(s)
cfg.wvar                = [];                     % number or list with indices, within-cell variable(s)
cfg.cvar                = [];                     % number or list with indices, control variable(s)
cfg.feedback            = 'textbar';              % 'gui', 'text', 'textbar' or 'no' (default = 'text')
cfg.randomseed          = 'yes';                  % 'yes', 'no' or a number (default = 'yes')

cfg.clusterstatistic    = 'maxsum';               % no, max, maxsize, maxsum, wcm init: maxsum
cfg.clusterthreshold    = 'nonparametric_common'; % parametric, nonparametric_individual, nonparametric_common
cfg.clusteralpha        = clusteralpha;                   %
cfg.clustercritval      = [] ;                    %
cfg.clustertail         = cfg.tail;               %

cfg.statistic           = 'depsamplesT';
cfg.avgoverchan         = 'no';
[stat, cfg]             = ft_statistics_montecarlo(cfg, dat, design);

%%% ploting result
TF_diff=squeeze(mean(tf_groupA,1)-mean(tf_groupB,1));
% mytval(~stat.mask)=0;
h=simpleTFplot(TF_diff,freqP,timeP,0,0); %colormap(cmap); %ylim([freqs(1) 8])
caxis([-1 1]*max(max(abs(TF_diff))));
newmask = 1*double(reshape(stat.mask, [length(timeP) length(freqP)]))';
newmask(newmask(:) == 0) = 0.5;
alpha(h, newmask); hold on;
if ~isempty(newmask==1)
    contour(timeP, freqP,newmask,0.5,'Color','k','LineWidth',3)
end