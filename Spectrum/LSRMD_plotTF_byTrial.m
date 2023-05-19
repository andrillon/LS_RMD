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

%DP - in case I want to split by cohort
% CohYoungM   =   {'036M_JK','037M_JD','059M_HP','061M_LG','068M_CB','091M_SW','093M_BR','103M_JK','106M_NH','108M_CY',...
%                  '114M_CS','118M_CS','133M_DC','134M_JM','136M_JC','147M_EB','186M_AF','189M_WM','191M_DM','203M_VF',...
%                  '205M_LE','220M_NB','221M_SJ','226M_SM','234M_SW','235M_JM','240M_FM','243M_JB','265M_EZ','279M_FT',...
%                  '289M_AS','291M_KS','301M_MO','302M_BS','303M_SA','314M_LK','323M_CZ','325M_KR','328M_EW','331M_CL',...
%                  '338M_SC','339M_YV','352M_MK','353M_AT','374M_PP','377M_BL','378M_MG','384M_PD','392M_PH','398M_LO',...
%                  '400M_ED','404M_RO','414M_LA','418M_AM','422M_MK','427M_SS','439M_TM','453M_LB','458M_AH','484M_AI'};
% CohYoungT   =   {'AA_15_04_14','AC_13_05_14','AR_08_04_14','AS_23_04_14','EL_24_05_14','GW_09_05_14','JC_23_05_14',...
%                  'LK_07_04_14','MH_14_04_14','ND_16_05_14','NT_16_04_14','OF_28_04_14','OM_07_05_14','OS_09_05_14',...
%                  'PR_20_04_14','RM_06_05_14','RO_25_04_14','SB_08_05_14','SF_20_05_14','SH_25_05_14','TL_23_05_14'};
% CohYoungHN  =   {'HN867','HN868','HN870','HN871','HN872','HN873','HN874','HN875','HN876','HN877','HN878','HN879',...
%                  'HN880','HN881','HN882','HN883','HN884','HN885','HN886','HN888','HN889','HN890','HN891','HN892',...
%                  'HN893','HN894','HN895','HN896','HN897','HN898','HN899'};
% CohOldA     =   {'A101','A102','A103','A104','A105','A106','A107','A108','A110','A111','A112','A113','A114','A115',...
%                  'A116','A117','A118','A119','A120'};
% CohOldHN    =   {'HN968','HN969','HN970','HN971','HN972','HN973','HN974','HN976','HN977','HN978','HN980','HN981',...
%                  'HN982','HN983','HN985','HN986','HN987','HN988','HN989','HN990','HN992','HN993','HN994','HN995',...
%                  'HN996','HN998','HN999'};

%%
res_mat=[];
redo=0; complete=0;
datatype=1; %0 = evoked, stim locked; 1 = evoked, resp locked; 2 = induced, resp locked; 3 = induced, stim locked

% m = 1; t = 1; h = 1; a = 1; hn = 1;
%
% YoungM_pow=[]; YoungM_SNR=[]; YoungM_SNRtag=[]; YoungT_pow=[]; YoungT_SNR=[]; YoungT_SNRtag=[]; YoungHN_pow=[]; YoungHN_SNR=[]; YoungHN_SNRtag=[];
% OldA_pow=[]; OldA_SNR=[]; OldA_SNRtag=[]; OldHN_pow=[]; OldHN_SNR=[]; OldHN_SNRtag=[];
absThr=250;
nFc=0;
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
    tic;
    fprintf('... working on %s (%g/%g)\n',SubID,nF,length(files))
    
    if datatype==0 %evoked, stim locked
        load([preproc_path SubID '_TF_perTrial_ICAcleaned.mat']); %Fixed window
        nFc=nFc+1;
        TFRhann.powspctrm_log=log10(TFRhann.powspctrm);
        TFRhann.powspctrm_bsl=TFRhann.powspctrm_log-repmat(nanmean(TFRhann.powspctrm_log(:,:,:,TFRhann.time<0),4),[1 1 1 length(TFRhann.time)]);
        temp_powspctrm_bsl=squeeze(nanmean((TFRhann.powspctrm_bsl(:,:,:,TFRhann.time>-0.5 & TFRhann.time<1.5)),1));
        if size(temp_powspctrm_bsl,3)>39
            temp_powspctrm_bsl=temp_powspctrm_bsl(:,:,2:end-1);
            warning('More than expected time points - removing first and last');
        end
        all_TFRhann(nFc,:,:,:)=temp_powspctrm_bsl;
        %     all_TFRhann(nFc,:,:,:)=squeeze(nanmean((TFRhann.powspctrm_bsl(:,:,:,TFRhann.time>-0.5 & TFRhann.time<1.5)),1));    
        behav_table=readtable([folder_name filesep 'behav_' file_name(findstr(file_name,'ft_')+3:end-4) '.csv']);
        meanRT(nFc)=nanmean(behav_table.RT);
        meanAcc(nFc)=(nnz(behav_table.RT>0))/(length(behav_table.RT));
        TFtimes=TFRhann.time(TFRhann.time>-0.5 & TFRhann.time<1.5);
        faxis=TFRhann.freq;
    elseif datatype==1 % evoked, resp locked
        load([preproc_path SubID '_TF_perTrial_ICAcleaned_RespLocked.mat']); %Fixed window
        nFc=nFc+1;
        all_TFRhann(nFc,:,:,:)=squeeze(nanmean(TFRhann_rlock.powspctrm_bsl,1));
        behav_table=readtable([folder_name filesep 'behav_' file_name(findstr(file_name,'ft_')+3:end-4) '.csv']);
        meanRT(nFc)=nanmean(behav_table.RT);
        meanAcc(nFc)=(nnz(behav_table.RT>0))/(length(behav_table.RT));
        TFtimes=TFRhann_rlock.time;
        faxis=TFRhann_rlock.freq;
    elseif datatype==2 % induced, resp locked
        load([preproc_path SubID '_TF_perTrial_ICAcleaned_RespLocked_ERPremoved.mat']); %Fixed window
        nFc=nFc+1;
        all_TFRhann(nFc,:,:,:)=squeeze(nanmean(TFRhann_rlock_induced.powspctrm_bsl,1));
        behav_table=readtable([folder_name filesep 'behav_' file_name(findstr(file_name,'ft_')+3:end-4) '.csv']);
        meanRT(nFc)=nanmean(behav_table.RT);
        meanAcc(nFc)=(nnz(behav_table.RT>0))/(length(behav_table.RT));
        TFtimes=TFRhann_rlock_induced.time;
        faxis=TFRhann_rlock_induced.freq;
    elseif datatype==3 % induced, resp locked
        load([preproc_path SubID '_TF_perTrial_ICAcleaned_ERPremoved.mat']); %Fixed window
        nFc=nFc+1;
        TFRhann_induced.powspctrm_bsl=TFRhann_induced.powspctrm_log-repmat(nanmean(TFRhann_induced.powspctrm_log(:,:,:,TFRhann_induced.time<0),4),[1 1 1 length(TFRhann_induced.time)]);
        temp_powspctrm_bsl=squeeze(nanmean((TFRhann_induced.powspctrm_bsl(:,:,:,TFRhann_induced.time>-0.5 & TFRhann_induced.time<1.5)),1));
        if size(temp_powspctrm_bsl,3)>39
            temp_powspctrm_bsl=temp_powspctrm_bsl(:,:,2:end-1);
            warning('More than expected time points - removing first and last');
        end
        all_TFRhann(nFc,:,:,:)=squeeze(nanmean(TFRhann_induced.powspctrm_bsl,1));
        behav_table=readtable([folder_name filesep 'behav_' file_name(findstr(file_name,'ft_')+3:end-4) '.csv']);
        meanRT(nFc)=nanmean(behav_table.RT);
        meanAcc(nFc)=(nnz(behav_table.RT>0))/(length(behav_table.RT));
        TFtimes=TFRhann_induced.time(TFRhann_induced.time>-0.5 & TFRhann_induced.time<1.5);;
        faxis=TFRhann_induced.freq;
    end
    
    if length(SubID)==4 && SubID(1)=='A' % OLD (MONASH) - UP & DOWN 90%COH
        thisgroup=2;
        thisagegroup=1;
    elseif length(SubID)==7 % YOUNG (MONASH) - DOWN 50%COH
        thisgroup=1;
        thisagegroup=2;
    elseif length(SubID)==11 % YOUNG (TRINITY) - DOWN 50%COH
        thisgroup=3;
        thisagegroup=2;
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
%     if thisagegroup>1
%         continue;
%     end
    all_group(nFc)=thisgroup; % STORE HERE THE INFO ABOUT THE GROUP
    all_agegroup(nFc)=thisagegroup; % STORE HERE THE INFO ABOUT THE AGE GROUP
%     if thisagegroup==0
%         count_y=count_y+1;
%         TFRhann.powspctrm=TFRhann.powspctrm_bsl;
%         data_young{count_y}=TFRhann;
%     elseif thisagegroup==1
%         count_o=count_o+1;
%         TFRhann.powspctrm=TFRhann.powspctrm_bsl;
%         data_old{count_o}=TFRhann;
%     end
    %     %     save([powerspec_path filesep 'cohorts_FFT_perBlock_byElec_ICAcleaned.mat'],'YoungM_pow','YoungM_SNR','YoungM_SNRtag','YoungT_pow','YoungT_SNR','YoungT_SNRtag','YoungHN_pow','YoungHN_SNR','YoungHN_SNRtag',...
    %     %          'OldA_pow','OldA_SNR','OldA_SNRtag','OldHN_pow','OldHN_SNR','OldHN_SNRtag');
end

%%
cfg = [];
cfg.layout = 'acticap-64ch-standard2.mat';
if datatype==0
    cfg.channel=TFRhann.label;
elseif datatype==1
    cfg.channel=TFRhann_rlock.label;
elseif datatype==2
    cfg.channel=TFRhann_rlock_induced.label;
elseif datatype==3
    cfg.channel=TFRhann_induced.label;
end
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);

% cfg=[];
%             cfg.method        = 'triangulation';
%             cfg.layout        = layout;
%             cfg.feedback      = 'no';
%             cfg.channel = layout.label;
%             [neighbours] = ft_prepare_neighbours(cfg);
%
%             data_all=data_young{1};
% data_all.powspctrm=nan(count_y+count_o,size(data_all.powspctrm,2),size(data_all.powspctrm,3),size(data_all.powspctrm,4));
% for k=1:count_y
%     data_all.powspctrm(k,:,:,:)=squeeze(nanmean(data_young{k}.powspctrm,1));
% end
% for k=1:count_o
%     data_all.powspctrm(count_y+k,:,:,:)=squeeze(nanmean(data_old{k}.powspctrm,1));
% end

%%
%             cfg         = [];
%             cfg.channel = 'Pz';
%             cfg.latency = [0 1.5];
%             cfg.frequency   = [11 20]; %, can be 'all'       (default = 'all')
%             cfg.avgoverchan = 'no';%'yes' or 'no'                   (default = 'no')
%             cfg.avgovertime =  'no';%'yes' or 'no'                   (default = 'no')
%             cfg.avgoverfreq =  'no';%'yes' or 'no'                   (default = 'no')
%             cfg.parameter   = 'powspctrm';
%
% cfg.method           = 'montecarlo';
% cfg.statistic        = 'indepsamplesT';
% cfg.correctm         = 'cluster';
% cfg.clusteralpha     = 0.05;
% cfg.clusterstatistic = 'maxsum';
% cfg.minnbchan        = 1;
% cfg.neighbours       = neighbours;  % same as defined for the between-trials experiment
% cfg.tail             = 0;
% cfg.clustertail      = 0;
% cfg.alpha            = 0.05;
% cfg.numrandomization = 500;
%
% design = zeros(1, count_y + count_o);
% % design(1,:) = [1:count_y 1:count_o];
% design(1,:) = [ones(1,count_y) ones(1,count_o)*2];
%
% cfg.design = design;
% % cfg.uvar   = 1;
% cfg.ivar   = 1;
%
%
% [stat] = ft_freqstatistics(cfg, data_all);

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
if datatype==0 || datatype==3
    all_TFRhann=all_TFRhann(:,:,:,TFtimes>=0 & TFtimes<=1.5);
    TFwindow=TFtimes(TFtimes>=0);
else
    all_TFRhann=all_TFRhann(:,:,:,TFtimes>=-.70 & TFtimes<=0.50);
    TFwindow=TFtimes(TFtimes>=-.70 & TFtimes<=0.50);
end
% faxis=TFRhann.freq;

%% Plotting - Time
% channels_to_plot={'Fz','Cz','Pz','Oz'};
channels_to_plot={'Oz'};
% channels_to_plot={'O1','O2'};
%Delta
f1=figure;
set(gcf,'Position',[ 2104         115         788         574]);
format_fig;
cmap=cbrewer('seq','Blues',5);
cmap2=cbrewer('seq','Oranges',5);
for nCh=1:length(channels_to_plot)
    hold on;
    simpleTplot(TFwindow,squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==0),match_str(newlabels,channels_to_plot{nCh}),faxis>2 & faxis<4,:),2),3)),0,cmap(nCh+1,:),[0],'-',0.5,1,0,0,1);
    simpleTplot(TFwindow,squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==1),match_str(newlabels,channels_to_plot{nCh}),faxis>2 & faxis<4,:),2),3)),0,cmap2(nCh+1,:),[0],'-',0.5,1,0,0,1);
%     simpleTplot(TFtimes(TFtimes>0),squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==0),match_str(newlabels,channels_to_plot{nCh}),faxis>2 & faxis<4,TFtimes>0),2),3)),0,cmap(nCh+1,:),[0],'-',0.5,1,0,0,1);
%     simpleTplot(TFtimes(TFtimes>0),squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==1),match_str(newlabels,channels_to_plot{nCh}),faxis>2 & faxis<4,TFtimes>0),2),3)),0,cmap2(nCh+1,:),[0],'-',0.5,1,0,0,1);
    fprintf('... ... running stats  on %s\n',channels_to_plot{nCh})
    %%%% Stats diff delta
    data=[squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==0),match_str(newlabels,channels_to_plot{nCh}),faxis>2 & faxis<4,:),2),3)) ;
        squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==1),match_str(newlabels,channels_to_plot{nCh}),faxis>2 & faxis<4,:),2),3))];
    group=[all_agegroup(all_agegroup==0)' ; all_agegroup(all_agegroup==1)'];
    [realpos realneg]=get_cluster_permutation_aov(data,group,0.05,0.05,1000,TFwindow,'full');
    if size(realpos.nclusters)==1
        j=1;
        xTime=TFwindow;
        plot(xTime(realpos.clusters==realpos.nclusters(j)),zeros(1,sum(realpos.clusters==realpos.nclusters(j))),'Color','k');
    else
        for j=1:realpos.nclusters
        xTime=TFwindow;
        plot(xTime(realpos.clusters==realpos.nclusters(j)),zeros(1,sum(realpos.clusters==realpos.nclusters(j))),'Color','k');
        end
    end
end
% xlim([0 1])
% ylim([-.7 2])
xlabel('Time (s)')
ylabel('Power')
title('Delta Power - Young 90% vs Old 90%','FontSize',10)
legend(['Young ' channels_to_plot{nCh}],['Old ' channels_to_plot{nCh}],'Location','eastoutside');

%Theta
f2=figure;
set(gcf,'Position',[ 2104         115         788         574]);
format_fig;
cmap=cbrewer('seq','Blues',5);
cmap2=cbrewer('seq','Oranges',5);
for nCh=1:length(channels_to_plot)
    hold on;
    simpleTplot(TFwindow,squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==0),match_str(newlabels,channels_to_plot{nCh}),faxis>4 & faxis<8,:),2),3)),0,cmap(nCh+1,:),[0],'-',0.1,1,0,0,1);
    simpleTplot(TFwindow,squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==1),match_str(newlabels,channels_to_plot{nCh}),faxis>4 & faxis<8,:),2),3)),0,cmap2(nCh+1,:),[0],'-',0.1,1,0,0,1);
    
    data=[squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==0),match_str(newlabels,channels_to_plot{nCh}),faxis>4 & faxis<8,:),2),3)) ;
        squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==1),match_str(newlabels,channels_to_plot{nCh}),faxis>4 & faxis<8,:),2),3))];
    group=[all_agegroup(all_agegroup==0)' ; all_agegroup(all_agegroup==1)'];
    [realpos realneg]=get_cluster_permutation_aov(data,group,0.05,0.05,1000,TFwindow,'full');
    for j=1:realpos.nclusters
        xTime=TFwindow;
%         plot(xTime(realpos.clusters==realpos.nclusters(j)),zeros(1,sum(realpos.clusters==realpos.nclusters(j))),'Color','k');
        plot(xTime(realpos.clusters==(j)),zeros(1,sum(realpos.clusters==(j))),'Color','k');

    end
    
end
% xlim([0 1])
% ylim([-.7 2])
xlabel('Time (s)')
ylabel('Power')
title('Theta Power - Young 90% vs Old 90%','FontSize',10)
legend(['Young ' channels_to_plot{nCh}],['Old ' channels_to_plot{nCh}],'Location','eastoutside');

%Alpha
f3=figure;
set(gcf,'Position',[ 2104         115         788         574]);
format_fig;
cmap=cbrewer('seq','Blues',5);
cmap2=cbrewer('seq','Oranges',5);
for nCh=1:length(channels_to_plot)
    hold on;
    simpleTplot(TFwindow,squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==0),match_str(newlabels,channels_to_plot{nCh}),faxis>8 & faxis<11,:),2),3)),0,cmap(nCh+1,:),[0],'-',0.1,1,0,0,1);
    simpleTplot(TFwindow,squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==1),match_str(newlabels,channels_to_plot{nCh}),faxis>8 & faxis<11,:),2),3)),0,cmap2(nCh+1,:),[0],'-',0.1,1,0,0,1);
    
    data=[squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==0),match_str(newlabels,channels_to_plot{nCh}),faxis>8 & faxis<11,:),2),3)) ;
        squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==1),match_str(newlabels,channels_to_plot{nCh}),faxis>8 & faxis<11,:),2),3))];
    group=[all_agegroup(all_agegroup==0)' ; all_agegroup(all_agegroup==1)'];
    [realpos realneg]=get_cluster_permutation_aov(data,group,0.05,0.05,1000,TFwindow,'full');
    for j=1:realpos.nclusters
        xTime=TFwindow;
        plot(xTime(realpos.clusters==(j)),zeros(1,sum(realpos.clusters==(j))),'Color','k');
    end
end
% xlim([0 1])
% ylim([-.7 2])
xlabel('Time (s)')
ylabel('Power')
title('Alpha Power - Young 90% vs Old 90%','FontSize',10)
legend(['Young ' channels_to_plot{nCh}],['Old ' channels_to_plot{nCh}],'Location','eastoutside');

%Mu
f4=figure;
hp=[];
set(gcf,'Position',[ 2104         115         788         574]);
format_fig;
cmap=cbrewer('seq','Blues',5);
cmap2=cbrewer('seq','Oranges',5);
% for nCh=1:length(channels_to_plot)
    hold on;
    [~,hp(1)]=simpleTplot(TFwindow,squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==0),match_str(newlabels,channels_to_plot{1}),faxis>12 & faxis<16,:),2),3)),0,cmap(2,:),[0],'-',0.1,1,0,0,1);
%     [~,hp(2)]=simpleTplot(TFwindow,squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==0),match_str(newlabels,channels_to_plot{2}),faxis>12 & faxis<16,:),2),3)),0,cmap(4,:),[0],'-',0.1,1,0,1,1);
    [~,hp(3)]=simpleTplot(TFwindow,squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==1),match_str(newlabels,channels_to_plot{1}),faxis>12 & faxis<16,:),2),3)),0,cmap2(2,:),[0],'-',0.1,1,0,0,1);
%     [~,hp(4)]=simpleTplot(TFwindow,squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==1),match_str(newlabels,channels_to_plot{2}),faxis>12 & faxis<16,:),2),3)),0,cmap2(4,:),[0],'-',0.1,1,0,1,1);
% Above is for plotting across multiple electrodes

%     data=[squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==0),match_str(newlabels,channels_to_plot{nCh}),faxis>12 & faxis<16,:),2),3)) ;
%         squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==1),match_str(newlabels,channels_to_plot{nCh}),faxis>12 & faxis<16,:),2),3))];
%     group=[all_agegroup(all_agegroup==0)' ; all_agegroup(all_agegroup==1)'];
%     [realpos realneg]=get_cluster_permutation_aov(data,group,0.05,0.05,1000,TFwindow,'full');
%     for j=1:realpos.nclusters
%         xTime=TFwindow;
%         plot(xTime(realpos.clusters==(j)),zeros(1,sum(realpos.clusters==(j))),'Color','k');
%     end
% end
% xlim([0 1])
ylim([-.3 .15])
xlabel('Time (s)')
ylabel('Power')
title('Mu Power, Young 90% vs Old 90%','FontSize',10)
% legend(['Young ' channels_to_plot{1}],['Old ' channels_to_plot{1}],['Young ' channels_to_plot{2}],['Old ' channels_to_plot{2}],'Location','eastoutside');
% legend(hp,{['Young ' channels_to_plot{1}],['Young ' channels_to_plot{2}],['Old ' channels_to_plot{1}],['Old ' channels_to_plot{2}]},'Location','eastoutside');
legend('Younger','Older');

%Beta
f5=figure;
set(gcf,'Position',[ 2104         115         788         574]);
format_fig;
cmap=cbrewer('seq','Blues',5);
cmap2=cbrewer('seq','Oranges',5);
for nCh=1:length(channels_to_plot)
    hold on;
    simpleTplot(TFwindow,squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==0),match_str(newlabels,channels_to_plot{nCh}),faxis>16 & faxis<29,:),2),3)),0,cmap(nCh+1,:),[0],'-',0.1,1,0,0,1);
    simpleTplot(TFwindow,squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==1),match_str(newlabels,channels_to_plot{nCh}),faxis>16 & faxis<29,:),2),3)),0,cmap2(nCh+1,:),[0],'-',0.1,1,0,0,1);
    data=[squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==0),match_str(newlabels,channels_to_plot{nCh}),faxis>16 & faxis<29,:),2),3)) ;
        squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==1),match_str(newlabels,channels_to_plot{nCh}),faxis>16 & faxis<29,:),2),3))];
    group=[all_agegroup(all_agegroup==0)' ; all_agegroup(all_agegroup==1)'];
    [realpos realneg]=get_cluster_permutation_aov(data,group,0.05,0.05,1000,TFwindow,'full');
    for j=1:realpos.nclusters
        xTime=TFwindow;
        plot(xTime(realpos.clusters==(j)),zeros(1,sum(realpos.clusters==(j))),'Color','k');
    end
end
% xlim([0 1])
% ylim([-.7 2])
xlabel('Time (s)')
ylabel('Power')
title('Beta Power - Young 90% vs Old 90%','FontSize',10)
legend(['Young ' channels_to_plot{nCh}],['Old ' channels_to_plot{nCh}],'Location','eastoutside');

%% Plotting - Frequency
% faxis=TFRhann.freq;

% channels_to_plot={'Fz','Cz','Pz','Oz'};
channels_to_plot={'Pz'};

%0-100ms
f1=figure;
set(gcf,'Position',[ 2104         115         788         574]);
format_fig;
cmap=cbrewer('seq','Blues',5);
cmap2=cbrewer('seq','Oranges',5);

for nCh=1:length(channels_to_plot)
    hold on;
    simpleTplot(faxis,squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==0),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.30 & TFtimes <.70),2),4)),0,cmap(nCh+1,:),[0],'-',0.1,1,0,0,1);
    simpleTplot(faxis,squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==1),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.30 & TFtimes <.70),2),4)),0,cmap2(nCh+1,:),[0],'-',0.1,1,0,0,1);
    
    data=[squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==0),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.30 & TFtimes <.70),2),4));
        squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==1),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.30 & TFtimes <.70),2),4))];
    group=[all_agegroup(all_agegroup==0)' ; all_agegroup(all_agegroup==1)'];
    [realpos realneg]=get_cluster_permutation_aov(data,group,0.05,0.05,1000,faxis,'full');
    for j=1:realpos.nclusters
        plot(faxis(realpos.clusters==(j)),zeros(1,sum(realpos.clusters==(j))),'Color','k');
    end
end
% xlim([0 1])
% ylim([-.7 2])
xlabel('Frequency (Hz)')
ylabel('Power')
title('300-700ms Spectral Power - Young 90% vs Old 90%','FontSize',10)
legend(['Young '  channels_to_plot{1}],['Old '  channels_to_plot{1}],'Location','eastoutside');

%100-200ms
f2=figure;
set(gcf,'Position',[ 2104         115         788         574]);
format_fig;
cmap=cbrewer('seq','Blues',5);
cmap2=cbrewer('seq','Oranges',5);

for nCh=1:length(channels_to_plot)
    hold on;
    simpleTplot(faxis,squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==0),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.10 & TFtimes <.20),2),4)),0,cmap(nCh+1,:),[0],'-',0.1,1,0,0,1);
    simpleTplot(faxis,squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==1),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.10 & TFtimes <.20),2),4)),0,cmap2(nCh+1,:),[0],'-',0.1,1,0,0,1);
    
    data=[squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==0),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.10 & TFtimes <.20),2),4));
        squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==1),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.10 & TFtimes <.20),2),4))];
    group=[all_agegroup(all_agegroup==0)' ; all_agegroup(all_agegroup==1)'];
    [realpos realneg]=get_cluster_permutation_aov(data,group,0.05,0.05,1000,faxis,'full');
    for j=1:realpos.nclusters
        plot(faxis(realpos.clusters==(j)),zeros(1,sum(realpos.clusters==(j))),'Color','k');
    end
end
% xlim([0 1])
% ylim([-.7 2])
xlabel('Frequency (Hz)')
ylabel('Power')
title('100-200ms Spectral Power - Young 90% vs Old 90%','FontSize',10)
legend(['Young '  channels_to_plot{1}],['Old '  channels_to_plot{1}],'Location','eastoutside');

%200-300ms
f3=figure;
set(gcf,'Position',[ 2104         115         788         574]);
format_fig;
cmap=cbrewer('seq','Blues',5);
cmap2=cbrewer('seq','Oranges',5);

for nCh=1:length(channels_to_plot)
    hold on;
    simpleTplot(faxis,squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==0),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.20 & TFtimes <.30),2),4)),0,cmap(nCh+1,:),[0],'-',0.1,1,0,0,1);
    simpleTplot(faxis,squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==1),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.20 & TFtimes <.30),2),4)),0,cmap2(nCh+1,:),[0],'-',0.1,1,0,0,1);
    
    data=[squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==0),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.20 & TFtimes <.30),2),4));
        squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==1),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.20 & TFtimes <.30),2),4))];
    group=[all_agegroup(all_agegroup==0)' ; all_agegroup(all_agegroup==1)'];
    [realpos realneg]=get_cluster_permutation_aov(data,group,0.05,0.05,1000,faxis,'full');
    for j=1:realpos.nclusters
        plot(faxis(realpos.clusters==(j)),zeros(1,sum(realpos.clusters==(j))),'Color','k');
    end
end
% xlim([0 1])
% ylim([-.7 2])
xlabel('Frequency (Hz)')
ylabel('Power')
title('200-300ms Spectral Power - Young 90% vs Old 90%','FontSize',10)
legend(['Young '  channels_to_plot{1}],['Old '  channels_to_plot{1}],'Location','eastoutside');

%300-400ms
f4=figure;
set(gcf,'Position',[ 2104         115         788         574]);
format_fig;
cmap=cbrewer('seq','Blues',5);
cmap2=cbrewer('seq','Oranges',5);

for nCh=1:length(channels_to_plot)
    hold on;
    simpleTplot(faxis,squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==0),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.30 & TFtimes <.40),2),4)),0,cmap(nCh+1,:),[0],'-',0.1,1,0,0,1);
    simpleTplot(faxis,squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==1),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.30 & TFtimes <.40),2),4)),0,cmap2(nCh+1,:),[0],'-',0.1,1,0,0,1);
    
    data=[squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==0),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.30 & TFtimes <.40),2),4));
        squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==1),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.30 & TFtimes <.40),2),4))];
    group=[all_agegroup(all_agegroup==0)' ; all_agegroup(all_agegroup==1)'];
    [realpos realneg]=get_cluster_permutation_aov(data,group,0.05,0.05,1000,faxis,'full');
    for j=1:realpos.nclusters
        plot(faxis(realpos.clusters==(j)),zeros(1,sum(realpos.clusters==(j))),'Color','k');
    end
end
% xlim([0 1])
% ylim([-.7 2])
xlabel('Frequency (Hz)')
ylabel('Power')
title('300-400ms Spectral Power - Young 90% vs Old 90%','FontSize',10)
legend(['Young '  channels_to_plot{1}],['Old '  channels_to_plot{1}],'Location','eastoutside');

%400-500ms
f5=figure;
set(gcf,'Position',[ 2104         115         788         574]);
format_fig;
cmap=cbrewer('seq','Blues',5);
cmap2=cbrewer('seq','Oranges',5);

for nCh=1:length(channels_to_plot)
    hold on;
    simpleTplot(faxis,squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==0),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.40 & TFtimes <.50),2),4)),0,cmap(nCh+1,:),[0],'-',0.1,1,0,0,1);
    simpleTplot(faxis,squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==1),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.40 & TFtimes <.50),2),4)),0,cmap2(nCh+1,:),[0],'-',0.1,1,0,0,1);
    
    data=[squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==0),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.40 & TFtimes <.50),2),4));
        squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==1),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.40 & TFtimes <.50),2),4))];
    group=[all_agegroup(all_agegroup==0)' ; all_agegroup(all_agegroup==1)'];
    [realpos realneg]=get_cluster_permutation_aov(data,group,0.05,0.05,1000,faxis,'full');
    for j=1:realpos.nclusters
        plot(faxis(realpos.clusters==(j)),zeros(1,sum(realpos.clusters==(j))),'Color','k');
    end
end
% xlim([0 1])
% ylim([-.7 2])
xlabel('Frequency (Hz)')
ylabel('Power')
title('400-500ms Spectral Power - Young 90% vs Old 90%','FontSize',10)
legend(['Young '  channels_to_plot{1}],['Old '  channels_to_plot{1}],'Location','eastoutside');

%500-600ms
f6=figure;
set(gcf,'Position',[ 2104         115         788         574]);
format_fig;
cmap=cbrewer('seq','Blues',5);
cmap2=cbrewer('seq','Oranges',5);

for nCh=1:length(channels_to_plot)
    hold on;
    simpleTplot(faxis,squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==0),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.50 & TFtimes <.60),2),4)),0,cmap(nCh+1,:),[0],'-',0.1,1,0,0,1);
    simpleTplot(faxis,squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==1),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.50 & TFtimes <.60),2),4)),0,cmap2(nCh+1,:),[0],'-',0.1,1,0,0,1);
    
    data=[squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==0),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.50 & TFtimes <.60),2),4));
        squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==1),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.50 & TFtimes <.60),2),4))];
    group=[all_agegroup(all_agegroup==0)' ; all_agegroup(all_agegroup==1)'];
    [realpos realneg]=get_cluster_permutation_aov(data,group,0.05,0.05,1000,faxis,'full');
    for j=1:realpos.nclusters
        plot(faxis(realpos.clusters==(j)),zeros(1,sum(realpos.clusters==(j))),'Color','k');
    end
end
% xlim([0 1])
% ylim([-.7 2])
xlabel('Frequency (Hz)')
ylabel('Power')
title('500-600ms Spectral Power - Young 90% vs Old 90%','FontSize',10)
legend(['Young '  channels_to_plot{1}],['Old '  channels_to_plot{1}],'Location','eastoutside');

%600-700ms
f7=figure;
set(gcf,'Position',[ 2104         115         788         574]);
format_fig;
cmap=cbrewer('seq','Blues',5);
cmap2=cbrewer('seq','Oranges',5);

for nCh=1:length(channels_to_plot)
    hold on;
    simpleTplot(faxis,squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==0),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.60 & TFtimes <.70),2),4)),0,cmap(nCh+1,:),[0],'-',0.1,1,0,0,1);
    simpleTplot(faxis,squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==1),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.60 & TFtimes <.70),2),4)),0,cmap2(nCh+1,:),[0],'-',0.1,1,0,0,1);
    
    data=[squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==0),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.60 & TFtimes <.70),2),4));
        squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==1),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.60 & TFtimes <.70),2),4))];
    group=[all_agegroup(all_agegroup==0)' ; all_agegroup(all_agegroup==1)'];
    [realpos realneg]=get_cluster_permutation_aov(data,group,0.05,0.05,1000,faxis,'full');
    for j=1:realpos.nclusters
        plot(faxis(realpos.clusters==(j)),zeros(1,sum(realpos.clusters==(j))),'Color','k');
    end
end
% xlim([0 1])
% ylim([-.7 2])
xlabel('Frequency (Hz)')
ylabel('Power')
title('600-700ms Spectral Power - Young 90% vs Old 90%','FontSize',10)
legend(['Young '  channels_to_plot{1}],['Old '  channels_to_plot{1}],'Location','eastoutside');

%700-800ms
f8=figure;
set(gcf,'Position',[ 2104         115         788         574]);
format_fig;
cmap=cbrewer('seq','Blues',5);
cmap2=cbrewer('seq','Oranges',5);

for nCh=1:length(channels_to_plot)
    hold on;
    simpleTplot(faxis,squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==0),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.70 & TFtimes <.80),2),4)),0,cmap(nCh+1,:),[0],'-',0.1,1,0,0,1);
    simpleTplot(faxis,squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==1),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.70 & TFtimes <.80),2),4)),0,cmap2(nCh+1,:),[0],'-',0.1,1,0,0,1);
    
    data=[squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==0),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.70 & TFtimes <.80),2),4));
        squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==1),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.70 & TFtimes <.80),2),4))];
    group=[all_agegroup(all_agegroup==0)' ; all_agegroup(all_agegroup==1)'];
    [realpos realneg]=get_cluster_permutation_aov(data,group,0.05,0.05,1000,faxis,'full');
    for j=1:realpos.nclusters
        plot(faxis(realpos.clusters==(j)),zeros(1,sum(realpos.clusters==(j))),'Color','k');
    end
end
% xlim([0 1])
% ylim([-.7 2])
xlabel('Frequency (Hz)')
ylabel('Power')
title('700-800ms Spectral Power - Young 90% vs Old 90%','FontSize',10)
legend(['Young '  channels_to_plot{1}],['Old '  channels_to_plot{1}],'Location','eastoutside');

%800-900ms
f9=figure;
set(gcf,'Position',[ 2104         115         788         574]);
format_fig;
cmap=cbrewer('seq','Blues',5);
cmap2=cbrewer('seq','Oranges',5);

for nCh=1:length(channels_to_plot)
    hold on;
    simpleTplot(faxis,squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==0),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.80 & TFtimes <.90),2),4)),0,cmap(nCh+1,:),[0],'-',0.1,1,0,0,1);
    simpleTplot(faxis,squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==1),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.80 & TFtimes <.90),2),4)),0,cmap2(nCh+1,:),[0],'-',0.1,1,0,0,1);
    
    data=[squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==0),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.80 & TFtimes <.90),2),4));
        squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==1),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.80 & TFtimes <.90),2),4))];
    group=[all_agegroup(all_agegroup==0)' ; all_agegroup(all_agegroup==1)'];
    [realpos realneg]=get_cluster_permutation_aov(data,group,0.05,0.05,1000,faxis,'full');
    for j=1:realpos.nclusters
        plot(faxis(realpos.clusters==(j)),zeros(1,sum(realpos.clusters==(j))),'Color','k');
    end
end
% xlim([0 1])
% ylim([-.7 2])
xlabel('Frequency (Hz)')
ylabel('Power')
title('800-900ms Spectral Power - Young 90% vs Old 90%','FontSize',10)
legend(['Young '  channels_to_plot{1}],['Old '  channels_to_plot{1}],'Location','eastoutside');

%900-100ms
f10=figure;
set(gcf,'Position',[ 2104         115         788         574]);
format_fig;
cmap=cbrewer('seq','Blues',5);
cmap2=cbrewer('seq','Oranges',5);

for nCh=1:length(channels_to_plot)
    hold on;
    simpleTplot(faxis,squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==0),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.90 & TFtimes <1.0),2),4)),0,cmap(nCh+1,:),[0],'-',0.1,1,0,0,1);
    simpleTplot(faxis,squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==1),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.90 & TFtimes <1.0),2),4)),0,cmap2(nCh+1,:),[0],'-',0.1,1,0,0,1);
    
    data=[squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==0),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.90 & TFtimes <1.0),2),4));
        squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==1),match_str(newlabels,channels_to_plot{nCh}),:,TFtimes>.90 & TFtimes <1.0),2),4))];
    group=[all_agegroup(all_agegroup==0)' ; all_agegroup(all_agegroup==1)'];
    [realpos realneg]=get_cluster_permutation_aov(data,group,0.05,0.05,1000,faxis,'full');
    for j=1:realpos.nclusters
        plot(faxis(realpos.clusters==(j)),zeros(1,sum(realpos.clusters==(j))),'Color','k');
    end
end
% xlim([0 1])
% ylim([-.7 2])
xlabel('Frequency (Hz)')
ylabel('Power')
title('900-1000ms Spectral Power - Young 90% vs Old 90%','FontSize',10)
legend(['Young '  channels_to_plot{1}],['Old '  channels_to_plot{1}],'Location','eastoutside');

%% Behaviour Correlation
% So what I need to do is create a var with avg power at time/freq of
% interest for each participant, then correlate (1) overall; and (2) within
% each group separately
% Of interest: 

%Trim participants
all_TFRhann=all_TFRhann(all_agegroup<=1,:,:,:);
meanRT=meanRT(:,all_agegroup<=1);
meanAcc=meanAcc(:,all_agegroup<=1);

% Separate groups
meanRT_young90=meanRT(:,all_agegroup==0);
meanRT_old90=meanRT(:,all_agegroup==1);
meanAcc_young90=meanAcc(:,all_agegroup==0);
meanAcc_old90=meanAcc(:,all_agegroup==1);

% % Stat by chan
% channels_to_stat={'Cz'};
% 
% for nCh=1:length(channels_to_stat)
%     [pV,~,stats]=ranksum(squeeze(nanmean(nanmean(all_TFRhann(:,match_str(newlabels,channels_to_stat{nCh}),faxis > 8 & faxis < 12,TFtimes>-.30 & TFtimes <.30),3),4)), ...
%         meanAcc);
% end

% Stat and plot topo
%% Accuracy Correlation - Spearman
cmap=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap=flipud(cmap);
zvalim=4; figure;
temp_topo=[]; temp_pV=[];
for nCh=1:length(layout.label)-2
    %EVOKED
    [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(:,match_str(newlabels,layout.label(nCh)),faxis>2 & faxis<4,TFtimes>-.10 & TFtimes <.150),3),4))), ...
        meanAcc','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(:,match_str(newlabels,layout.label(nCh)),faxis>4 & faxis<8,TFtimes>-.30 & TFtimes <.10),3),4))), ...
%         meanAcc','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(:,match_str(newlabels,layout.label(nCh)),faxis>8 & faxis<11,TFtimes>-.50 & TFtimes <.0),3),4))), ...
%         meanAcc','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(:,match_str(newlabels,layout.label(nCh)),faxis>12 & faxis<16,TFtimes>-.60 & TFtimes <.10),3),4))), ...
%         meanAcc','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(:,match_str(newlabels,layout.label(nCh)),faxis>16 & faxis<29,TFtimes>-.60 & TFtimes <.0),3),4))), ...
%         meanAcc','type','spearman','rows','pairwise');

%   INDUCED
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(:,match_str(newlabels,layout.label(nCh)),faxis>2 & faxis<4,TFtimes>.0 & TFtimes <.20),3),4))), ...
%         meanAcc','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(:,match_str(newlabels,layout.label(nCh)),faxis>4 & faxis<8,TFtimes>-.30 & TFtimes <.10),3),4))), ...
%         meanAcc','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(:,match_str(newlabels,layout.label(nCh)),faxis>8 & faxis<11,TFtimes>-.60 & TFtimes <.0),3),4))), ...
%         meanAcc','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(:,match_str(newlabels,layout.label(nCh)),faxis>12 & faxis<16,TFtimes>-.70 & TFtimes <.0),3),4))), ...
%         meanAcc','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(:,match_str(newlabels,layout.label(nCh)),faxis>16 & faxis<29,TFtimes>-.60 & TFtimes <.0),3),4))), ...
%         meanAcc','type','spearman','rows','pairwise');
    temp_topo(nCh)=rho;
    temp_pV(nCh)=pV;
end
temp_topo(temp_pV>0.05)=0;
temp_topo(match_str(layout.label,{'TP7','TP8'}))=NaN;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
%EVOKED
title(['Spearmans Rho - Delta*Accuracy, Evoked Resp-Locked, -100 - 150ms'])
% title(['Spearmans Rho - Theta*Accuracy, Evoked Resp-Locked, -300 - 100ms'])
% title(['Spearmans Rho - Alpha*Accuracy, Evoked Resp-Locked, -500 - 0ms'])
% title(['Spearmans Rho - Mu*Accuracy, Evoked Resp-Locked, -600 - 100ms'])
% title(['Spearmans Rho - Beta*Accuracy, Evoked Resp-Locked, -600 - 0ms'])
%INDUCED
% title(['Spearmans Rho - Delta*Accuracy, Induced Resp-Locked, 0 - 200ms'])
% title(['Spearmans Rho - Theta*Accuracy, Induced Resp-Locked, -300 - 100ms'])
% title(['Spearmans Rho - Alpha*Accuracy, Induced Resp-Locked, -600 - 0ms'])
% title(['Spearmans Rho - Mu*Accuracy, Induced Resp-Locked, -700 - 0ms'])
% title(['Spearmans Rho - Beta*Accuracy, Induced Resp-Locked, -600 - 0ms'])
colormap(cmap);
colorbar; caxis([-1 1]*zvalim)

%% RT Correlation - Spearman
cmap=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap=flipud(cmap);
zvalim=4; figure;
temp_topo=[]; temp_pV=[];
for nCh=1:length(layout.label)-2
    %EVOKED
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(:,match_str(newlabels,layout.label(nCh)),faxis>1 & faxis<3,TFtimes>-.10 & TFtimes <.150),3),4))), ...
%         meanRT','type','spearman','rows','pairwise');
    [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(:,match_str(newlabels,layout.label(nCh)),faxis>4 & faxis<7,TFtimes>0.1 & TFtimes <.6),3),4))), ...
        meanRT','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(:,match_str(newlabels,layout.label(nCh)),faxis>8 & faxis<11,TFtimes>-.50 & TFtimes <.0),3),4))), ...
%         meanRT','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(:,match_str(newlabels,layout.label(nCh)),faxis>12 & faxis<16,TFtimes>-.60 & TFtimes <.10),3),4))), ...
%         meanRT','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(:,match_str(newlabels,layout.label(nCh)),faxis>16 & faxis<29,TFtimes>-.60 & TFtimes <.0),3),4))), ...
%         meanRT','type','spearman','rows','pairwise');
    
%   INDUCED
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(:,match_str(newlabels,layout.label(nCh)),faxis>2 & faxis<4,TFtimes>.0 & TFtimes <.20),3),4))), ...
%         meanRT','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(:,match_str(newlabels,layout.label(nCh)),faxis>4 & faxis<8,TFtimes>-.30 & TFtimes <.10),3),4))), ...
%         meanRT','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(:,match_str(newlabels,layout.label(nCh)),faxis>8 & faxis<11,TFtimes>-.60 & TFtimes <.0),3),4))), ...
%         meanRT','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(:,match_str(newlabels,layout.label(nCh)),faxis>12 & faxis<16,TFtimes>-.70 & TFtimes <.0),3),4))), ...
%         meanRT','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(:,match_str(newlabels,layout.label(nCh)),faxis>16 & faxis<29,TFtimes>-.60 & TFtimes <.0),3),4))), ...
%         meanRT','type','spearman','rows','pairwise');
    temp_topo(nCh)=rho;
    temp_pV(nCh)=pV;
end
temp_topo(temp_pV>0.05)=0;
temp_topo(match_str(layout.label,{'TP7','TP8'}))=NaN;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
%EVOKED
% title(['Spearmans Rho - Delta*RT, Evoked Resp-Locked, -100 - 150ms'])
title(['Spearmans Rho - Theta*RT, Evoked Stim-Locked'])
% title(['Spearmans Rho - Alpha*RT, Evoked Resp-Locked, -500 - 0ms'])
% title(['Spearmans Rho - Mu*RT, Evoked Resp-Locked, -600 - 100ms'])
% title(['Spearmans Rho - Beta*RT, Evoked Resp-Locked, -600 - 0ms'])

%INDUCED
% title(['Spearmans Rho - Delta*RT, Induced Resp-Locked, 0 - 200ms'])
% title(['Spearmans Rho - Theta*RT, Induced Resp-Locked, -300 - 100ms'])
% title(['Spearmans Rho - Alpha*RT, Induced Resp-Locked, -600 - 0ms'])
% title(['Spearmans Rho - Mu*RT, Induced Resp-Locked, -700 - 0ms'])
% title(['Spearmans Rho - Beta*RT, Induced Resp-Locked, -600 - 0ms'])
colormap(cmap);
colorbar; caxis([-1 1])

%% Accuracy Correlation by Group - Spearman
cmap=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap=flipud(cmap);
zvalim=4; figure;
temp_topo=[]; temp_pV=[];
subplot(1,2,1);
for nCh=1:length(layout.label)-2
    %EVOKED
    [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==0,match_str(newlabels,layout.label(nCh)),faxis>2 & faxis<4,TFtimes>-.10 & TFtimes <.150),3),4))), ...
        meanAcc_young90','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==0,match_str(newlabels,layout.label(nCh)),faxis>4 & faxis<8,TFtimes>-.30 & TFtimes <.10),3),4))), ...
%         meanAcc_young90','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==0,match_str(newlabels,layout.label(nCh)),faxis>8 & faxis<11,TFtimes>-.50 & TFtimes <.0),3),4))), ...
%         meanAcc_young90','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==0,match_str(newlabels,layout.label(nCh)),faxis>12 & faxis<16,TFtimes>-.60 & TFtimes <.10),3),4))), ...
%         meanAcc_young90','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==0,match_str(newlabels,layout.label(nCh)),faxis>16 & faxis<29,TFtimes>-.60 & TFtimes <.0),3),4))), ...
%         meanAcc_young90','type','spearman','rows','pairwise');
    
    %   INDUCED
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==0,match_str(newlabels,layout.label(nCh)),faxis>2 & faxis<4,TFtimes>.0 & TFtimes <.20),3),4))), ...
%         meanAcc_young90','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==0,match_str(newlabels,layout.label(nCh)),faxis>4 & faxis<8,TFtimes>-.30 & TFtimes <.10),3),4))), ...
%         meanAcc_young90','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==0,match_str(newlabels,layout.label(nCh)),faxis>8 & faxis<11,TFtimes>-.60 & TFtimes <.0),3),4))), ...
%         meanAcc_young90','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==0,match_str(newlabels,layout.label(nCh)),faxis>12 & faxis<16,TFtimes>-.70 & TFtimes <.0),3),4))), ...
%         meanAcc_young90','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==0,match_str(newlabels,layout.label(nCh)),faxis>16 & faxis<29,TFtimes>-.60 & TFtimes <.0),3),4))), ...
%         meanAcc_young90','type','spearman','rows','pairwise');
    temp_topo(nCh)=rho;
    temp_pV(nCh)=pV;
end
temp_topo(temp_pV>0.05)=0;
temp_topo(match_str(layout.label,{'TP7','TP8'}))=NaN;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
title('Younger 90%');
colormap(cmap);
colorbar; caxis([-1 1]*zvalim)

subplot(1,2,2);
for nCh=1:length(layout.label)-2
    %EVOKED
    [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==1,match_str(newlabels,layout.label(nCh)),faxis>2 & faxis<4,TFtimes>-.10 & TFtimes <.150),3),4))), ...
        meanAcc_old90','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==1,match_str(newlabels,layout.label(nCh)),faxis>4 & faxis<8,TFtimes>-.30 & TFtimes <.10),3),4))), ...
%         meanAcc_old90','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==1,match_str(newlabels,layout.label(nCh)),faxis>8 & faxis<11,TFtimes>-.50 & TFtimes <.0),3),4))), ...
%         meanAcc_old90','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==1,match_str(newlabels,layout.label(nCh)),faxis>12 & faxis<16,TFtimes>-.60 & TFtimes <.10),3),4))), ...
%         meanAcc_old90','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==1,match_str(newlabels,layout.label(nCh)),faxis>16 & faxis<29,TFtimes>-.60 & TFtimes <.0),3),4))), ...
%         meanAcc_old90','type','spearman','rows','pairwise');
    
    %   INDUCED
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==1,match_str(newlabels,layout.label(nCh)),faxis>2 & faxis<4,TFtimes>.0 & TFtimes <.20),3),4))), ...
%         meanAcc_old90','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==1,match_str(newlabels,layout.label(nCh)),faxis>4 & faxis<8,TFtimes>-.30 & TFtimes <.10),3),4))), ...
%         meanAcc_old90','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==1,match_str(newlabels,layout.label(nCh)),faxis>8 & faxis<11,TFtimes>-.60 & TFtimes <.0),3),4))), ...
%         meanAcc_old90','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==1,match_str(newlabels,layout.label(nCh)),faxis>12 & faxis<16,TFtimes>-.70 & TFtimes <.0),3),4))), ...
%         meanAcc_old90','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==1,match_str(newlabels,layout.label(nCh)),faxis>16 & faxis<29,TFtimes>-.60 & TFtimes <.0),3),4))), ...
%         meanAcc_old90','type','spearman','rows','pairwise');
    temp_topo(nCh)=rho;
    temp_pV(nCh)=pV;
end
temp_topo(temp_pV>0.05)=0;
temp_topo(match_str(layout.label,{'TP7','TP8'}))=NaN;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
title('Older 90%');

colormap(cmap);
colorbar; caxis([-1 1]*zvalim)

%% RT Correlation by Group - Spearman
cmap=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap=flipud(cmap);
zvalim=4; figure;
temp_topo=[]; temp_pV=[];
subplot(1,2,1);
for nCh=1:length(layout.label)-2
    %EVOKED
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==0,match_str(newlabels,layout.label(nCh)),faxis>2 & faxis<4,TFtimes>-.10 & TFtimes <.150),3),4))), ...
%         meanRT_young90','type','spearman','rows','pairwise');
    [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==0,match_str(newlabels,layout.label(nCh)),faxis>4 & faxis<7,TFtimes>0 & TFtimes <.60),3),4))), ...
        meanRT_young90','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==0,match_str(newlabels,layout.label(nCh)),faxis>8 & faxis<11,TFtimes>-.50 & TFtimes <.0),3),4))), ...
%         meanRT_young90','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==0,match_str(newlabels,layout.label(nCh)),faxis>12 & faxis<16,TFtimes>-.60 & TFtimes <.10),3),4))), ...
%         meanRT_young90','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==0,match_str(newlabels,layout.label(nCh)),faxis>16 & faxis<29,TFtimes>-.60 & TFtimes <.0),3),4))), ...
%         meanRT_young90','type','spearman','rows','pairwise');
    
    %INDUCED
    %   INDUCED
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==0,match_str(newlabels,layout.label(nCh)),faxis>2 & faxis<4,TFtimes>.0 & TFtimes <.20),3),4))), ...
%         meanRT_young90','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==0,match_str(newlabels,layout.label(nCh)),faxis>4 & faxis<8,TFtimes>-.30 & TFtimes <.10),3),4))), ...
%         meanRT_young90','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==0,match_str(newlabels,layout.label(nCh)),faxis>8 & faxis<11,TFtimes>-.60 & TFtimes <.0),3),4))), ...
%         meanRT_young90','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==0,match_str(newlabels,layout.label(nCh)),faxis>12 & faxis<16,TFtimes>-.70 & TFtimes <.0),3),4))), ...
%         meanRT_young90','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==0,match_str(newlabels,layout.label(nCh)),faxis>16 & faxis<29,TFtimes>-.60 & TFtimes <.0),3),4))), ...
%         meanRT_young90','type','spearman','rows','pairwise');
    temp_topo(nCh)=rho;
    temp_pV(nCh)=pV;
end
temp_topo(temp_pV>0.05)=0;
temp_topo(match_str(layout.label,{'TP7','TP8'}))=NaN;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
title('Younger 90%');
colormap(cmap);
colorbar; caxis([-1 1]*zvalim)

subplot(1,2,2);
for nCh=1:length(layout.label)-2
    %EVOKED
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==1,match_str(newlabels,layout.label(nCh)),faxis>2 & faxis<4,TFtimes>-.10 & TFtimes <.150),3),4))), ...
%         meanRT_old90','type','spearman','rows','pairwise');
    [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==1,match_str(newlabels,layout.label(nCh)),faxis>4 & faxis<7,TFtimes>0 & TFtimes <.60),3),4))), ...
        meanRT_old90','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==1,match_str(newlabels,layout.label(nCh)),faxis>8 & faxis<11,TFtimes>-.50 & TFtimes <.0),3),4))), ...
%         meanRT_old90','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==1,match_str(newlabels,layout.label(nCh)),faxis>12 & faxis<16,TFtimes>-.60 & TFtimes <.10),3),4))), ...
%         meanRT_old90','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==1,match_str(newlabels,layout.label(nCh)),faxis>16 & faxis<29,TFtimes>-.60 & TFtimes <.0),3),4))), ...
%         meanRT_old90','type','spearman','rows','pairwise');
    
    %   INDUCED
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==1,match_str(newlabels,layout.label(nCh)),faxis>2 & faxis<4,TFtimes>.0 & TFtimes <.20),3),4))), ...
%         meanRT_old90','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==1,match_str(newlabels,layout.label(nCh)),faxis>4 & faxis<8,TFtimes>-.30 & TFtimes <.10),3),4))), ...
%         meanRT_old90','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==1,match_str(newlabels,layout.label(nCh)),faxis>8 & faxis<11,TFtimes>-.60 & TFtimes <.0),3),4))), ...
%         meanRT_old90','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==1,match_str(newlabels,layout.label(nCh)),faxis>12 & faxis<16,TFtimes>-.70 & TFtimes <.0),3),4))), ...
%         meanRT_old90','type','spearman','rows','pairwise');
%     [rho,pV]=corr((squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==1,match_str(newlabels,layout.label(nCh)),faxis>16 & faxis<29,TFtimes>-.60 & TFtimes <.0),3),4))), ...
%         meanRT_old90','type','spearman','rows','pairwise');
    temp_topo(nCh)=rho;
    temp_pV(nCh)=pV;
end
temp_topo(temp_pV>0.05)=0;
temp_topo(match_str(layout.label,{'TP7','TP8'}))=NaN;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
title('Older 90%');

colormap(cmap);
colorbar; caxis([-1 1]*zvalim)

%% Time / frequency window topoplots

% Theta -200ms to 100ms around response
cmap=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap=flipud(cmap); maxmin=[];
temp_topo=[];

for nG=1:2
    for nCh=1:length(layout.label)-2
        temp_topo(nCh,:)=squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==nG-1,match_str(newlabels,layout.label(nCh)),faxis>4 & faxis<7,TFtimes>-.20 & TFtimes <.10),3),4));
        maxmin(nG,nCh,1)=min((temp_topo(nCh)));
        maxmin(nG,nCh,2)=max((temp_topo(nCh)));
    end
    temp_topo=[];
end

figure;
subplot(1,2,1); temp_topo=[];
for nCh=1:length(layout.label)-2
    temp_topo(nCh,:)=squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==0,match_str(newlabels,layout.label(nCh)),faxis>4 & faxis<7,TFtimes>-.20 & TFtimes <.10),3),4));
end
temp_topo(match_str(layout.label,{'TP7','TP8'}),:)=NaN;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1); title('Younger 90%');
caxis([0-max(max(max(abs(maxmin)))) max(max(max(abs(maxmin))))])
colorbar;
 
subplot(1,2,2);  temp_topo=[];
for nCh=1:length(layout.label)-2
    temp_topo(nCh,:)=squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==1,match_str(newlabels,layout.label(nCh)),faxis>4 & faxis<7,TFtimes>-.20 & TFtimes <.10),3),4));
end
temp_topo(match_str(layout.label,{'TP7','TP8'}),:)=NaN;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1); title('Older 90%');
caxis([0-max(max(max(abs(maxmin)))) max(max(max(abs(maxmin))))])
colorbar;

% Beta -200ms to 100ms around response
cmap=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap=flipud(cmap); maxmin=[];
temp_topo=[];

for nG=1:2
    for nCh=1:length(layout.label)-2
        temp_topo(nCh,:)=squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==nG-1,match_str(newlabels,layout.label(nCh)),faxis>12 & faxis<29,TFtimes>-.20 & TFtimes <.10),3),4));
        maxmin(nG,nCh,1)=min((temp_topo(nCh)));
        maxmin(nG,nCh,2)=max((temp_topo(nCh)));
    end
    temp_topo=[];
end

figure;
subplot(1,2,1); temp_topo=[];
for nCh=1:length(layout.label)-2
    temp_topo(nCh,:)=squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==0,match_str(newlabels,layout.label(nCh)),faxis>12 & faxis<29,TFtimes>-.20 & TFtimes <.10),3),4));
end
temp_topo(match_str(layout.label,{'TP7','TP8'}),:)=NaN;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1); title('Younger 90%');
caxis([0-max(max(max(abs(maxmin)))) max(max(max(abs(maxmin))))])
colorbar;
 
subplot(1,2,2);  temp_topo=[];
for nCh=1:length(layout.label)-2
    temp_topo(nCh,:)=squeeze(nanmean(nanmean(all_TFRhann(all_agegroup==1,match_str(newlabels,layout.label(nCh)),faxis>12 & faxis<29,TFtimes>-.20 & TFtimes <.10),3),4));
end
temp_topo(match_str(layout.label,{'TP7','TP8'}),:)=NaN;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1); title('Older 90%');
caxis([0-max(max(max(abs(maxmin)))) max(max(max(abs(maxmin))))])
colorbar;