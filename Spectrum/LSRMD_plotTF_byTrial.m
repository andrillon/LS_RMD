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

% m = 1; t = 1; h = 1; a = 1; hn = 1;
%
% YoungM_pow=[]; YoungM_SNR=[]; YoungM_SNRtag=[]; YoungT_pow=[]; YoungT_SNR=[]; YoungT_SNRtag=[]; YoungHN_pow=[]; YoungHN_SNR=[]; YoungHN_SNRtag=[];
% OldA_pow=[]; OldA_SNR=[]; OldA_SNRtag=[]; OldHN_pow=[]; OldHN_SNR=[]; OldHN_SNRtag=[];
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
    
    load([preproc_path filesep SubID '_TF_perTrial_ICAcleaned.mat']); %Fixed window
    %
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
    TFtimes=TFRhann.time(TFRhann.time>-0.5 & TFRhann.time<1.5);
    
    
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
    all_group(nFc)=thisgroup; % STORE HERE THE INFO ABOUT THE GROUP
    all_agegroup(nFc)=thisagegroup; % STORE HERE THE INFO ABOUT THE AGE GROUP
    %     %     save([powerspec_path filesep 'cohorts_FFT_perBlock_byElec_ICAcleaned.mat'],'YoungM_pow','YoungM_SNR','YoungM_SNRtag','YoungT_pow','YoungT_SNR','YoungT_SNRtag','YoungHN_pow','YoungHN_SNR','YoungHN_SNRtag',...
    %     %          'OldA_pow','OldA_SNR','OldA_SNRtag','OldHN_pow','OldHN_SNR','OldHN_SNRtag');
end

%% Topographies 25Hz tag
cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.channel=newlabels;
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);
correspCh=[];
for nCh=1:length(layout.label)-2
    correspCh(nCh)=match_str(newlabels,layout.label(nCh));
end

%% Plotting
faxis=TFRhann.freq
faxis=faxis(faxis>=min(TFRhann.freq) & faxis<=max(TFRhann.freq));

channels_to_plot={'Fz','Cz','Pz','Oz'};
figure; 
set(gcf,'Position',[ 2104         115         788         574]);
format_fig;
cmap=cbrewer('seq','YlOrRd',64);
    for nCh=1:length(channels_to_plot)
        hold on;
        simpleTplot(TFtimes(TFtimes>0 & TFtimes<.650),squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==0),match_str(newlabels,channels_to_plot{nCh}),faxis>4 & faxis<8,TFtimes>0 & TFtimes<.650),1),3)),0,cmap(nCh+1,:),[0],'-',0.1,1,0,1,1);
        simpleTplot(TFtimes(TFtimes>0 & TFtimes<.650),squeeze(nanmean(nanmean(all_TFRhann((all_agegroup==1),match_str(newlabels,channels_to_plot{nCh}),faxis>4 & faxis<8,TFtimes>0 & TFtimes<.650),1),3)),0,cmap(nCh+30,:),[0],'-',0.1,1,0,1,1);
    end
% xlim([0 1])
% ylim([-.7 2])
xlabel('Time')
ylabel('Power')