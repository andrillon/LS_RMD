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
files=dir([preproc_path filesep 'ICAcleaned_eblock_ft_*.mat']);

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

res_mat=[];
redo=1; complete=0;

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
    
    if redo==1 || exist([preproc_path filesep SubID '_FFT_perBlock_byElec_ICAcleaned_v2.mat'])==0
        subj_pow=[];
        %         subj_SNR=[];
        %         subj_SNRtag=[];
        
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
        if size(data.trial,2)<9
            continue;
        end
        for nBl=1:9 %size(data.trial,2) %DP capping at 9 blocks
            
            temp_data=data.trial{nBl}(1:58,data.time{nBl}>0);
            temp_data=temp_data-repmat(mean(temp_data,2),1,size(temp_data,2));
            
            
            for nEl=1:58
                fprintf('\b\b\b\b\b\b%2.0f-%2.0f\n',nBl,nEl)
                block_data=temp_data(nEl,:);
                
                %                 w_window=10*data_clean.fsample;
                %                 w_overlap=w_window/2;
                %                 df=0.1;
                %                 freqV=2:0.1:40;
                %                 [pow,faxis] = pwelch(block_data,w_window,w_overlap,freqV,data_clean.fsample,'psd');
                
                param=[];
                
                param.method='fft';
                param.w_window=10*data.fsample;
                param.w_overlap=param.w_window/2;
                param.df=0.1;
                param.freqV=2:0.1:30;
                param.mindist=0.5;
                
                
%                 addpath(([path_fieldtrip filesep '/external/signal/dpss_hack']));
                block_data_win=[];
                for start=1:param.w_window/2:length(block_data)-param.w_window
                    block_data_win=[block_data_win ; block_data(start:start+param.w_window)-mean(block_data(start:start+param.w_window))];
                end
                block_data_win_cleaned=block_data_win(max(abs(block_data_win),[],2)<absThr,:);
%                 block_data_cleaned=reshape(block_data_win_cleaned',1,numel(block_data_win_cleaned));
                if isempty(block_data_win_cleaned)==0
                [logSNR, faxis, logpow]=get_logSNR(block_data_win_cleaned,data.fsample,param);
                logpow=mean(logpow,1);
                                
                Frac = amri_sig_fractal(block_data_win',data.fsample,'detrend',1,'frange',[param.freqV(1) param.freqV(end)]);
                
                subj_logpow(nBl,nEl,:)=logpow(:,faxis>=min(param.freqV) & faxis<=max(param.freqV));
                subj_logsnr(nBl,nEl,:)=logSNR(:,faxis>=min(param.freqV) & faxis<=max(param.freqV));
                ssvep_freq=faxis(faxis>=min(param.freqV) & faxis<=max(param.freqV));
                else
                    subj_logpow(nBl,nEl,:)=nan(1,280);
                    subj_logsnr(nBl,nEl,:)=nan(1,280);
                end
            end
        end
        
        save([preproc_path filesep SubID '_SSVEP_perBlock_byElec_ICAcleaned_v2.mat'],'subj_logpow','newlabels','ssvep_freq','subj_logsnr');
    else
        load([preproc_path filesep SubID '_SSVEP_perBlock_byElec_ICAcleaned_v2.mat']);
    end
   
    nFc=nFc+1;
    all_logpow(nFc,:,:,:)=subj_logpow;
    all_logsnr(nFc,:,:,:)=subj_logsnr;


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
    %     save([powerspec_path filesep 'cohorts_FFT_perBlock_byElec_ICAcleaned.mat'],'YoungM_pow','YoungM_SNR','YoungM_SNRtag','YoungT_pow','YoungT_SNR','YoungT_SNRtag','YoungHN_pow','YoungHN_SNR','YoungHN_SNRtag',...
    %          'OldA_pow','OldA_SNR','OldA_SNRtag','OldHN_pow','OldHN_SNR','OldHN_SNRtag');
    
end
save([preproc_path filesep 'all_SSVEP_perBlock_byElec_ICAcleaned_v2'],'all_logpow','all_logsnr','all_group','all_agegroup','newlabels','ssvep_freq');

%% Topographies
cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.channel=newlabels;
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);
correspCh=[];
for nCh=1:length(layout.label)-2
    correspCh(nCh)=match_str(newlabels,layout.label(nCh));
end


%%
cmap2=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap2=flipud(cmap2);

figure; set(gcf,'Position',[213         173        1027/4         805/3]);
% subplot(1,4,1);
% jbfill([24.5 25.5],[-.7 -.7],[2 2],[50,205,50]/256,[50,205,50]/256,1,0.2);
format_fig;
hold on;
simpleTplot(ssvep_freq,squeeze(nanmean(all_logpow(:,:,match_str(newlabels,'Oz'),:),2)),0,'k',[0 0.05 0.0001 1000],'-',0.1,1,0,1,1);
xlim([2 30])
% ylim([-.7 2])
xlabel('Frequency (Hz)')
ylabel('Power')


figure; set(gcf,'Position',[213         173        1027/4         805/3]);
% subplot(1,4,1);
% jbfill([24.5 25.5],[-.7 -.7],[2 2],[50,205,50]/256,[50,205,50]/256,1,0.2);
format_fig;
hold on;
simpleTplot(ssvep_freq,squeeze(nanmean(all_logsnr(:,:,match_str(newlabels,'Oz'),:),2)),0,'k',[0 0.05 0.0001 1000],'-',0.1,1,0,1,1);
xlim([2 30])
% ylim([-.7 2])
xlabel('Frequency (Hz)')
ylabel('SNR')