%%
clear all;
close all;

%%
run LS_RMD_localdef.m
addpath((path_fieldtrip));
ft_defaults;
addpath(genpath(path_LSCPtools));

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
myfreqs=1:0.1:40;
redo=1; complete=0;

% m = 1; t = 1; h = 1; a = 1; hn = 1;
% 
% YoungM_pow=[]; YoungM_SNR=[]; YoungM_SNRtag=[]; YoungT_pow=[]; YoungT_SNR=[]; YoungT_SNRtag=[]; YoungHN_pow=[]; YoungHN_SNR=[]; YoungHN_SNRtag=[];
% OldA_pow=[]; OldA_SNR=[]; OldA_SNRtag=[]; OldHN_pow=[]; OldHN_SNR=[]; OldHN_SNRtag=[];

% for nF=1:length(files)
nFc=0;
for nF=1:5
    file_name = files(nF).name;
    folder_name = files(nF).folder;
    SubID=file_name(1:end-4);
    seps=findstr(SubID,'ft_');
    SubID=SubID(seps(1)+3:end);
    tic;
    fprintf('... working on %s (%g/%g)\n',SubID,nF,length(files))
    
    if redo==1 || exist([preproc_path filesep SubID '_FFT_perBlock_byElec_ICAcleaned.mat'])==0
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
        for nBl=1:9 %size(data.trial,2)
            
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
                
                param.method='welch';
                param.w_window=10*data.fsample;
                param.w_overlap=param.w_window/2;
                param.df=0.1;
                param.freqV=2:0.1:40;
                param.mindist=0.5;
                
                [logSNR, faxis, logpow]=get_logSNR(block_data,data.fsample,param);
                
                subj_pow(nBl,nEl,:)=logpow(faxis<40);
%                 subj_SNR(nBl,nEl,:)=logSNR(faxis<40);
%                 [~,closestidx]=findclosest(faxis,25);
%                 subj_SNRtag(nBl,nEl)=logSNR(closestidx);
                faxis=faxis(faxis<40);
            end
        end
               
        save([preproc_path filesep SubID '_FFT_perBlock_byElec_ICAcleaned.mat'],'subj_pow','newlabels','faxis');
             
    end
    
    %DP - in case I want to split by cohort
%     if ismember(SubID,CohYoungM)
%         YoungM_pow(m,:,:,:)=subj_pow;
%         YoungM_SNR(m,:,:,:)=subj_SNR;
%         YoungM_SNRtag(m,:,:)=subj_SNRtag;
%         m=m+1;
%     elseif ismember(SubID,CohYoungT)
%         YoungT_pow(t,:,:,:)=subj_pow;
%         YoungT_SNR(t,:,:,:)=subj_SNR;
%         YoungT_SNRtag(t,:,:)=subj_SNRtag;
%         t=t+1;
%     elseif ismember(SubID,CohYoungHN)
%         YoungHN_pow(h,:,:,:)=subj_pow;
%         YoungHN_SNR(h,:,:,:)=subj_SNR;
%         YoungHN_SNRtag(h,:,:)=subj_SNRtag;
%         h=h+1;
%     elseif ismember(SubID,CohOldA)
%         OldA_pow(a,:,:,:)=subj_pow;
%         OldA_SNR(a,:,:,:)=subj_SNR;
%         OldA_SNRtag(a,:,:)=subj_SNRtag;
%         a=a+1;
%     elseif ismember(SubID,CohOldHN)
%         OldHN_pow(hn,:,:,:)=subj_pow;
%         OldHN_SNR(hn,:,:,:)=subj_SNR;
%         OldHN_SNRtag(hn,:,:)=subj_SNRtag;
%         hn=hn+1;
%     end
    nFc=nFc+1;
    all_pow(nFc,:,:,:)=subj_pow;
%     all_SNR(nF,:,:,:)=subj_SNR;
%     all_SNRtag(nF,:,:)=subj_SNRtag;


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
    save([powerspec_path filesep 'all_FFT_perBlock_byElec_ICAcleaned.mat'],'all_pow','all_group','all_agegroup');

%% Topographies 25Hz tag  - need to check channel layout etc
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
jbfill([24.5 25.5],[-.7 -.7],[2 2],[50,205,50]/256,[50,205,50]/256,1,0.2);
format_fig;
hold on;
simpleTplot(faxis',squeeze(nanmean(all_pow(:,:,match_str(chLabels,'Cz'),:),3)),0,Colors(nD,:),[0 0.05 0.0001 1000],'-',0.1,1,0,1,1);
xlim([2 30])
ylim([-.7 2])
xlabel('Frequency (Hz)')
ylabel('Power')
% print([powerspec_path filesep 'Topo_FreqTag_Clusters_byFreq_v5.eps'],'-dpng', '-r300');


%% Alpha 8-11Hz - needs adapting of nD etc
% fix channel names
    mylabels=data.label';
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
    
figure; set(gcf,'Position',[213         173        1027         805/3]);
cmap=cbrewer('seq','YlOrRd',64); % select a sequential colorscale from yellow to red (64)
cmap(cmap<0)=0;
for nD=1:2
    subplot(1,2,nD); format_fig;
    temp_topo=squeeze(nanmean(nanmean(nanmean(all_pow(all_agegroup==nD-1,:,correspCh,faxis>8 & faxis<11),1),2),4));
    simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
    colormap(cmap);
    if nD==4
        hb=colorbar('Position',[0.9195    0.6373    0.0143    0.2881]);
    end
    caxis([-2.5 -1]);
%     title(ColorsDlabels{nD});
end
% print([powerspec_path filesep 'Topo_Alpha_v5.eps'],'-dpng', '-r300');

%%
% % alphaFreqs=find((faxis>8 & faxis<8.6) | (faxis>8.95 & faxis<9.8) | (faxis>10.125 & faxis<11));
% alphaFreqs=find((faxis>8 & faxis<11));
% figure; set(gcf,'Position',[213         173        1027/4         805/3]);
% % subplot(1,4,1); 
% format_fig;
% jbfill([8 11],[-5 -5],[-.5 -.5],[50,205,50]/256,[50,205,50]/256,1,0.2);
% hold on;
% for nD=1:4
%     simpleTplot(faxis',squeeze(nanmean(all_pow(:,nD,:,match_str(chLabels,'Cz'),:),3)),0,Colors(nD,:),[0 0.05 0.0001 1000],'-',0.1,1,0,1,1);
% end
% xlim([2 30])
% ylim([[-5 -0.5]])
% xlabel('Frequency (Hz)')
% ylabel('Power (dB)')
% print([powerspec_path filesep 'Topo_Alpha_Clusters_byFreq_v5.eps'],'-dpng', '-r300');
% 
% figure; set(gcf,'Position',[1 1 880 880]);
% winTime=[0.05 0.3];
% for nD=1:size(PosDrugs,1)
%     for nD2=1:size(PosDrugs,2)
%         if isempty(PosDrugs{nD,nD2})
%             continue;
%         end
%     subplot(3,3,3*(nD-1)+(nD2)); format_fig;
%     temp_topo=[];
%     for nCh=1:length(layout.label)-2
%         temp1=squeeze(nanmean(nanmean(all_pow(:,PosDrugs{nD,nD2}(1),:,match_str(chLabels,layout.label(nCh)),alphaFreqs),3),5));
%         temp0=squeeze(nanmean(nanmean(all_pow(:,PosDrugs{nD,nD2}(2),:,match_str(chLabels,layout.label(nCh)),alphaFreqs),3),5));
%         [h, pV, ~ , stats]=ttest(temp1,temp0);
%         temp_topo(nCh)=stats.tstat;%-...
%         %             squeeze(nanmean(nanmean(nanmean(all_ERP_NT_offset(:,nD,match_str(chLabels,layout.label(nCh)),xTime_offset>winTime(1) & xTime_offset<winTime(2)),1),2),4));
%     end
%     temp_topo1=squeeze(nanmean(nanmean(all_pow(:,PosDrugs{nD,nD2}(2),:,correspCh,alphaFreqs),3),5));
%     temp_topo2=squeeze(nanmean(nanmean(all_pow(:,PosDrugs{nD,nD2}(1),:,correspCh,alphaFreqs),3),5));
%     [stat] = compute_Topo_clusterPerm_v2(temp_topo1,temp_topo2,0,chLabels(correspCh),data_clean.fsample,0.05,0.05,1000,layout);
%     
%     
%     simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
%     colormap(cmap2);
%     
%     if isfield(stat,'posclusters') && ~isempty(stat.posclusters)
%         sigClusters=find([stat.posclusters.prob]<0.05);
%         for k=1:length(sigClusters)
%             ft_plot_lay_me(layout, 'chanindx',find(stat.posclusterslabelmat==sigClusters(k)),'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
%         end
%     end
%     if isfield(stat,'negclusters') && ~isempty(stat.negclusters)
%         sigClusters=find([stat.negclusters.prob]<0.05);
%         for k=1:length(sigClusters)
%             ft_plot_lay_me(layout, 'chanindx',find(stat.negclusterslabelmat==sigClusters(k)),'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
%         end
%     end
%     if nD==4
%         hb=colorbar('Position',[0.9195    0.6373    0.0143    0.2881]);
%     end
%     caxis([-1 1]*5);
%     title(sprintf('%s vs %s',ColorsDlabels{PosDrugs{nD,nD2}(1)},ColorsDlabels{PosDrugs{nD,nD2}(2)}));
%     end
% end
% 
% print([powerspec_path filesep 'Topo_Alpha_Clusters_v5.eps'],'-dpng', '-r300');

%%
% figure; set(gcf,'Position',[213         173        1027/4         805/3]);
% 
% print('-dpng', '-r300', '../../Figures/Cz_PowSpec_v5.png')
% 
% figure; set(gcf,'Position',[213         173        1027/4         805/3]);
% 
% % ylim([[-5 -0.5]])