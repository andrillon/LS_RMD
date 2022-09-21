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

nF=1;
    file_name = files(nF).name;
    folder_name = files(nF).folder;
    SubID=file_name(1:end-4);
    seps=findstr(SubID,'ft_');
    SubID=SubID(seps(1)+3:end);
    tic;
    fprintf('... working on %s (%g/%g)\n',SubID,nF,length(files))
    
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
    
    
load([preproc_path filesep 'all_FFT_perBlock_byElec_ICAcleaned_v2.mat']);
% all_aperiodic_component should be offset then slope
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


%%
cmap2=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap2=flipud(cmap2);
faxis=fractal.freq;
% subplot(1,4,1);   
% jbfill([24.5 25.5],[-.7 -.7],[2 2],[50,205,50]/256,[50,205,50]/256,1,0.2);

channels_to_plot={'Fz','Cz','Pz','Oz'};
figure; set(gcf,'Position',[ 2104         115         788         574]);
format_fig;
cmap=cbrewer('seq','OrRd',5);
for nCh=1:length(channels_to_plot)
simpleTplot(faxis,squeeze(nanmean(log10(all_pow(:,match_str(newlabels,channels_to_plot{nCh}),:)),2)),0,cmap(nCh+1,:),[0],'-',0.1,1,0,1,1);
hold on;
end
xlim([2 30])
% ylim([-.7 2])
xlabel('Frequency (Hz)')
ylabel('Power')

figure; set(gcf,'Position',[ 2104         115         788         574]);
format_fig;
cmap=cbrewer('seq','OrRd',5);
for nCh=1:length(channels_to_plot)
hold on;
simpleTplot(faxis,squeeze(nanmean(log10(all_pow(all_group==4,match_str(newlabels,channels_to_plot{nCh}),:)),2)),0,cmap(nCh+1,:),[0],'-',0.1,1,0,1,1);
simpleTplot(faxis,squeeze(nanmean(log10(all_pow(all_group==5,match_str(newlabels,channels_to_plot{nCh}),:)),2)),0,cmap(nCh+1,:),[0],':',0.1,1,0,1,1);
end
xlim([2 30])
% ylim([-.7 2])
xlabel('Frequency (Hz)')
ylabel('Power')

%% Aperiodic activity
figure; set(gcf,'Position',[ 2104         115         788         574]);
format_fig;
cmap=cbrewer('seq','OrRd',5);
for nCh=1:length(channels_to_plot)
simpleTplot(faxis,squeeze(nanmean(log10(all_frac(:,match_str(newlabels,channels_to_plot{nCh}),:)),2)),0,cmap(nCh+1,:),[0],'-',0.1,1,0,1,1);
hold on;
end
xlim([2 30])
% ylim([-.7 2])
xlabel('Frequency (Hz)')
ylabel('Power')

figure; set(gcf,'Position',[ 2104         115         788         574]);
format_fig;
cmap=cbrewer('seq','OrRd',5);
for nCh=1:length(channels_to_plot)
hold on;
simpleTplot(faxis,squeeze(nanmean(log10(all_frac(all_group==4,match_str(newlabels,channels_to_plot{nCh}),:)),2)),0,cmap(nCh+1,:),[0],'-',0.1,1,0,1,1);
simpleTplot(faxis,squeeze(nanmean(log10(all_frac(all_group==5,match_str(newlabels,channels_to_plot{nCh}),:)),2)),0,cmap(nCh+1,:),[0],':',0.1,1,0,1,1);
end
xlim([2 30])
% ylim([-.7 2])
xlabel('Frequency (Hz)')
ylabel('Power')


%% Slope
figure; set(gcf,'Position',[213         173        1027         805/3]);
cmap=cbrewer('seq','YlOrRd',64); % select a sequential colorscale from yellow to red (64)
cmap(cmap<0)=0;

for nD=1:2
    subplot(1,2,nD); format_fig; %On nD==2 this line deletes prev graph & starts new figure
    temp_topo=[];
    for nCh=1:length(layout.label)-2
     temp_topo(nCh)=nanmean(all_aperiodic_component(all_aperiodic_component(:,2)==nD+3 & all_aperiodic_component(:,4)==match_str(newlabels,layout.label{nCh}),end)); %Slope
    end
    simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
    colormap(cmap);
%     if nD==1
        
%     end
maxmin(nD,1)=min((temp_topo));
maxmin(nD,2)=max((temp_topo));
end
caxis([min(maxmin(:,1)) max(maxmin(:,2))])
hb=colorbar('Position',[0.9195    0.6373    0.0143    0.2881]);
% print([powerspec_path filesep 'Topo_Alpha_v5.eps'],'-dpng', '-r300');
title('Slope Aperiodic Young 90% vs Old 90%')
%% Offset
figure; set(gcf,'Position',[213         173        1027         805/3]);
cmap=cbrewer('seq','YlOrRd',64); % select a sequential colorscale from yellow to red (64)
cmap(cmap<0)=0;

for nD=1:2
    subplot(1,2,nD); format_fig; %On nD==2 this line deletes prev graph & starts new figure
    temp_topo=[];
    for nCh=1:length(layout.label)-2
     temp_topo(nCh)=nanmean(all_aperiodic_component(all_aperiodic_component(:,2)==nD+3 & all_aperiodic_component(:,4)==match_str(newlabels,layout.label{nCh}),end-1)); %Slope
    end
    simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
    colormap(cmap);
%     if nD==1
        
%     end
maxmin(nD,1)=min((temp_topo));
maxmin(nD,2)=max((temp_topo));
end
caxis([min(maxmin(:,1)) max(maxmin(:,2))])
hb=colorbar('Position',[0.9195    0.6373    0.0143    0.2881]);
% print([powerspec_path filesep 'Topo_Alpha_v5.eps'],'-dpng', '-r300');
title('Offset Aperiodic Young 90% vs Old 90%')

%% Alpha Peak Frequency
all_alphapeak=all_peak_component(all_peak_component(:,5)>8 & all_peak_component(:,5)<13,:);

figure; set(gcf,'Position',[213         173        1027         805/3]);
cmap=cbrewer('seq','YlOrRd',64); % select a sequential colorscale from yellow to red (64)
cmap(cmap<0)=0;

for nD=1:2
    subplot(1,2,nD); format_fig; %On nD==2 this line deletes prev graph & starts new figure
    temp_topo=nan(1,length(layout.label)-2);
    for nCh=1:length(layout.label)-2
     temp_topo(nCh)=nanmean(all_alphapeak(all_alphapeak(:,2)==nD+3 & all_alphapeak(:,4)==match_str(newlabels,layout.label{nCh}),5)); %Slope
    end
    simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
    colormap(cmap);
%     if nD==1
        
%     end
maxmin(nD,1)=min((temp_topo));
maxmin(nD,2)=max((temp_topo));
end
caxis([min(maxmin(:,1)) max(maxmin(:,2))])
hb=colorbar('Position',[0.9195    0.6373    0.0143    0.2881]);
% print([powerspec_path filesep 'Topo_Alpha_v5.eps'],'-dpng', '-r300');
title('Alpha Peak Frequency Young 90% vs Old 90%')

%% Alpha Peak Amplitude
all_alphapeak=all_peak_component(all_peak_component(:,5)>8 & all_peak_component(:,5)<13,:);

figure; set(gcf,'Position',[213         173        1027         805/3]);
cmap=cbrewer('seq','YlOrRd',64); % select a sequential colorscale from yellow to red (64)
cmap(cmap<0)=0;

for nD=1:2
    subplot(1,2,nD); format_fig; %On nD==2 this line deletes prev graph & starts new figure
    temp_topo=nan(1,length(layout.label)-2);
    for nCh=1:length(layout.label)-2
     temp_topo(nCh)=nanmean(all_alphapeak(all_alphapeak(:,2)==nD+3 & all_alphapeak(:,4)==match_str(newlabels,layout.label{nCh}),6)); %
    end
    simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
    colormap(cmap);
%     if nD==1
        
%     end
maxmin(nD,1)=min((temp_topo));
maxmin(nD,2)=max((temp_topo));
end
caxis([min(maxmin(:,1)) max(maxmin(:,2))])
hb=colorbar('Position',[0.9195    0.6373    0.0143    0.2881]);
% print([powerspec_path filesep 'Topo_Alpha_v5.eps'],'-dpng', '-r300');
title('Alpha Peak Amplitude Young 90% vs Old 90%')

%% Beta Peak Amplitude
all_betapeak=all_peak_component(all_peak_component(:,5)>15 & all_peak_component(:,5)<20,:);

figure; set(gcf,'Position',[213         173        1027         805/3]);
cmap=cbrewer('seq','YlOrRd',64); % select a sequential colorscale from yellow to red (64)
cmap(cmap<0)=0;

for nD=1:2
    subplot(1,2,nD); format_fig; %On nD==2 this line deletes prev graph & starts new figure
    temp_topo=nan(1,length(layout.label)-2);
    for nCh=1:length(layout.label)-2
     temp_topo(nCh)=nanmean(all_betapeak(all_betapeak(:,2)==nD+3 & all_betapeak(:,4)==match_str(newlabels,layout.label{nCh}),6)); %
    end
    simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
    colormap(cmap);
%     if nD==1
        
%     end
maxmin(nD,1)=min((temp_topo));
maxmin(nD,2)=max((temp_topo));
end
caxis([min(maxmin(:,1)) max(maxmin(:,2))])
hb=colorbar('Position',[0.9195    0.6373    0.0143    0.2881]);
% print([powerspec_path filesep 'Topo_Alpha_v5.eps'],'-dpng', '-r300');
title('Alpha Peak Amplitude Young 90% vs Old 90%')

