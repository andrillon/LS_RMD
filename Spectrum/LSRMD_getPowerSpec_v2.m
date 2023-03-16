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
all_peak_component=[];
all_aperiodic_component=[];
all_detectedpeaks=[];
if redo==1
    for nF=1:length(files)
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
        
        fprintf('%2.0f-%2.0f\n',0,0)
        if size(data.trial,2)<8
            continue;
        end
        
        
        % chunk into 10-second segments
        cfg               = [];
        cfg.length        = 10;
        cfg.overlap       = 0.5;
        data2              = ft_redefinetrial(cfg, data);
        
        %%%%% TRY TO REMOVE ARTEFACTS?
        maxSignal=[];
        for nTr=1:length(data2.trial)
            maxSignal(nTr,:)=max(abs(data2.trial{nTr}),[],2);
        end
        cfg               = [];
        absThr=250;
        pptionThr=0.4;
        fprintf('... ... discarding %g trials (/%g) because of %g channels above threshold of %g\n',...
            sum(mean(maxSignal>absThr,2)>=pptionThr),length(data2.trial),100*pptionThr,absThr)
        cfg.trials        = find(mean(maxSignal>absThr,2)<pptionThr);
        data2              = ft_redefinetrial(cfg, data2);
        
        % compute the fractal and original spectra
        cfg               = [];
        cfg.foilim        = [2 30];
        cfg.pad           = 'nextpow2';
        cfg.tapsmofrq     = 0.5;
        cfg.method        = 'mtmfft';
        cfg.output        = 'fooof_aperiodic';
        cfg.fooof.max_peaks = 4;
        cfg.fooof.proximity_threshold = 1;
        fractal = ft_freqanalysis(cfg, data2);
        cfg.output        = 'pow';
        pow = ft_freqanalysis(cfg, data2);
        cfg.output        = 'fooof_peaks';
        pow_peaks = ft_freqanalysis(cfg, data2);
        
        nFc=nFc+1;
        all_pow(nFc,:,:,:)=log10(pow.powspctrm);
        %         if exist('subj_osci')~=0
        %             all_osci(nFc,:,:,:)=subj_osci;
        %         end
        %         if exist('subj_mixed')~=0
        %             all_mixed(nFc,:,:,:)=subj_mixed;
        %         end
        %     if exist('subj_frac')~=0
        all_frac(nFc,:,:,:)=fractal.powspctrm;
        all_osci(nFc,:,:,:)=pow_peaks.powspctrm;
        %     end
        
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
        
        table_fractal=fractal.fooofparams;
        for k=1:length(fractal.label)
            all_aperiodic_component=[all_aperiodic_component ; [nFc all_group(nFc) all_agegroup(nFc) k fractal.fooofparams(k).aperiodic_params]];
            all_peak_component=[all_peak_component ; [repmat([nFc all_group(nFc) all_agegroup(nFc) k],size(fractal.fooofparams(k).peak_params,1),1) fractal.fooofparams(k).peak_params]];
            
            if ~isempty(fractal.fooofparams(k).peak_params)
                is_alphapeak=sum(fractal.fooofparams(k).peak_params(:,1)>8 & fractal.fooofparams(k).peak_params(:,1)<12);
                is_slowpeak=sum(fractal.fooofparams(k).peak_params(:,1)<7);
                
                all_detectedpeaks=[all_detectedpeaks ; [nFc all_group(nFc) all_agegroup(nFc) k is_alphapeak is_slowpeak]];
            else
                all_detectedpeaks=[all_detectedpeaks ; [nFc all_group(nFc) all_agegroup(nFc) k 0 0]];
            end
        end
        
        
    end
    save([preproc_path filesep 'all_FFT_perBlock_byElec_ICAcleaned_v2.mat'],'all_pow','all_group','all_agegroup','all_peak_component','all_aperiodic_component','all_frac','fractal','newlabels','all_osci');
    
else
    load([preproc_path filesep 'all_FFT_perBlock_byElec_ICAcleaned_v2.mat']);
end

%% Topographies 25Hz tag
cfg = [];
cfg.layout = 'acticap-64ch-standard2.mat';
cfg.channel=newlabels;
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);
correspCh=[];
correspCh2=[];
for nCh=1:length(layout.label)-2
    correspCh(nCh)=match_str(newlabels,layout.label(nCh));
    correspCh2(nCh)=match_str(newlabels,fractal.label(nCh));
end


%%
cmap2=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap2=flipud(cmap2);
faxis=fractal.freq;
figure; set(gcf,'Position',[213         173        1027/4         805/3]);
% subplot(1,4,1);
% jbfill([24.5 25.5],[-.7 -.7],[2 2],[50,205,50]/256,[50,205,50]/256,1,0.2);
format_fig;
hold on;
% simpleTplot(faxis,squeeze(nanmean(all_pow(:,match_str(fractal.label,'Cz'),:),2)),0,'k',[0 0.05 0.0001 1000],'-',0.1,1,0,1,1); % DP - original, don't think use of fractal correct
% simpleTplot(faxis,squeeze(nanmean(all_pow(:,match_str(layout.label,'Cz'),:),2)),0,'k',[0 0.05 0.0001 1000],'-',0.1,1,0,1,1);
simpleTplot(faxis,squeeze(nanmean(all_osci(:,match_str(layout.label,'Cz'),:),2)),0,'k',[0 0.05 0.0001 1000],'-',0.1,1,0,1,1);

xlim([2 30])
% ylim([-.7 2])
xlabel('Frequency (Hz)')
ylabel('Power')

%% Delta 1-3Hz

figure; set(gcf,'Position',[213         173        1027         805/3]);
cmap=cbrewer('seq','YlOrRd',64); % select a sequential colorscale from yellow to red (64)
cmap(cmap<0)=0;

for nD=1:2
    subplot(1,2,nD); format_fig; %On nD==2 this line deletes prev graph & starts new figure
    temp_topo=squeeze(nanmean(nanmean(all_pow(all_agegroup==nD-1,correspCh2,faxis>1 & faxis<3),1),3)); %Delta
    %     temp_topo=squeeze(nanmean(nanmean(all_osci(all_agegroup==nD-1,correspCh2,faxis>1 & faxis<3),1),3)); %Delta
    simpleTopoPlot_ft(temp_topo(correspCh), layout,'on',[],0,1);
    %     simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1); %DP - original, think wrong
    colormap(cmap);
    if nD==1
        hb=colorbar('Position',[0.9195    0.6373    0.0143    0.2881]);
        title('Delta Young')
    else
        title('Delta Old')
    end
    maxmin(nD,1)=min((temp_topo));
    maxmin(nD,2)=max((temp_topo));
end
caxis([min(maxmin(:,1)) max(maxmin(:,2))])

figure; format_fig;
youngtopo=squeeze(nanmean(nanmean(all_pow(all_agegroup==0,correspCh2,faxis>1 & faxis<3),1),3));
oldtopo=squeeze(nanmean(nanmean(all_pow(all_agegroup==1,correspCh2,faxis>1 & faxis<3),1),3));
% youngtopo=squeeze(nanmean(nanmean(all_osci(all_agegroup==0,correspCh2,faxis>1 & faxis<3),1),3));
% oldtopo=squeeze(nanmean(nanmean(all_osci(all_agegroup==1,correspCh2,faxis>1 & faxis<3),1),3));
difftopo=oldtopo-youngtopo;
simpleTopoPlot_ft(difftopo(correspCh), layout,'on',[],0,1);
%     simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1); %DP - original, think wrong
colormap(cmap);
hb=colorbar('Position',[0.9195    0.6373    0.0143    0.2881]);

%% Theta 4-7Hz

figure; set(gcf,'Position',[213         173        1027         805/3]);
cmap=cbrewer('seq','YlOrRd',64); % select a sequential colorscale from yellow to red (64)
cmap(cmap<0)=0;

for nD=1:2
    subplot(1,2,nD); format_fig; %On nD==2 this line deletes prev graph & starts new figure
    temp_topo=squeeze(nanmean(nanmean(all_pow(all_agegroup==nD-1,correspCh2,faxis>4 & faxis<7),1),3)); %Theta
    %     temp_topo=squeeze(nanmean(nanmean(all_osci(all_agegroup==nD-1,correspCh2,faxis>4 & faxis<7),1),3)); %Theta
    simpleTopoPlot_ft(temp_topo(correspCh), layout,'on',[],0,1);
    %     simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1); %DP - original, think wrong
    colormap(cmap);
    if nD==1
        hb=colorbar('Position',[0.9195    0.6373    0.0143    0.2881]);
        title('Theta Young')
    else
        title('Theta Old')
    end
    maxmin(nD,1)=min((temp_topo));
    maxmin(nD,2)=max((temp_topo));
end
caxis([min(maxmin(:,1)) max(maxmin(:,2))])

figure; format_fig;
youngtopo=squeeze(nanmean(nanmean(all_pow(all_agegroup==0,correspCh2,faxis>4 & faxis<7),1),3));
oldtopo=squeeze(nanmean(nanmean(all_pow(all_agegroup==1,correspCh2,faxis>4 & faxis<7),1),3));
% youngtopo=squeeze(nanmean(nanmean(all_osci(all_agegroup==0,correspCh2,faxis>4 & faxis<7),1),3));
% oldtopo=squeeze(nanmean(nanmean(all_osci(all_agegroup==1,correspCh2,faxis>4 & faxis<7),1),3));
difftopo=oldtopo-youngtopo;
simpleTopoPlot_ft(difftopo(correspCh), layout,'on',[],0,1);
%     simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1); %DP - original, think wrong
colormap(cmap);
hb=colorbar('Position',[0.9195    0.6373    0.0143    0.2881]);

%% Alpha 8-11Hz

figure; set(gcf,'Position',[213         173        1027         805/3]);
cmap=cbrewer('seq','YlOrRd',64); % select a sequential colorscale from yellow to red (64)
cmap(cmap<0)=0;

for nD=1:2
    subplot(1,2,nD); format_fig; %On nD==2 this line deletes prev graph & starts new figure
    temp_topo=squeeze(nanmean(nanmean(all_pow(all_agegroup==nD-1,correspCh2,faxis>8 & faxis<11),1),3)); %Alpha
    %     temp_topo=squeeze(nanmean(nanmean(all_osci(all_agegroup==nD-1,correspCh2,faxis>8 & faxis<11),1),3)); %Alpha
    simpleTopoPlot_ft(temp_topo(correspCh), layout,'on',[],0,1);
    %     simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1); %DP - original, think wrong
    colormap(cmap);
    if nD==1
        hb=colorbar('Position',[0.9195    0.6373    0.0143    0.2881]);
        title('Alpha Young')
    else
        title('Alpha Old')
    end
    maxmin(nD,1)=min((temp_topo));
    maxmin(nD,2)=max((temp_topo));
end
caxis([min(maxmin(:,1)) max(maxmin(:,2))])

figure; format_fig;
youngtopo=squeeze(nanmean(nanmean(all_pow(all_agegroup==0,correspCh2,faxis>8 & faxis<11),1),3));
oldtopo=squeeze(nanmean(nanmean(all_pow(all_agegroup==1,correspCh2,faxis>8 & faxis<11),1),3));
% youngtopo=squeeze(nanmean(nanmean(all_osci(all_agegroup==0,correspCh2,faxis>8 & faxis<11),1),3));
% oldtopo=squeeze(nanmean(nanmean(all_osci(all_agegroup==1,correspCh2,faxis>8 & faxis<11),1),3));
difftopo=oldtopo-youngtopo;
simpleTopoPlot_ft(difftopo(correspCh), layout,'on',[],0,1);
%     simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1); %DP - original, think wrong
colormap(cmap);
hb=colorbar('Position',[0.9195    0.6373    0.0143    0.2881]);

%% Beta 12-29Hz

figure; set(gcf,'Position',[213         173        1027         805/3]);
cmap=cbrewer('seq','YlOrRd',64); % select a sequential colorscale from yellow to red (64)
cmap(cmap<0)=0;

for nD=1:2
    subplot(1,2,nD); format_fig; %On nD==2 this line deletes prev graph & starts new figure
    temp_topo=squeeze(nanmean(nanmean(all_pow(all_agegroup==nD-1,correspCh2,faxis>12 & faxis<29),1),3)); %Beta
    %     temp_topo=squeeze(nanmean(nanmean(all_osci(all_agegroup==nD-1,correspCh2,faxis>12 & faxis<29),1),3)); %Beta
    simpleTopoPlot_ft(temp_topo(correspCh), layout,'on',[],0,1);
    %     simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1); %DP - original, think wrong
    colormap(cmap);
    if nD==1
        hb=colorbar('Position',[0.9195    0.6373    0.0143    0.2881]);
        title('Beta Young')
    else
        title('Beta Old')
    end
    maxmin(nD,1)=min((temp_topo));
    maxmin(nD,2)=max((temp_topo));
end
caxis([min(maxmin(:,1)) max(maxmin(:,2))])

figure; format_fig;
youngtopo=squeeze(nanmean(nanmean(all_pow(all_agegroup==0,correspCh2,faxis>12 & faxis<29),1),3));
oldtopo=squeeze(nanmean(nanmean(all_pow(all_agegroup==1,correspCh2,faxis>12 & faxis<29),1),3));
% youngtopo=squeeze(nanmean(nanmean(all_osci(all_agegroup==0,correspCh2,faxis>12 & faxis<29),1),3));
% oldtopo=squeeze(nanmean(nanmean(all_osci(all_agegroup==1,correspCh2,faxis>12 & faxis<29),1),3));
difftopo=oldtopo-youngtopo;
simpleTopoPlot_ft(difftopo(correspCh), layout,'on',[],0,1);
%     simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1); %DP - original, think wrong
colormap(cmap);
hb=colorbar('Position',[0.9195    0.6373    0.0143    0.2881]);

%% Aperiodic component
figure; set(gcf,'Position',[213         173        1027         805/3]);
cmap=cbrewer('seq','YlOrRd',64); % select a sequential colorscale from yellow to red (64)
cmap(cmap<0)=0;


agegroups={[4],[5]}; % Young vs Old
for nD=1:2
    temp_topo=[];
    for nCh=1:length(fractal.label)
        temp_topo(nCh)=squeeze(nanmean(all_aperiodic_component(all_aperiodic_component(:,2)==agegroups{nD} & all_aperiodic_component(:,4)==nCh,5),1));
    end
    subplot(1,2,nD); format_fig; %On nD==2 this line deletes prev graph & starts new figure
    %     simpleTopoPlot_ft(temp_topo(correspCh2),layout,'on',[],0,1); %DP - this was the original, but seems to plot incorrect layout. If I change to fractal, error
    simpleTopoPlot_ft(temp_topo(correspCh),layout,'on',[],0,1); %DP - looks better but double check
    colormap(cmap);
    %     if nD==1
    %         hb=colorbar('Position',[0.9195    0.6373    0.0143    0.2881]);
    %     end
    colorbar;
    if nD==1
        title('Offset Young')
    elseif nD==2
        title('Offset Old')
    end
    
    maxmin(nD,1)=min((temp_topo));
    maxmin(nD,2)=max((temp_topo));
end
caxis([min(maxmin(:,1)) max(maxmin(:,2))])

figure;
youngtopo=[];oldtopo=[];
for nCh=1:length(fractal.label)
    youngtopo(nCh)=squeeze(nanmean(all_aperiodic_component(all_aperiodic_component(:,2)==agegroups{1} & all_aperiodic_component(:,4)==nCh,5),1));
    oldtopo(nCh)=squeeze(nanmean(all_aperiodic_component(all_aperiodic_component(:,2)==agegroups{2} & all_aperiodic_component(:,4)==nCh,5),1));
end
difftopo=oldtopo-youngtopo;
simpleTopoPlot_ft(difftopo(correspCh),layout,'on',[],0,1); %DP - looks better but double check
colormap(cmap);
%     if nD==1
%         hb=colorbar('Position',[0.9195    0.6373    0.0143    0.2881]);
%     end
colorbar;

figure;
for nD=1:2
    temp_topo=[];
    for nCh=1:length(fractal.label)
        temp_topo(nCh)=squeeze(nanmean(all_aperiodic_component(all_aperiodic_component(:,2)==agegroups{nD} & all_aperiodic_component(:,4)==nCh,6),1));
    end
    subplot(1,2,nD); format_fig; %On nD==2 this line deletes prev graph & starts new figure
    %     simpleTopoPlot_ft(temp_topo(correspCh2),layout,'on',[],0,1); %DP - this was the original, but seems to plot incorrect layout
    simpleTopoPlot_ft(temp_topo(correspCh),layout,'on',[],0,1); %DP - looks better but double check
    colormap(cmap);
    %     if nD==1
    %         hb=colorbar('Position',[0.9195    0.6373    0.0143    0.2881]);
    %     end
    colorbar;
    if nD==1
        title('Slope Young')
    elseif nD==2
        title('Slope Old')
    end
    maxmin(nD,1)=min((temp_topo));
    maxmin(nD,2)=max((temp_topo));
end
caxis([min(maxmin(:,1)) max(maxmin(:,2))])

figure;
youngtopo=[];oldtopo=[];
for nCh=1:length(fractal.label)
    youngtopo(nCh)=squeeze(nanmean(all_aperiodic_component(all_aperiodic_component(:,2)==agegroups{1} & all_aperiodic_component(:,4)==nCh,6),1));
    oldtopo(nCh)=squeeze(nanmean(all_aperiodic_component(all_aperiodic_component(:,2)==agegroups{2} & all_aperiodic_component(:,4)==nCh,6),1));
end
difftopo=oldtopo-youngtopo;
simpleTopoPlot_ft(difftopo(correspCh),layout,'on',[],0,1); %DP - looks better but double check
colormap(cmap);
%     if nD==1
%         hb=colorbar('Position',[0.9195    0.6373    0.0143    0.2881]);
%     end
colorbar;

%% Analysis of peaks
figure;
myLabels={'Fz','Cz','Pz','Oz'};
for nCh=1:length(myLabels)
    subplot(2,2,nCh)
    hold on;
    for nD=1:2
        histogram(all_peak_component(all_aperiodic_component(:,2)==agegroups{nD} & all_aperiodic_component(:,4)==match_str(fractal.label,myLabels{nCh}),5),2:1:30);
    end
    xlabel('Freq (Hz)')
    ylabel('Count peaks')
    format_fig;
    legend({'Young','Old'})
    title(['Peaks ' myLabels{nCh}])
end

%% Slow peaks (1-7Hz)

all_slowpeak_component=all_peak_component(all_peak_component(:,5)>1 & all_peak_component(:,5)<7,:);

figure;
agegroups={[4],[5]}; % Young vs Old
for nD=1:2
    temp_topo=[];
    for nCh=1:length(fractal.label)
        temp_topo(nCh)=squeeze(nanmean(all_slowpeak_component(all_slowpeak_component(:,2)==agegroups{nD} & all_slowpeak_component(:,4)==nCh,5),1)); %Alpha
    end
    subplot(1,2,nD); format_fig; %On nD==2 this line deletes prev graph & starts new figure
    %     simpleTopoPlot_ft(temp_topo(correspCh2),layout,'on',[],0,1); %DP - this was the original, but seems to plot incorrect layout
    simpleTopoPlot_ft(temp_topo(correspCh),layout,'on',[],0,1); %DP - looks better but double check
    colormap(cmap);
    %     if nD==1
    %         hb=colorbar('Position',[0.9195    0.6373    0.0143    0.2881]);
    %     end
    if nD==1
        title('Slow Freq Young')
    elseif nD==2
        title('Slow Freq Old')
    end
    maxmin(nD,1)=min((temp_topo));
    maxmin(nD,2)=max((temp_topo));
end
caxis([min(maxmin(:,1)) max(maxmin(:,2))])
colorbar;

figure;
for nD=1:2
    temp_topo=[];
    for nCh=1:length(fractal.label)
        temp_topo(nCh)=squeeze(nanmean(all_slowpeak_component(all_slowpeak_component(:,2)==agegroups{nD} & all_slowpeak_component(:,4)==nCh,6),1)); %Alpha
    end
    subplot(1,2,nD); format_fig;
    %     simpleTopoPlot_ft(temp_topo(correspCh2),layout,'on',[],0,1); %DP - this was the original, but seems to plot incorrect layout
    simpleTopoPlot_ft(temp_topo(correspCh),layout,'on',[],0,1); %DP - looks better but double check
    colormap(cmap);
    %     if nD==1
    %         hb=colorbar('Position',[0.9195    0.6373    0.0143    0.2881]);
    %     end
    if nD==1
        title('Slow Peak Amplitude Young')
        colorbar;
    elseif nD==2
        title('Slow Peak Amplitude Old')
    end
    maxmin(nD,1)=min((temp_topo));
    maxmin(nD,2)=max((temp_topo));
end
caxis([min(maxmin(:,1)) max(maxmin(:,2))])

figure; oldtopo=[]; youngtopo=[];
for nCh=1:length(fractal.label)
    youngtopo(nCh)=squeeze(nanmean(all_slowpeak_component(all_slowpeak_component(:,2)==agegroups{1} & all_slowpeak_component(:,4)==nCh,5),1));
    oldtopo(nCh)=squeeze(nanmean(all_slowpeak_component(all_slowpeak_component(:,2)==agegroups{2} & all_slowpeak_component(:,4)==nCh,5),1));
end
difftopo=oldtopo-youngtopo;
simpleTopoPlot_ft(difftopo(correspCh),layout,'on',[],0,1); %DP - looks better but double check
colormap(cmap);
colorbar;

figure; oldtopo=[]; youngtopo=[];
for nCh=1:length(fractal.label)
    youngtopo(nCh)=squeeze(nanmean(all_slowpeak_component(all_slowpeak_component(:,2)==agegroups{1} & all_slowpeak_component(:,4)==nCh,6),1));
    oldtopo(nCh)=squeeze(nanmean(all_slowpeak_component(all_slowpeak_component(:,2)==agegroups{2} & all_slowpeak_component(:,4)==nCh,6),1));
end
difftopo=oldtopo-youngtopo;
simpleTopoPlot_ft(difftopo(correspCh),layout,'on',[],0,1); %DP - looks better but double check
colormap(cmap);
colorbar;

%% Delta peaks (1-3Hz)

all_deltapeak_component=all_peak_component(all_peak_component(:,5)>1 & all_peak_component(:,5)<3,:);

figure;
agegroups={[4],[5]}; % Young vs Old
for nD=1:2
    temp_topo=[];
    for nCh=1:length(fractal.label)
        temp_topo(nCh)=squeeze(nanmean(all_deltapeak_component(all_deltapeak_component(:,2)==agegroups{nD} & all_deltapeak_component(:,4)==nCh,5),1)); %Alpha
    end
    subplot(1,2,nD); format_fig; %On nD==2 this line deletes prev graph & starts new figure
    %     simpleTopoPlot_ft(temp_topo(correspCh2),layout,'on',[],0,1); %DP - this was the original, but seems to plot incorrect layout
    simpleTopoPlot_ft(temp_topo(correspCh),layout,'on',[],0,1); %DP - looks better but double check
    colormap(cmap);
    %     if nD==1
    %         hb=colorbar('Position',[0.9195    0.6373    0.0143    0.2881]);
    %     end
    if nD==1
        title('Delta Freq Young')
        colorbar;
    elseif nD==2
        title('Delta Freq Old')
    end
    maxmin(nD,1)=min((temp_topo));
    maxmin(nD,2)=max((temp_topo));
end
caxis([min(maxmin(:,1)) max(maxmin(:,2))])

figure;
for nD=1:2
    temp_topo=[];
    for nCh=1:length(fractal.label)
        temp_topo(nCh)=squeeze(nanmean(all_deltapeak_component(all_deltapeak_component(:,2)==agegroups{nD} & all_deltapeak_component(:,4)==nCh,6),1)); %Alpha
    end
    subplot(1,2,nD); format_fig;
    %     simpleTopoPlot_ft(temp_topo(correspCh2),layout,'on',[],0,1); %DP - this was the original, but seems to plot incorrect layout
    simpleTopoPlot_ft(temp_topo(correspCh),layout,'on',[],0,1); %DP - looks better but double check
    colormap(cmap);
    %     if nD==1
    %         hb=colorbar('Position',[0.9195    0.6373    0.0143    0.2881]);
    %     end
    if nD==1
        title('Delta Peak Amplitude Young')
        colorbar;
    elseif nD==2
        title('Delta Peak Amplitude Old')
    end
    maxmin(nD,1)=min((temp_topo));
    maxmin(nD,2)=max((temp_topo));
end
caxis([min(maxmin(:,1)) max(maxmin(:,2))])

figure; oldtopo=[]; youngtopo=[];
for nCh=1:length(fractal.label)
    youngtopo(nCh)=squeeze(nanmean(all_deltapeak_component(all_deltapeak_component(:,2)==agegroups{1} & all_deltapeak_component(:,4)==nCh,5),1));
    oldtopo(nCh)=squeeze(nanmean(all_deltapeak_component(all_deltapeak_component(:,2)==agegroups{2} & all_deltapeak_component(:,4)==nCh,5),1));
end
difftopo=oldtopo-youngtopo;
simpleTopoPlot_ft(difftopo(correspCh),layout,'on',[],0,1); %DP - looks better but double check
colormap(cmap);
colorbar;
title('Old-Young Freq');

figure; oldtopo=[]; youngtopo=[];
for nCh=1:length(fractal.label)
    youngtopo(nCh)=squeeze(nanmean(all_deltapeak_component(all_deltapeak_component(:,2)==agegroups{1} & all_deltapeak_component(:,4)==nCh,6),1));
    oldtopo(nCh)=squeeze(nanmean(all_deltapeak_component(all_deltapeak_component(:,2)==agegroups{2} & all_deltapeak_component(:,4)==nCh,6),1));
end
difftopo=oldtopo-youngtopo;
simpleTopoPlot_ft(difftopo(correspCh),layout,'on',[],0,1); %DP - looks better but double check
colormap(cmap);
colorbar;
title('Old-Young Amp');

%% Theta peaks (4-7Hz)

all_thetapeak_component=all_peak_component(all_peak_component(:,5)>4 & all_peak_component(:,5)<7,:);

figure;
agegroups={[4],[5]}; % Young vs Old
for nD=1:2
    temp_topo=[];
    for nCh=1:length(fractal.label)
        temp_topo(nCh)=squeeze(nanmean(all_thetapeak_component(all_thetapeak_component(:,2)==agegroups{nD} & all_thetapeak_component(:,4)==nCh,5),1)); %Alpha
    end
    subplot(1,2,nD); format_fig; %On nD==2 this line deletes prev graph & starts new figure
    %     simpleTopoPlot_ft(temp_topo(correspCh2),layout,'on',[],0,1); %DP - this was the original, but seems to plot incorrect layout
    simpleTopoPlot_ft(temp_topo(correspCh),layout,'on',[],0,1); %DP - looks better but double check
    colormap(cmap);
    %     if nD==1
    %         hb=colorbar('Position',[0.9195    0.6373    0.0143    0.2881]);
    %     end
    if nD==1
        title('Theta Freq Young')
        colorbar;
    elseif nD==2
        title('Theta Freq Old')
    end
    maxmin(nD,1)=min((temp_topo));
    maxmin(nD,2)=max((temp_topo));
end
caxis([min(maxmin(:,1)) max(maxmin(:,2))])

figure;
for nD=1:2
    temp_topo=[];
    for nCh=1:length(fractal.label)
        temp_topo(nCh)=squeeze(nanmean(all_thetapeak_component(all_thetapeak_component(:,2)==agegroups{nD} & all_thetapeak_component(:,4)==nCh,6),1)); %Alpha
    end
    subplot(1,2,nD); format_fig;
    %     simpleTopoPlot_ft(temp_topo(correspCh2),layout,'on',[],0,1); %DP - this was the original, but seems to plot incorrect layout
    simpleTopoPlot_ft(temp_topo(correspCh),layout,'on',[],0,1); %DP - looks better but double check
    colormap(cmap);
    %     if nD==1
    %         hb=colorbar('Position',[0.9195    0.6373    0.0143    0.2881]);
    %     end
    if nD==1
        title('Theta Peak Amplitude Young')
        colorbar;
    elseif nD==2
        title('Theta Peak Amplitude Old')
    end
    maxmin(nD,1)=min((temp_topo));
    maxmin(nD,2)=max((temp_topo));
end
caxis([min(maxmin(:,1)) max(maxmin(:,2))])

figure; oldtopo=[]; youngtopo=[];
for nCh=1:length(fractal.label)
    youngtopo(nCh)=squeeze(nanmean(all_thetapeak_component(all_thetapeak_component(:,2)==agegroups{1} & all_thetapeak_component(:,4)==nCh,5),1));
    oldtopo(nCh)=squeeze(nanmean(all_thetapeak_component(all_thetapeak_component(:,2)==agegroups{2} & all_thetapeak_component(:,4)==nCh,5),1));
end
difftopo=oldtopo-youngtopo;
simpleTopoPlot_ft(difftopo(correspCh),layout,'on',[],0,1); %DP - looks better but double check
colormap(cmap);
colorbar;
title('Old-Young Freq');

figure; oldtopo=[]; youngtopo=[];
for nCh=1:length(fractal.label)
    youngtopo(nCh)=squeeze(nanmean(all_thetapeak_component(all_thetapeak_component(:,2)==agegroups{1} & all_thetapeak_component(:,4)==nCh,6),1));
    oldtopo(nCh)=squeeze(nanmean(all_thetapeak_component(all_thetapeak_component(:,2)==agegroups{2} & all_thetapeak_component(:,4)==nCh,6),1));
end
difftopo=oldtopo-youngtopo;
simpleTopoPlot_ft(difftopo(correspCh),layout,'on',[],0,1); %DP - looks better but double check
colormap(cmap);
colorbar;
title('Old-Young Amp');


%% Alpha peaks (8-11Hz)

all_alphapeak_component=all_peak_component(all_peak_component(:,5)>8 & all_peak_component(:,5)<11,:);

figure;
agegroups={[4],[5]}; % Young vs Old
for nD=1:2
    temp_topo=[];
    for nCh=1:length(fractal.label)
        temp_topo(nCh)=squeeze(nanmean(all_alphapeak_component(all_alphapeak_component(:,2)==agegroups{nD} & all_alphapeak_component(:,4)==nCh,5),1)); %Alpha
    end
    subplot(1,2,nD); format_fig; %On nD==2 this line deletes prev graph & starts new figure
    %     simpleTopoPlot_ft(temp_topo(correspCh2),layout,'on',[],0,1); %DP - this was the original, but seems to plot incorrect layout
    simpleTopoPlot_ft(temp_topo(correspCh),layout,'on',[],0,1); %DP - looks better but double check
    colormap(cmap);
    %     if nD==1
    %         hb=colorbar('Position',[0.9195    0.6373    0.0143    0.2881]);
    %     end
    if nD==1
        title('Alpha Freq Young')
        colorbar;
    elseif nD==2
        title('Alpha Freq Old')
    end
    maxmin(nD,1)=min((temp_topo));
    maxmin(nD,2)=max((temp_topo));
end
caxis([min(maxmin(:,1)) max(maxmin(:,2))])

figure;
for nD=1:2
    temp_topo=[];
    for nCh=1:length(fractal.label)
        temp_topo(nCh)=squeeze(nanmean(all_alphapeak_component(all_alphapeak_component(:,2)==agegroups{nD} & all_alphapeak_component(:,4)==nCh,6),1)); %Alpha
    end
    subplot(1,2,nD); format_fig;
    %     simpleTopoPlot_ft(temp_topo(correspCh2),layout,'on',[],0,1); %DP - this was the original, but seems to plot incorrect layout
    simpleTopoPlot_ft(temp_topo(correspCh),layout,'on',[],0,1); %DP - looks better but double check
    colormap(cmap);
    %     if nD==1
    %         hb=colorbar('Position',[0.9195    0.6373    0.0143    0.2881]);
    %     end
    
    if nD==1
        title('Alpha Peak Amplitude Young')
        colorbar;
    elseif nD==2
        title('Alpha Peak Amplitude Old')
    end
    maxmin(nD,1)=min((temp_topo));
    maxmin(nD,2)=max((temp_topo));
end
caxis([min(maxmin(:,1)) max(maxmin(:,2))])

figure; oldtopo=[]; youngtopo=[];
for nCh=1:length(fractal.label)
    youngtopo(nCh)=squeeze(nanmean(all_alphapeak_component(all_alphapeak_component(:,2)==agegroups{1} & all_alphapeak_component(:,4)==nCh,5),1));
    oldtopo(nCh)=squeeze(nanmean(all_alphapeak_component(all_alphapeak_component(:,2)==agegroups{2} & all_alphapeak_component(:,4)==nCh,5),1));
end
difftopo=oldtopo-youngtopo;
simpleTopoPlot_ft(difftopo(correspCh),layout,'on',[],0,1); %DP - looks better but double check
colormap(cmap);
colorbar;
title('Old-Young Freq');

figure; oldtopo=[]; youngtopo=[];
for nCh=1:length(fractal.label)
    youngtopo(nCh)=squeeze(nanmean(all_alphapeak_component(all_alphapeak_component(:,2)==agegroups{1} & all_alphapeak_component(:,4)==nCh,6),1));
    oldtopo(nCh)=squeeze(nanmean(all_alphapeak_component(all_alphapeak_component(:,2)==agegroups{2} & all_alphapeak_component(:,4)==nCh,6),1));
end
difftopo=oldtopo-youngtopo;
simpleTopoPlot_ft(difftopo(correspCh),layout,'on',[],0,1); %DP - looks better but double check
colormap(cmap);
colorbar;
title('Old-Young Amp');


%% Beta peaks (12-29Hz)

all_betapeak_component=all_peak_component(all_peak_component(:,5)>12 & all_peak_component(:,5)<29,:);

figure;
agegroups={[4],[5]}; % Young vs Old
for nD=1:2
    temp_topo=[];
    for nCh=1:length(fractal.label)
        temp_topo(nCh)=squeeze(nanmean(all_betapeak_component(all_betapeak_component(:,2)==agegroups{nD} & all_betapeak_component(:,4)==nCh,5),1)); %Alpha
    end
    subplot(1,2,nD); format_fig; %On nD==2 this line deletes prev graph & starts new figure
    %     simpleTopoPlot_ft(temp_topo(correspCh2),layout,'on',[],0,1); %DP - this was the original, but seems to plot incorrect layout
    simpleTopoPlot_ft(temp_topo(correspCh),layout,'on',[],0,1); %DP - looks better but double check
    colormap(cmap);
    %     if nD==1
    %         hb=colorbar('Position',[0.9195    0.6373    0.0143    0.2881]);
    %     end
    if nD==1
        title('Beta Freq Young')
        colorbar;
    elseif nD==2
        title('Beta Freq Old')
        colorbar;
    end
    maxmin(nD,1)=min((temp_topo));
    maxmin(nD,2)=max((temp_topo));
end
caxis([min(maxmin(:,1)) max(maxmin(:,2))])

figure;
for nD=1:2
    temp_topo=[];
    for nCh=1:length(fractal.label)
        temp_topo(nCh)=squeeze(nanmean(all_betapeak_component(all_betapeak_component(:,2)==agegroups{nD} & all_betapeak_component(:,4)==nCh,6),1)); %Alpha
    end
    subplot(1,2,nD); format_fig;
    %     simpleTopoPlot_ft(temp_topo(correspCh2),layout,'on',[],0,1); %DP - this was the original, but seems to plot incorrect layout
    simpleTopoPlot_ft(temp_topo(correspCh),layout,'on',[],0,1); %DP - looks better but double check
    colormap(cmap);
    %     if nD==1
    %         hb=colorbar('Position',[0.9195    0.6373    0.0143    0.2881]);
    %     end
    if nD==1
        title('Beta Peak Amplitude Young')
        %         colorbar;
    elseif nD==2
        title('Beta Peak Amplitude Old')
        %         colorbar;
    end
    maxmin(nD,1)=min((temp_topo));
    maxmin(nD,2)=max((temp_topo));
end
caxis([min(maxmin(:,1)) max(maxmin(:,2))])
colorbar;

figure; oldtopo=[]; youngtopo=[];
for nCh=1:length(fractal.label)
    youngtopo(nCh)=squeeze(nanmean(all_betapeak_component(all_betapeak_component(:,2)==agegroups{1} & all_betapeak_component(:,4)==nCh,5),1));
    oldtopo(nCh)=squeeze(nanmean(all_betapeak_component(all_betapeak_component(:,2)==agegroups{2} & all_betapeak_component(:,4)==nCh,5),1));
end
difftopo=oldtopo-youngtopo;
simpleTopoPlot_ft(difftopo(correspCh),layout,'on',[],0,1); %DP - looks better but double check
colormap(cmap);
colorbar;
title('Old-Young Freq');

figure; oldtopo=[]; youngtopo=[];
for nCh=1:length(fractal.label)
    youngtopo(nCh)=squeeze(nanmean(all_betapeak_component(all_betapeak_component(:,2)==agegroups{1} & all_betapeak_component(:,4)==nCh,6),1));
    oldtopo(nCh)=squeeze(nanmean(all_betapeak_component(all_betapeak_component(:,2)==agegroups{2} & all_betapeak_component(:,4)==nCh,6),1));
end
difftopo=oldtopo-youngtopo;
simpleTopoPlot_ft(difftopo(correspCh),layout,'on',[],0,1); %DP - looks better but double check
colormap(cmap);
colorbar;
title('Old-Young Amp');

%% Proportion of Peaks
%all_detectedpeaks: participant #, group #, agegroup #, electrode, alpha peak, slow peak
young_alphapeaks=0; old_alphapeaks=0; young_slowpeaks=0; old_slowpeaks=0;
for pps=1:length(all_agegroup)
    if sum(all_detectedpeaks(all_detectedpeaks(:,1)==pps,5))>0
        if all_agegroup(1,pps)==0
            young_alphapeaks=young_alphapeaks+1;
        elseif all_agegroup(1,pps)==1
            old_alphapeaks=old_alphapeaks+1;
        end
    else
    end
    if sum(all_detectedpeaks(all_detectedpeaks(:,1)==pps,6))>0
        if all_agegroup(1,pps)==0
            young_slowpeaks=young_slowpeaks+1;
        elseif all_agegroup(1,pps)==1
            old_slowpeaks=old_slowpeaks+1;
        end
    else
    end
end

prop_peaks=[];
prop_peaks(1,1)=young_alphapeaks/sum(all_agegroup(1,:)==0);
prop_peaks(1,2)=young_slowpeaks/sum(all_agegroup(1,:)==0);
prop_peaks(2,1)=old_alphapeaks/sum(all_agegroup(1,:)==1);
prop_peaks(2,2)=old_slowpeaks/sum(all_agegroup(1,:)==1);

figure; x=categorical({'Young','Old'});
bar(x,prop_peaks);
xlabel('Group')
    ylabel('Proportion of Participants with Peak')
    format_fig;
    legend({'Alpha','Slow'})
    title(['Younger Older Peak Proportion'])

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