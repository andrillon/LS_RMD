%%
clear all
close all

run ../LS_RMD_localdef.m

addpath(path_fieldtrip);
ft_defaults;
addpath(genpath(path_LSCPtools));

files=dir([preproc_path filesep 'ICAcleaned_eblock_ft_*.mat']);

%%
all_nSW_perE=[];
all_P2P_perE=[];
all_DnS_perE=[];
all_UpS_perE=[];

all_waves_byTr=[];
all_P2P_byTr=[];
all_DnS_byTr=[];
all_UpS_byTr=[];
            
all_nSW=[];
all_newlabels=[];
num_ch_blocks=[];
nFc=0;
group=[];
agegroup=[];
erplabels=[];
for nF=1:length(files)
    file_name = files(nF).name;
    folder_name = files(nF).folder;
    SubID=file_name(1:end-4);
    seps=findstr(SubID,'ft_');
    SubID=SubID(seps(1)+3:end);
    
    load([preproc_path filesep '..' filesep 'SWdetection' filesep 'allSW_' SubID]); %,'slow_Waves','paramSW')
    if max(all_Waves(:,2))<8
        warning('subj has less than 8 blocks so skipping')
        continue;
    else
        fprintf('subj has %g blocks, taking the first 8\n',max(all_Waves(:,2)))
        all_Waves(all_Waves(:,2)>8,:)=[];
    end
    
    load([folder_name filesep file_name]);
%         etrial_data =    load([folder_name filesep 'ICf_etrial_ft_f_etrial_ft_' SubID],'data');

    tic;
    fprintf('... working on %s (%g/%g)\n',SubID,nF,length(files))
    nFc=nFc+1;
    if length(SubID)==4 && SubID(1)=='A' % OLD (MONASH) - UP & DOWN 90%COH
        group(nFc)=2;
        agegroup(nFc)=1;
    elseif length(SubID)==7 % YOUNG (MONASH) - DOWN 50%COH
        group(nFc)=1;
        agegroup(nFc)=0;
    elseif length(SubID)==11 % YOUNG (TRINITY) - DOWN 50%COH
        group(nFc)=3;
        agegroup(nFc)=0;
    elseif length(SubID)==5 && SubID(3)=='8' % YOUNG - UP & DOWN 90%COH
        group(nFc)=4;
        agegroup(nFc)=0;
    elseif length(SubID)==5 && SubID(3)=='9'% OLD - UP & DOWN 90%COH
        group(nFc)=5;
        agegroup(nFc)=1;
    else
        group(nFc)=NaN;
        agegroup(nFc)=NaN;
    end
    
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
    
    
    %%% clean detection
    paramSW.prticle_Thr=90; % 80 or 90 or 95
    paramSW.AmpCriterionIdx=4; % 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
    paramSW.fixThr=[];
    paramSW.art_ampl=150;
    paramSW.max_posampl=75;
    paramSW.max_Freq=8;
    paramSW.min_Freq=4;
    
    all_Waves=double(all_Waves);
    all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./data.fsample);
    fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq>paramSW.max_Freq | all_freq<paramSW.min_Freq)*100)
    fprintf('... ... %g %% waves discarded because of max P2P ampl\n',mean(all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl)*100)
    fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(all_Waves(:,11)>paramSW.max_posampl | all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl)*100)
    all_Waves(all_freq>paramSW.max_Freq | all_freq<paramSW.min_Freq | all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl | all_Waves(:,11)>paramSW.max_posampl| all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl,:)=[];
    slow_Waves=[];
    for nE=1:size(data.label,1)
        thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
        temp_p2p=thisE_Waves(:,paramSW.AmpCriterionIdx);
        
        if ~isempty(paramSW.fixThr)
            thr_Wave(nFc,nE)=paramSW.fixThr;
        else
            thr_Wave(nFc,nE)=prctile(thisE_Waves(:,paramSW.AmpCriterionIdx),paramSW.prticle_Thr);
        end
    
        thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
        temp_p2p=thisE_Waves(:,paramSW.AmpCriterionIdx);
        
        slow_Waves=[slow_Waves ; thisE_Waves(temp_p2p>thr_Wave(nFc,nE),:)];
    end
    
    if exist([preproc_path filesep 'behav_' SubID '.csv'])==0
        continue;
    end
        table=readtable([preproc_path filesep 'behav_' SubID '.csv']);
    table.RT(table.RT>3.5)=nan;
    
    fprintf('... ... %g channels and %g blocks\n',length(newlabels),length(data.trial))
    num_ch_blocks=[num_ch_blocks ; [nF length(newlabels) length(data.trial) group(nFc) agegroup(nFc) mean(cellfun(@length,data.time)/data.fsample/60) NaN]];
    all_newlabels=[all_newlabels ; newlabels];
    
    nout=histc(slow_Waves(:,3),1:length(newlabels));
    all_nSW_perE=[all_nSW_perE ; [nF 0 group(nFc) agegroup(nFc) nout'./(sum(cellfun(@length,data.time(1:5)))/data.fsample/60) nanmean(table.RT) mean(isnan(table.RT))]];
    all_nSW(nFc,1)=mean(nout);
    
    P2P_out=[];
    DnS_out =[];
    UpS_out =[];
    for nE=1:length(newlabels)
        P2P_out(nE)=nanmean(slow_Waves(slow_Waves(:,3)==nE,4));
        DnS_out(nE)=nanmean(slow_Waves(slow_Waves(:,3)==nE,12));
        UpS_out(nE)=nanmean(slow_Waves( slow_Waves(:,3)==nE,13));
    end
    all_P2P_perE=[all_P2P_perE ; [nF 0 group(nFc) agegroup(nFc) P2P_out]];
    all_DnS_perE=[all_DnS_perE ; [nF 0 group(nFc) agegroup(nFc) DnS_out]];
    all_UpS_perE=[all_UpS_perE ; [nF 0 group(nFc) agegroup(nFc) UpS_out]];
    
        
    fprintf('B: %g - T: %2.0f\n',0,0)
    for nB=1:8 %length(data.time)
        temp_Waves=slow_Waves(slow_Waves(:,2)==nB,:);
        temp_Behav=table(table.Block==nB,:);
        
        for nTr=1:size(temp_Behav,1)
            fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\bB: %g - T: %2.0f\n',nB,nTr)
            start_Tr=temp_Behav.Stat(nTr)/2; % divided by two because of the resampling
            end_Tr=temp_Behav.End(nTr)/2;
            
            waves_Tr=temp_Waves(temp_Waves(:,5)>=start_Tr & temp_Waves(:,5)<=end_Tr,:);
            out_waves_Tr=histc(waves_Tr(:,3),1:length(newlabels));
            out_waves_Tr=out_waves_Tr./((end_Tr-start_Tr)/250);
            if size(out_waves_Tr,1)>1
                out_waves_Tr=out_waves_Tr';
            end
            all_waves_byTr=[all_waves_byTr ; [nB nTr group(nFc) agegroup(nFc) temp_Behav.Cond(nTr) temp_Behav.RT(nTr) out_waves_Tr]];
            
            
            temp_P2P_out = nan(1,length(newlabels));
            temp_DnS_out = nan(1,length(newlabels));
            temp_UpS_out = nan(1,length(newlabels));
            for nE=1:length(newlabels)
                temp_P2P_out(nE)=nanmean(waves_Tr(waves_Tr(:,3)==nE,4));
                temp_DnS_out(nE)=nanmean(waves_Tr(waves_Tr(:,3)==nE,12));
                temp_UpS_out(nE)=nanmean(waves_Tr( waves_Tr(:,3)==nE,13));
            end
            all_P2P_byTr=[all_P2P_byTr ; [nB nTr group(nFc) agegroup(nFc) temp_Behav.Cond(nTr) temp_Behav.RT(nTr) temp_P2P_out]];
            all_DnS_byTr=[all_DnS_byTr ; [nB nTr group(nFc) agegroup(nFc) temp_Behav.Cond(nTr) temp_Behav.RT(nTr) temp_P2P_out]];
            all_UpS_byTr=[all_UpS_byTr ; [nB nTr group(nFc) agegroup(nFc) temp_Behav.Cond(nTr) temp_Behav.RT(nTr) temp_P2P_out]];

        end
    end
end

%% Layout
av_thr_Wave2=mean(thr_Wave);

cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.channel=newlabels;
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);

figure;
temp_topo=[];
for nCh=1:length(layout.label)-2
    temp_topo(nCh)=av_thr_Wave2(match_str(newlabels,layout.label(nCh)));
end
temp_topo(match_str(layout.label,{'TP7','TP8'}))=NaN;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
title('Threshold')
colorbar;

figure;
agegroup_labels={'Young','Old'};
for nage=1:2
    subplot(2,1,nage)
    temp_topo=[];
    for nCh=1:length(layout.label)-2
        temp_topo(nCh)=nanmean(thr_Wave(agegroup==nage-1 & group~=2,match_str(newlabels,layout.label(nCh))),1);
    end
    temp_topo(match_str(layout.label,{'TP7','TP8'}))=NaN;
    simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
    title(['Threshold - ' agegroup_labels{nage}])
    colorbar;
caxis([15 35])
end
%% Plot topo of SW dens, amplitude and slope
offset=4;
clims=[9.1 11.8; 19 37; 250 450; 250 450];
figure;
subplot(2,2,1);
temp_topo=[];
for nCh=1:length(layout.label)-2
    temp_topo(nCh)=nanmean(all_nSW_perE(all_UpS_perE(:,3)~=2,offset+match_str(newlabels,layout.label(nCh))),1);
end
temp_topo(match_str(layout.label,{'TP7','TP8'}))=NaN;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
title('Density')
colorbar; caxis(clims(1,:))

subplot(2,2,2);
temp_topo=[];
for nCh=1:length(layout.label)-2
    temp_topo(nCh)=nanmean(all_P2P_perE(all_UpS_perE(:,3)~=2,offset+match_str(newlabels,layout.label(nCh))),1);
end
temp_topo(match_str(layout.label,{'TP7','TP8'}))=NaN;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
title('P2P')
colorbar; caxis(clims(2,:))

subplot(2,2,3);
temp_topo=[];
for nCh=1:length(layout.label)-2
    temp_topo(nCh)=nanmean(all_DnS_perE(all_UpS_perE(:,3)~=2,offset+match_str(newlabels,layout.label(nCh))),1);
end
temp_topo(match_str(layout.label,{'TP7','TP8'}))=NaN;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
title('Downward Slope')
colorbar; caxis(clims(3,:))

subplot(2,2,4);
temp_topo=[];
for nCh=1:length(layout.label)-2
    temp_topo(nCh)=nanmean(all_UpS_perE(all_UpS_perE(:,3)~=2,offset+match_str(newlabels,layout.label(nCh))),1);
end
temp_topo(match_str(layout.label,{'TP7','TP8'}))=NaN;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
title('Upward Slope')
colorbar; caxis(clims(4,:))

%%%% per group
figure;
for nage=1:2
    subplot(2,4,1+(nage-1)*4);
    temp_topo=[];
    for nCh=1:length(layout.label)-2
        temp_topo(nCh)=nanmean(all_nSW_perE(all_nSW_perE(:,4)==nage-1 & all_nSW_perE(:,3)~=2,offset+match_str(newlabels,layout.label(nCh))),1);
    end
    temp_topo(match_str(layout.label,{'TP7','TP8'}))=NaN;
    simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
    title(['Density' agegroup_labels{nage}])
colorbar; caxis(clims(1,:))
    
    subplot(2,4,2+(nage-1)*4);
    temp_topo=[];
    for nCh=1:length(layout.label)-2
        temp_topo(nCh)=nanmean(all_P2P_perE(all_P2P_perE(:,4)==nage-1 & all_nSW_perE(:,3)~=2,offset+match_str(newlabels,layout.label(nCh))),1);
    end
    temp_topo(match_str(layout.label,{'TP7','TP8'}))=NaN;
    simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
    title(['P2P' agegroup_labels{nage}])
colorbar; caxis(clims(2,:))

    subplot(2,4,3+(nage-1)*4);
    temp_topo=[];
    for nCh=1:length(layout.label)-2
        temp_topo(nCh)=nanmean(all_DnS_perE(all_DnS_perE(:,4)==nage-1 & all_nSW_perE(:,3)~=2,offset+match_str(newlabels,layout.label(nCh))),1);
    end
    temp_topo(match_str(layout.label,{'TP7','TP8'}))=NaN;
    simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
    title(['Downward Slope' agegroup_labels{nage}])
colorbar; caxis(clims(3,:))

    subplot(2,4,4+(nage-1)*4);
    temp_topo=[];
    for nCh=1:length(layout.label)-2
        temp_topo(nCh)=nanmean(all_UpS_perE(all_UpS_perE(:,4)==nage-1 & all_nSW_perE(:,3)~=2,offset+match_str(newlabels,layout.label(nCh))),1);
    end
    temp_topo(match_str(layout.label,{'TP7','TP8'}))=NaN;
    simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
    title(['Upward Slope' agegroup_labels{nage}])
colorbar; caxis(clims(4,:))
end

%% stats
%%%%% per group
cmap=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap=flipud(cmap);
zvalim=4;
figure;
subplot(1,4,1);
temp_topo=[];
temp_pV=[];
for nCh=1:length(layout.label)-2
    [pV,~,stats]=ranksum(all_nSW_perE(all_nSW_perE(:,4)==0 & all_nSW_perE(:,3)~=2,offset+match_str(newlabels,layout.label(nCh))),...
        all_nSW_perE(all_nSW_perE(:,4)==1 & all_nSW_perE(:,3)~=2,offset+match_str(newlabels,layout.label(nCh))));
    temp_topo(nCh)=stats.zval;
    temp_pV(nCh)=pV;
end
temp_topo(temp_pV>0.05)=0;
temp_topo(match_str(layout.label,{'TP7','TP8'}))=NaN;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
title(['Density' agegroup_labels{nage}])
colormap(cmap);
colorbar; caxis([-1 1]*zvalim)

subplot(1,4,2);
temp_topo=[];
for nCh=1:length(layout.label)-2
    [pV,~,stats]=ranksum(all_P2P_perE(all_P2P_perE(:,4)==0 & all_P2P_perE(:,3)~=2,offset+match_str(newlabels,layout.label(nCh))),...
        all_P2P_perE(all_P2P_perE(:,4)==1 & all_P2P_perE(:,3)~=2,offset+match_str(newlabels,layout.label(nCh))));
    temp_topo(nCh)=stats.zval;
    temp_pV(nCh)=pV;
end
temp_topo(temp_pV>0.05)=0;temp_topo(match_str(layout.label,{'TP7','TP8'}))=NaN;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
title(['P2P' agegroup_labels{nage}])
colormap(cmap);
colorbar; caxis([-1 1]*zvalim)

subplot(1,4,3);
temp_topo=[];
for nCh=1:length(layout.label)-2
    [pV,~,stats]=ranksum(all_DnS_perE(all_DnS_perE(:,4)==0 & all_DnS_perE(:,3)~=2,offset+match_str(newlabels,layout.label(nCh))),...
        all_DnS_perE(all_DnS_perE(:,4)==1 & all_DnS_perE(:,3)~=2,offset+match_str(newlabels,layout.label(nCh))));
    temp_topo(nCh)=stats.zval;
    temp_pV(nCh)=pV;
end
temp_topo(temp_pV>0.05)=0;temp_topo(match_str(layout.label,{'TP7','TP8'}))=NaN;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
title(['Downward Slope' agegroup_labels{nage}])
colormap(cmap);
colorbar; caxis([-1 1]*zvalim)

subplot(1,4,4);
temp_topo=[];
for nCh=1:length(layout.label)-2
    [pV,~,stats]=ranksum(all_UpS_perE(all_UpS_perE(:,4)==0 & all_UpS_perE(:,3)~=2,offset+match_str(newlabels,layout.label(nCh))),...
        all_UpS_perE(all_UpS_perE(:,4)==1 & all_UpS_perE(:,3)~=2,offset+match_str(newlabels,layout.label(nCh))));
    temp_topo(nCh)=stats.zval;
    temp_pV(nCh)=pV;
end
temp_topo(temp_pV>0.05)=0;temp_topo(match_str(layout.label,{'TP7','TP8'}))=NaN;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
title(['Upward Slope' agegroup_labels{nage}])
colormap(cmap);
colorbar; caxis([-1 1]*zvalim)

%% Plot by-block temporal dynamic of SW dens, amplitude and slope
temp_plot=[];
for nB=1:8
    temp_plot(:,nB)=(all_nSW_perE_perB(all_nSW_perE_perB(:,2)==nB,offset+match_str(newlabels,'Fz')));
end
figure;
subplot(1,3,1);
simpleTplot(1:8,temp_plot,0,'k',0,'-',0.5,1,0,1,2);
xlabel('Block Number')
ylabel('SW density')
format_fig;

subplot(1,3,2);
hp=[];
[~,hp(1)]=simpleTplot(1:8,temp_plot(agegroup==0,:),0,[0 0 1],0,'-',0.5,1,0,1,2);
[~,hp(2)]=simpleTplot(1:8,temp_plot(agegroup==1,:),0,[1 0 0],0,'-',0.5,1,0,1,2);
xlabel('Block Number')
ylabel('SW density')
legend(hp,{'Young','Old'})
format_fig;

subplot(1,3,3);
hp=[];
[~,hp(1)]=simpleTplot(1:8,temp_plot(group==1 | group==3,:),0,[0 1 0],0,'-',0.5,1,0,1,2);
[~,hp(2)]=simpleTplot(1:8,temp_plot(group==4,:),0,[0 0 1],0,'-',0.5,1,0,1,2);
[~,hp(3)]=simpleTplot(1:8,temp_plot(agegroup==1,:),0,[1 0 0],0,'-',0.5,1,0,1,2);
xlabel('Block Number')
ylabel('SW density')
legend(hp,{'Y-50COH','Y-90COH','O-90COH'})
format_fig;

% OLD (MONASH) - UP & DOWN 90%COH | group(nFc)=2; agegroup(nFc)=1;
% YOUNG (MONASH) - DOWN 50%COH | group(nFc)=1;agegroup(nFc)=0;
% YOUNG (TRINITY) - DOWN 50%COH | group(nFc)=3;agegroup(nFc)=0;
% YOUNG - UP & DOWN 90%COH | group(nFc)=4;agegroup(nFc)=0;
% OLD - UP & DOWN 90%COH | group(nFc)=5; agegroup(nFc)=1;

    
%%
figure
mylabels={'Fz','Cz','Pz','Oz'};
for nplot=1:length(mylabels)
    subplot(1,length(mylabels),nplot);
    for ngroup=1:2
        plot((-1:1/data.fsample:1),mean(all_ERP_perB(all_ERP_perB(:,5)==ngroup-1 & all_ERP_perB(:,3)==match_str(newlabels,mylabels{nplot}),7:end)))
        hold on;
    end
    xlim([-0.5 1])
    if nplot==1
        legend({'Young','Old'})
    end
    format_fig;
    xlabel('Time from onset (s)')
    ylabel('Voltage (\muV)')
end

%% Correlation with behaviour
[r_RT_perB, pV_RT_perB]=corr(all_nSW_perE_perB(:,5:62),all_nSW_perE_perB(:,63),'rows','pairwise','type','spearman');
[r_RT, pV_RT]=corr(all_nSW_perE(:,5:62),all_nSW_perE(:,63),'rows','pairwise','type','spearman');


figure;
temp_topo=r_RT_perB; temp_topo(match_str(layout.label,{'TP7','TP8'}))=NaN;
simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
title(['Correlation SW * RT'])
% colormap(cmap);
colorbar; caxis([-1 1]*max(abs(temp_topo)))


[r_Miss_perB, pV_Miss_perB]=corr(all_nSW_perE_perB(:,5:62),all_nSW_perE_perB(:,64),'rows','pairwise','type','spearman');
[r_Miss, pV_Miss]=corr(all_nSW_perE(:,5:62),all_nSW_perE(:,64),'rows','pairwise','type','spearman');


figure;
temp_topo=r_Miss_perB; temp_topo(match_str(layout.label,{'TP7','TP8'}))=NaN;
simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
title(['Correlation SW * Misses'])
% colormap(cmap);
colorbar; caxis([-1 1]*max(abs(temp_topo)))

%% Ploting by condition
mySubs=unique(all_waves_byTr(:,1));
offset=6;
clims=[9.1 11.8; 19 37; 250 450; 250 450];
figure;
subplot(1,4,1);
temp_topo=[];
for nCh=1:length(layout.label)-2
    for nS=1:length(mySubs)
        temp_topo(nCh,nS)=nanmean(all_waves_byTr(all_waves_byTr(:,1)==mySubs(nS),offset+match_str(newlabels,layout.label(nCh))),1)*60;
    end
end
temp_topo=mean(temp_topo,2);
temp_topo(match_str(layout.label,{'TP7','TP8'}))=NaN;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
title('Density - Cond 1&2')
colorbar; 

subplot(1,4,2);
temp_topo1=[];
for nCh=1:length(layout.label)-2
    for nS=1:length(mySubs)
        temp_topo1(nCh,nS)=nanmean(all_waves_byTr(all_waves_byTr(:,1)==mySubs(nS) & all_waves_byTr(:,5)==1,offset+match_str(newlabels,layout.label(nCh))),1)*60;
    end
end
temp_topo=mean(temp_topo1,2);
temp_topo(match_str(layout.label,{'TP7','TP8'}))=NaN;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
title('Density - Cond 1')
colorbar; 

subplot(1,4,3);
temp_topo2=[];
for nCh=1:length(layout.label)-2
    for nS=1:length(mySubs)
        temp_topo2(nCh,nS)=nanmean(all_waves_byTr(all_waves_byTr(:,1)==mySubs(nS) & all_waves_byTr(:,5)==2,offset+match_str(newlabels,layout.label(nCh))),1)*60;
    end
end
temp_topo=mean(temp_topo2,2);
temp_topo(match_str(layout.label,{'TP7','TP8'}))=NaN;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
title('Density - Cond 2')
colorbar; 


subplot(1,4,4);
[h, pV, ~, stats]=ttest(temp_topo1',temp_topo2');
temp_topo=stats.tstat;
temp_topo(match_str(layout.label,{'TP7','TP8'}))=NaN;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
title('t-valie - Cond 1 vs 2')
colorbar; 