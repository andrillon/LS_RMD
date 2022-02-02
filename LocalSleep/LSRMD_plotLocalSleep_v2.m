%%
clear all
close all

run ../LS_RMD_localdef.m

addpath(path_fieldtrip);
ft_defaults;
addpath(genpath(path_LSCPtools));

files=dir([preproc_path filesep 'ICAcleaned_eblock_ft_*.mat']);

%%
nFc=0;
for nF=1:length(files)
    file_name = files(nF).name;
    folder_name = files(nF).folder;
    SubID=file_name(1:end-4);
    seps=findstr(SubID,'ft_');
    SubID=SubID(seps(1)+3:end);
    load([folder_name filesep file_name]);
    
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
    
    load([preproc_path filesep '..' filesep 'SWdetection' filesep 'allSW_' SubID]); %,'slow_Waves','paramSW')
    
    %%% clean detection
    paramSW.prticle_Thr=90; % 80 or 90 or 95
    paramSW.LimFrqW=[1 4]; % [1 4] or [4 10]
    paramSW.AmpCriterionIdx=4; % 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
    paramSW.fixThr=[];
    paramSW.art_ampl=150;
    paramSW.max_posampl=75;
    paramSW.max_Freq=7;
    
    all_Waves=double(all_Waves);
    all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./data.fsample);
    fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq>paramSW.max_Freq)*100)
    fprintf('... ... %g %% waves discarded because of max P2P ampl\n',mean(all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl)*100)
    fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(all_Waves(:,11)>paramSW.max_posampl | all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl)*100)
    all_Waves(all_freq>paramSW.max_Freq | all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl | all_Waves(:,11)>paramSW.max_posampl| all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl,:)=[];
    
    for nE=1:size(data.label,1)
        thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
        temp_p2p=thisE_Waves(:,paramSW.AmpCriterionIdx);
        
        if ~isempty(paramSW.fixThr)
            thr_Wave(nFc,nE)=paramSW.fixThr;
        else
            thr_Wave(nFc,nE)=prctile(thisE_Waves(:,paramSW.AmpCriterionIdx),paramSW.prticle_Thr);
        end
    end
end
av_thr_Wave=mean(thr_Wave(group==1,:));
%%
all_nSW_perE=[];
all_P2P_perE=[];
all_DnS_perE=[];
all_UpS_perE=[];
all_ERP_perB=[];

all_nSW=[];
all_newlabels=[];
num_ch_blocks=[];
nFc=0;
group=[];
agegroup=[];
for nF=1:length(files)
    file_name = files(nF).name;
    folder_name = files(nF).folder;
    SubID=file_name(1:end-4);
    seps=findstr(SubID,'ft_');
    SubID=SubID(seps(1)+3:end);
    load([folder_name filesep file_name]);
    
    etrial_data =    load([folder_name filesep 'ICf_etrial_ft_f_etrial_ft_' SubID],'data');
    
    if length(etrial_data.data.trial)/length(data.trial)<15 || length(data.trial)<9
        fprintf('... skipping %s (%g/%g) (less than 15 trials per block on average or less than 9 blocks)\n',SubID,nF,length(files))
        continue;
    end
    
    tic;
    fprintf('... working on %s (%g/%g)\n',SubID,nF,length(files))
    nFc=nFc+1;
    if length(SubID)==4 && SubID(1)=='A' % OLD - UP & DOWN 90%COH - MONASH
        group(nFc)=2;
        agegroup(nFc)=1;
    elseif length(SubID)==7 % YOUNG - DOWN 50%COH - MONASH
        group(nFc)=1;
        agegroup(nFc)=0;
    elseif length(SubID)==11 % YOUNG - DOWN 50%COH - TRINITY
        group(nFc)=3;
        agegroup(nFc)=0;
    elseif length(SubID)==5 && SubID(3)=='8' % YOUNG - UP & DOWN 90%COH - MONASH
        group(nFc)=4;
        agegroup(nFc)=0;
    elseif length(SubID)==5 && SubID(3)=='9'% OLD - UP & DOWN 90%COH - MONASH
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
    
    load([preproc_path filesep '..' filesep 'SWdetection' filesep 'allSW_' SubID]); %,'slow_Waves','paramSW')
    
    %%% clean detection
    paramSW.prticle_Thr=90; % 80 or 90 or 95
    paramSW.LimFrqW=[1 4]; % [1 4] or [4 10]
    paramSW.AmpCriterionIdx=4; % 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
    paramSW.fixThr=[];
    paramSW.art_ampl=150;
    paramSW.max_posampl=75;
    paramSW.max_Freq=7;
    
    all_Waves=double(all_Waves);
    all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./data.fsample);
    fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq>paramSW.max_Freq)*100)
    fprintf('... ... %g %% waves discarded because of max P2P ampl\n',mean(all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl)*100)
    fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(all_Waves(:,11)>paramSW.max_posampl | all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl)*100)
    all_Waves(all_freq>paramSW.max_Freq | all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl | all_Waves(:,11)>paramSW.max_posampl| all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl,:)=[];
    
    slow_Waves=[];
    for nE=1:size(data.label,1)
        thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
        temp_p2p=thisE_Waves(:,paramSW.AmpCriterionIdx);
        
        slow_Waves=[slow_Waves ; thisE_Waves(temp_p2p>av_thr_Wave(nE),:)];
    end
    
    fprintf('... ... %g channels and %g blocks\n',length(newlabels),length(data.trial))
    num_ch_blocks=[num_ch_blocks ; [nF length(newlabels) length(data.trial) group(nFc) agegroup(nFc) mean(cellfun(@length,data.time)/data.fsample/60) length(etrial_data.data.trial)]];
    all_newlabels=[all_newlabels ; newlabels];
    for nB=1:length(data.time)
        nout=histc(slow_Waves(slow_Waves(:,2)==nB,3),1:length(newlabels));
        all_nSW_perE=[all_nSW_perE ; [nF nB nout'./(length(data.time{nB})/data.fsample/60)]];
        all_nSW(nFc,nB)=mean(nout);
        
        P2P_out=[];
        DnS_out =[];
        UpS_out =[];
        for nE=1:length(newlabels)
            P2P_out(nE)=nanmean(slow_Waves(slow_Waves(:,2)==nB & slow_Waves(:,3)==nE,4));
            DnS_out(nE)=nanmean(slow_Waves(slow_Waves(:,2)==nB & slow_Waves(:,3)==nE,12));
            UpS_out(nE)=nanmean(slow_Waves(slow_Waves(:,2)==nB & slow_Waves(:,3)==nE,13));
        end
        all_P2P_perE=[all_P2P_perE ; [nF nB P2P_out]];
        all_DnS_perE=[all_DnS_perE ; [nF nB DnS_out]];
        all_UpS_perE=[all_UpS_perE ; [nF nB UpS_out]];
        
        for nE=1:length(newlabels)
            these_SW=slow_Waves(slow_Waves(:,2)==nB & slow_Waves(:,3)==nE,5);
            temp_ERP=[];
            for nSW=1:length(these_SW)
                if these_SW(nSW)-data.fsample>0 && these_SW(nSW)+data.fsample<size(data.trial{nB},2)
                    tempEEG=data.trial{nB}(nE,these_SW(nSW)-data.fsample:these_SW(nSW)+data.fsample);
                    tempEEG=tempEEG-mean(tempEEG(1:data.fsample/2));
                    temp_ERP=[temp_ERP ; tempEEG];
                end
            end
            if size(temp_ERP,1)>6
                all_ERP_perB=[all_ERP_perB ; [nF nB nE group(nFc) agegroup(nFc) size(data.trial{nB},2)/data.fsample/60 mean(temp_ERP,1)]];
            end
        end
    end
end

%% Layout
av_thr_Wave2=mean(thr_Wave(group==5,:));

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
simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
title('Threshold')
colorbar;
%% Plot topo of SW dens, amplitude and slope
figure;
subplot(2,2,1);
temp_topo=[];
for nCh=1:length(layout.label)-2
    temp_topo(nCh)=nanmean(all_nSW_perE(:,2+match_str(newlabels,layout.label(nCh))),1);
end
temp_topo(match_str(layout.label,{'TP7','TP8'}))=NaN;
simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
title('Density')

subplot(2,2,2);
temp_topo=[];
for nCh=1:length(layout.label)-2
    temp_topo(nCh)=nanmean(all_P2P_perE(:,2+match_str(newlabels,layout.label(nCh))),1);
end
temp_topo(match_str(layout.label,{'TP7','TP8'}))=NaN;
simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
title('P2P')

subplot(2,2,3);
temp_topo=[];
for nCh=1:length(layout.label)-2
    temp_topo(nCh)=nanmean(all_DnS_perE(:,2+match_str(newlabels,layout.label(nCh))),1);
end
temp_topo(match_str(layout.label,{'TP7','TP8'}))=NaN;
simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
title('Downward Slope')

subplot(2,2,4);
temp_topo=[];
for nCh=1:length(layout.label)-2
    temp_topo(nCh)=nanmean(all_UpS_perE(:,2+match_str(newlabels,layout.label(nCh))),1);
end
temp_topo(match_str(layout.label,{'TP7','TP8'}))=NaN;
simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
title('Upward Slope')

%% Plot by-block temporal dynamic of SW dens, amplitude and slope
temp_plot=[];
temp_plot2=[];
for nB=1:8
    temp_plot(nB)=nanmean(all_P2P_perE(all_P2P_perE(:,2)==nB,2+match_str(newlabels,'Cz')));
    temp_plot2(nB)=sem(all_P2P_perE(all_P2P_perE(:,2)==nB,2+match_str(newlabels,'Cz')));
end
figure;
plot(1:8,temp_plot);
hold on;
plot(1:8,temp_plot-temp_plot2);
plot(1:8,temp_plot+temp_plot2);


%%
figure;
temp_topo=[];
for ngroup=1:5
    subplot(1,5,ngroup);
    temp_topo{ngroup}=[];
    for nCh=1:length(layout.label)-2
        temp_topo{ngroup}(nCh)=nanmedian(all_nSW_perE(ismember(all_nSW_perE(:,1),num_ch_blocks(group==ngroup,1)) & all_nSW_perE(:,2)<10,2+match_str(newlabels,layout.label(nCh))),1);
    end
    temp_topo{ngroup}(match_str(layout.label,{'TP7','TP8'}))=NaN;
    if sum(isnan(temp_topo{ngroup}))==length(temp_topo{ngroup})
        continue;
    end
    simpleTopoPlot_ft(temp_topo{ngroup}', layout,'on',[],0,1);
    title('Density')
    colorbar;
    %     caxis([0 12])
end

%%
figure;
temp_topo=[];
for ngroup=1:2
    subplot(1,2,ngroup);
    temp_topo{ngroup}=[];
    for nCh=1:length(layout.label)-2
        temp_topo{ngroup}(nCh)=nanmean(all_nSW_perE(ismember(all_nSW_perE(:,1),num_ch_blocks(agegroup==ngroup-1,1)) & all_nSW_perE(:,2)<10,2+match_str(newlabels,layout.label(nCh))),1);
    end
    temp_topo{ngroup}(match_str(layout.label,{'TP7','TP8'}))=NaN;
    simpleTopoPlot_ft(temp_topo{ngroup}', layout,'on',[],0,1);
    title('Density')
    colorbar;
    %     caxis([0 12])
end

%%
figure
mylabels={'Fz','Cz','Pz','Oz'};
for nplot=1:length(mylabels)
    subplot(1,length(mylabels),nplot);
    for ngroup=1:5
        plot((-1:1/data.fsample:1),mean(all_ERP_perB(all_ERP_perB(:,4)==ngroup & all_ERP_perB(:,3)==match_str(newlabels,mylabels{nplot}),7:end)))
        hold on;
    end
    xlim([-0.5 1])
    if nplot==1
        legend({'1','2','3','4','5'})
    end
end