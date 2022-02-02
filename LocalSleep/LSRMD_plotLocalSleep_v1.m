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
        
all_nSW=[];
all_newlabels=[];
num_ch_blocks=[];
nFc=0;
for nF=1:length(files)
    file_name = files(nF).name;
    folder_name = files(nF).folder;
    SubID=file_name(1:end-4);
    seps=findstr(SubID,'ft_');
    SubID=SubID(seps(1)+3:end);
    
    tic;
    fprintf('... working on %s (%g/%g)\n',SubID,nF,length(files))
    nFc=nFc+1;
    if length(SubID)==4 && SubID(1)=='A'
        group(nFc)=2;
        agegroup(nFc)=1;
    elseif length(SubID)==7
        group(nFc)=1;
        agegroup(nFc)=0;
    elseif length(SubID)==11
        group(nFc)=3;
        agegroup(nFc)=1;
    elseif length(SubID)==5 && SubID(3)=='8'
        group(nFc)=4;
        agegroup(nFc)=0;
    elseif length(SubID)==5 && SubID(3)=='9'
        group(nFc)=5;
        agegroup(nFc)=1;
    else
        group(nFc)=NaN;
        agegroup(nFc)=NaN;
    end
    
    load([preproc_path filesep '..' filesep 'SWdetection' filesep 'SW_' SubID]); %,'slow_Waves','paramSW')
    
    load([folder_name filesep file_name]);
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
    fprintf('... ... %g channels and %g blocks\n',length(newlabels),length(data.trial))
    num_ch_blocks=[num_ch_blocks ; [nF length(newlabels) length(data.trial) group(nFc) agegroup(nFc)]];
    all_newlabels=[all_newlabels ; newlabels];
    for nB=1:16
        nout=histc(slow_Waves(slow_Waves(:,2)==nB,3),1:length(newlabels));
        all_nSW_perE=[all_nSW_perE ; [nF nB nout']];
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

    end
end

%% Layout
cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.channel=newlabels;
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);

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
