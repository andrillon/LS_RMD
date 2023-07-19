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
redo=1; complete=0;

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
    
    if redo==1 || exist([preproc_path SubID '_CPP_perTrial_ICAcleaned.mat'])==0
        
        subj_pow=[]; nFc=nFc+1;
        
        
        load([folder_name filesep file_name]);
        behav_table=readtable([folder_name filesep 'behav_' file_name(findstr(file_name,'ft_')+3:end-4) '.csv']);
        
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
        
        % remove NaN values
        temp_data=data;
        cfg=[];
        cfg.trials=find(~isnan(behav_table.RT));
        data=ft_redefinetrial(cfg,data);
        cfg=[];
        cfg.trials=find(behav_table.Cond==1);
        data_left=ft_redefinetrial(cfg,temp_data);
        cfg=[];
        cfg.trials=find(behav_table.Cond==2);
        data_right=ft_redefinetrial(cfg,temp_data);
        
        % remove weird extra field
        if isfield(data,'elec')
            data=rmfield(data,'elec')
        end
        
        if isfield(data_left,'elec')
            data_left=rmfield(data_left,'elec')
        end
        
        if isfield(data_right,'elec')
            data_right=rmfield(data_right,'elec')
        end
        
        fprintf('%2.0f-%2.0f\n',0,0)
        
        % re-locking
        cfg=[];
        cfg.offset =-round(behav_table.RT(~isnan(behav_table.RT))*data.fsample);
        data_rlock=ft_redefinetrial(cfg,data);
        data_left_rlock=ft_redefinetrial(cfg,data_left);
        data_right_rlock=ft_redefinetrial(cfg,data_right);
        
        % crop trials
        cfg=[];
        cfg.toilim =[0 1.0];
        data=ft_redefinetrial(cfg,data);
        data_left=ft_redefinetrial(cfg,data_left);
        data_right=ft_redefinetrial(cfg,data_right);
        
        cfg=[];
        cfg.toilim =[-1 1];
        data_rlock=ft_redefinetrial(cfg,data_rlock);
        data_left_rlock=ft_redefinetrial(cfg,data_left_rlock);
        data_right_rlock=ft_redefinetrial(cfg,data_right_rlock);
        
        % estimate the ERP
        %         trial_side=behav_table.Cond(~isnan(behav_table.RT));
        tlck(nF) = ft_timelockanalysis([],data);
        tlck_left(nF) = ft_timelockanalysis([],data_left);
        tlck_right(nF) = ft_timelockanalysis([],data_right);
        %         for side=1:2
        %             tlck_side(nF,side) = ft_timelockanalysis([],data(data.trial(find(trial_side==side))));
        %         end
        rlck(nF) = ft_timelockanalysis([],data_rlock);
        rlck_left(nF) = ft_timelockanalysis([],data_left_rlock);
        rlck_right(nF) = ft_timelockanalysis([],data_right_rlock);
        %         for side=1:2
        %             rlck_side(nF,side) = ft_timelockanalysis([],data(data.trial(find(trial_side==side))));
        %         end
        
        % pull out Pz and Pz*side
        PzERP(nF,:)=tlck(nF).avg(strcmp(tlck(nF).label,'Pz'),:);
        PzERP_side(nF,:,1)=tlck_left(nF).avg(strcmp(tlck_left(nF).label,'Pz'),:);
        PzERP_side(nF,:,2)=tlck_right(nF).avg(strcmp(tlck_right(nF).label,'Pz'),:);
        PzERP_rlock(nF,:)=rlck(nF).avg(strcmp(rlck(nF).label,'Pz'),:);
        PzERP_rlock_side(nF,:,1)=rlck_left(nF).avg(strcmp(rlck_left(nF).label,'Pz'),:);
        PzERP_rlock_side(nF,:,2)=rlck_right(nF).avg(strcmp(rlck_right(nF).label,'Pz'),:);
        
        % groups
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
        all_agegroup(nFc)=thisagegroup;
    end
end

% for pps=1:nFc
%     PzERP(pps,:)=tlck(pps).avg(strcmp(tlck(pps).label,'Pz'),:);
%     PzERP_rlock(pps,:)=rlck(pps).avg(strcmp(tlck(pps).label,'Pz'),:);
%
%     for side=1:2
%         PzERP_side{side,pps}=tlck(pps).avg(strcmp(tlck(pps).label,'Pz'),behav_table.Cond==side);
%     end
% end

% young_CPP=nanmean(PzERP(all_agegroup==0,:));
% old_CPP=nanmean(PzERP(all_agegroup==1,:));

t=tlck(1).time;
tr=rlck(1).time;

figure;
[~,hp(1)]=simpleTplot(t,PzERP(all_agegroup==0,:),0,'r',[0 0.05 0.0001 1000],'-',0.1,1,0,1,1);
[~,hp(2)]=simpleTplot(t,PzERP(all_agegroup==1,:),0,'b',[0 0.05 0.0001 1000],'-',0.1,1,0,1,1);
legend(hp,'Younger','Older','Location','northeast'); xlim([0 1]); hp(1).LineWidth=3; hp(2).LineWidth=3; xticks([0 .25 .50 .75 1.00]);
xlabel('Post-target (s)'); set(gcf,'Color','white'); title('CPP (Pz)');

figure;
[~,hp(1)]=simpleTplot(tr,PzERP_rlock(all_agegroup==0,:),0,'r',[0 0.05 0.0001 1000],'-',0.1,1,0,1,1);
[~,hp(2)]=simpleTplot(tr,PzERP_rlock(all_agegroup==1,:),0,'b',[0 0.05 0.0001 1000],'-',0.1,1,0,1,1);
legend(hp,'Younger','Older','Location','northeast'); hp(1).LineWidth=3; hp(2).LineWidth=3; xlim([-.4 .1]);
xlabel('Time From Response (s)'); set(gcf,'Color','white'); title('CPP (Pz)');

%% Onset
fs=500;

t1 = .150; t2 = 1.000;
t1p = .400; t2p = 1.000;

[H_LPF,G_LPF]=butter(4,2*(2/fs));

thresh=3;
CPP_onsets=[];
sliding_windows = 0;
window_size = 100; % in ms
for s = 1:length(all_group)
    %     figure, hold on
    if sliding_windows==0
        for side=1:2
            CPP_temp = squeeze(PzERP_side(s,:,side)); % time
            CPP_temp = filtfilt(H_LPF,G_LPF,double(CPP_temp));
            [peak_amp,peak_ts] = max(CPP_temp(find(t>=t1p & t<=t2p)));
            peak_t = t(peak_ts+length(find(t<t1p)));
            baseline_std = std(CPP_temp(find(t<t1)));
            baseline_mean = mean(CPP_temp(find(t<t1)));
            %             CPP_onset = t(find(CPP_temp(find(t>=t1 & t<=t2))>baseline_mean+(thresh*baseline_std),1)+length(find(t<t1)));
            half_amp = baseline_mean+((peak_amp-baseline_mean)/2);
                CPP_onset(side) = t(find(CPP_temp(find(t>=t1 & t<=peak_t))>half_amp,1)+length(find(t<t1)));
        end
        CPP_temp = squeeze(PzERP(s,:)); % time
        CPP_temp = filtfilt(H_LPF,G_LPF,double(CPP_temp));
        [peak_amp,peak_ts] = max(CPP_temp(find(t>=t1p & t<=t2p)));
        peak_t = t(peak_ts+length(find(t<t1p)));
        baseline_std = std(CPP_temp(find(t<t1)));
        baseline_mean = mean(CPP_temp(find(t<t1)));
        %             CPP_onset = t(find(CPP_temp(find(t>=t1 & t<=t2))>baseline_mean+(thresh*baseline_std),1)+length(find(t<t1)));
        half_amp = baseline_mean+((peak_amp-baseline_mean)/2);
            CPP_onset(3) = t(find(CPP_temp(find(t>=t1 & t<=peak_t))>half_amp,1)+length(find(t<t1)));
        
            %         else
            %             CPP_temp = movmean(squeeze(mean(ERP_conip_side(s,ch_CPP,:,side),2)),ceil(window_size/2));
            %             CPP_temp_t = movmean(t,ceil(window_size*2));
            %             [peak_amp,peak_ts] = max(CPP_temp(find(CPP_temp_t>=t1p & CPP_temp_t<=t2p)));
            %             peak_t = CPP_temp_t(peak_ts+length(find(CPP_temp_t<t1p)));
            %             baseline_std = std(CPP_temp(find(CPP_temp_t<t1)));
            %             baseline_mean = mean(CPP_temp(find(CPP_temp_t<t1)));
            % %             CPP_onset = CPP_temp_t(find(CPP_temp(find(CPP_temp_t>=t1 & CPP_temp_t<=t2))>baseline_mean+(thresh*baseline_std),1)+length(find(CPP_temp_t<t1)));
            %             half_amp = baseline_mean+((peak_amp-baseline_mean)/2);
            %             CPP_onset = CPP_temp_t(find(CPP_temp(find(CPP_temp_t>=t1 & CPP_temp_t<=peak_t))>half_amp,1)+length(find(CPP_temp_t<t1)));
        end
        
        if isempty(CPP_onset), CPP_onset = NaN; end
        CPP_onsets(s,:) = CPP_onset;
    end
    
    %% Slope
    
    t1=-.150; t2=-.05;
    CPPr_slopes = [];
    for s = 1:size(PzERP_rlock,1)
        coef = polyfit(tr(find(tr>=t1 & tr<=t2)),squeeze(PzERP_rlock(s,find(tr>=t1 & tr<=t2))),1);
        CPPr_slopes(s) = coef(1);
    end
    
    %% Save
    
    save([preproc_path 'CPP_stats'],'CPP_onsets','CPPr_slopes');
    
    %% Group Comparisons
    load([preproc_path 'CPP_stats']);
    
    
    
    %% Pull out participant data
    
    for pps=1:length(all_group)
        TF_all(pps,:,:)=tlck(pps).avg;
    end
    
    %% CPP Topo
    cfg = [];
    cfg.layout = 'acticap-64ch-standard2.mat';
    cfg.channel=newlabels;
    cfg.center      = 'yes';
    layout=ft_prepare_layout(cfg);
    correspCh=[];
    for nCh=1:length(layout.label)-2
        correspCh(nCh)=match_str(newlabels,layout.label(nCh));
    end
    
    t1 = .375; t2 = .475;
    
    figure; zvalim=6;
    temp_topo=[]; temp_pV=[];
    for nCh=1:length(newlabels)
        temp_topo(nCh)=squeeze(nanmean(nanmean(TF_all(all_agegroup==0,nCh,find(t>=t1 & t<=t2)),1),3));
    end
    temp_topo=temp_topo(correspCh);
    simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
    title(['CPP - Young'])
    colormap(parula);
    colorbar;
    caxis([-1 1]*zvalim);
    
    
    t1 = .6; t2 = .7;
    
    figure; zvalim=6;
    temp_topo=[]; temp_pV=[];
    for nCh=1:length(newlabels)
        temp_topo(nCh)=squeeze(nanmean(nanmean(TF_all(all_agegroup==1,nCh,find(t>=t1 & t<=t2)),1),3));
    end
    temp_topo=temp_topo(correspCh);
    simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
    title(['CPP - Old'])
    colormap(parula);
    colorbar;
    caxis([-1 1]*zvalim);