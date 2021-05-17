%%
clear all;
close all;

%%
run ../LS_RMD_localdef.m
addpath((path_fieldtrip));
ft_defaults;

files=dir([preproc_path filesep 'f_etrial_ft_*.mat']);

load(['..' filesep 'LS_RMD_Bad Trials_Channels.mat']);

%% loop on subjects
redo=0;
for nF=1:length(files)
    file_name=files(nF).name;
    SubID=file_name(13:end-4);
    fprintf('Processing %s...',SubID);
    tic;
    if redo==1 || exist([preproc_path filesep 'ICf_etrial_ft_' SubID '.mat'])==0
        %%% load epoched data
        load([preproc_path filesep 'f_etrial_ft_' SubID]);
        
        %%% take out trial
        thisF=match_str(LSRMDBadTrialsChannels.SubID,SubID);
        if isempty(thisF)
            continue;
        end
%         tempChannels=LSRMDBadTrialsChannels.ExcludedTrials{thisF}; tempChannels(tempChannels==' ')=[];
%         eval(sprintf('badChannels={%s};',tempChannels));
%         badChannels(cellfun(@isempty,badChannels))=[];
        badChannels=[];
        eval(['badTrials=[' LSRMDBadTrialsChannels.ExcludedTrials{thisF} '];']);
        
        cfg=[];
        cfg.trials          = setdiff(1:length(data.trial),badTrials);
        data = ft_preprocessing(cfg, data);
        
        if ~isempty(badChannels)
            fprintf('... ... interpolating %g channels\n',length(badChannels))
            % find neighbours
            cfg=[];
            cfg.method        = 'triangulation';
            cfg.layout        = layout;
            cfg.feedback      = 'no';
            cfg.channel = layout.label;
            [neighbours] = ft_prepare_neighbours(cfg);
            
            % interpolate channels
            cfg=[];
            cfg.method         = 'weighted';
            cfg.badchannel     = badChannels;
            cfg.missingchannel = [];
            cfg.neighbours     = neighbours;
            cfg.trials         = 'all';
            cfg.layout         = layout;
            cfg.channel = layout.label;
            [data] = ft_channelrepair(cfg, data);
        end
        
        cfg=[];
        cfg.reref      = 'yes';
        cfg.refchannel = 'all';
        data = ft_preprocessing(cfg,data);
        
        %%% run ICA
        rankICA = rank(data.trial{1,1});
        cfg        = [];
        cfg.method = 'runica'; % this is the default and uses the implementation from EEGLAB
        cfg.numcomponent = rankICA;
        comp = ft_componentanalysis(cfg, data);
        save([preproc_path filesep 'ICf_etrial_ft_' file_name(1:end-4)],'data','comp','rankICA');
        
    end
    toc;
end