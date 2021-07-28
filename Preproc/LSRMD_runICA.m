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
redo=1;
for nF=1:length(files)
    file_name=files(nF).name;
    SubID=file_name(13:end-4);
    fprintf('Processing %s...\n',SubID);
    tic;
    if redo==1 || exist([preproc_path filesep 'ICf_etrial_ft_' SubID '.mat'])==0
        
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
        eval(['badChannels=[' LSRMDBadTrialsChannels.ExcludedChannels{thisF} '];']);
        if isempty(badChannels)
            continue;
        end
        %%% load epoched data
        load([preproc_path filesep 'f_etrial_ft_' SubID]);
        
        cfg=[];
        cfg.trials          = setdiff(1:length(data.trial),badTrials);
        data = ft_preprocessing(cfg, data);
        
        if ~isempty(badChannels)
            fprintf('... ... interpolating %g channels\n',length(badChannels))
            % retrieve layout
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
            data.label=newlabels;
            cfg = [];
            cfg.layout = 'biosemi64.lay';
            cfg.channel=newlabels';
            cfg.center      = 'yes';
            layout=ft_prepare_layout(cfg);

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
            cfg.badchannel     = newlabels(badChannels);
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