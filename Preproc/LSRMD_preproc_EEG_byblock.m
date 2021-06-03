%%
clear all;
close all;

%%
run ../LS_RMD_localdef.m
addpath((path_fieldtrip));
ft_defaults;
addpath(genpath(path_LSCPtools));

folders=dir([data_path filesep '*_*']);

load(['..' filesep 'LS_RMD_Bad Trials_Channels.mat']);
load(['..' filesep 'LS_RMD_Bad_Components.mat']);


%% loop on subjects
redo=0;
for nF=1:length(folders)
    files=dir([folders(nF).folder filesep folders(nF).name filesep '*.eeg']);
    type_File=1;
    these_names={files.name};
    files(find(~(cellfun(@isempty,regexp(these_names,'RS')))))=[];
    SubID=folders(nF).name;
    if isempty(files)
        files=dir([folders(nF).folder filesep folders(nF).name filesep '*.bdf']);
        type_File=2;
    end
    
    if strcmp(SubID,'AA_15_04_14')
        warning('Skipping AA_15_04_14... missing event information');
        continue;
    elseif strcmp(SubID,'AR_08_04_14') || strcmp(SubID,'LK_07_04_14') || strcmp(SubID,'MH_14_04_14') || strcmp(SubID,'OM_07_05_14') ...
            || strcmp(SubID,'RM_06_05_14')
        warning(sprintf('Skipping %s... unsupported headers',SubID));
        continue;
    end
    %     if length(files)~=numBlocks
    %         warning(sprintf('Skipping %s... wrong number of blocks (%g vs %g)',SubID,length(files),numBlocks));
    %         continue;
    %     end
    if exist([preproc_path filesep 'ICf_etrial_ft_f_etrial_ft_' SubID '.mat'])==0
        warning(sprintf('... missing %s ICA file\n',SubID));
        continue;
    end
    
    numBlocks=length(files);
    fprintf('Processing %s...',SubID);
    tic;
    if redo==1 || exist([preproc_path filesep 'ICAcleaned_eblock_ft_' SubID '.mat'])==0
        %%% Loop on individual block files
        all_channels=[];
        for k=1:numBlocks
            if type_File==1
                this_file=dir([folders(nF).folder filesep folders(nF).name filesep '*M' num2str(k) '.eeg']);
                %                 if isempty(this_file)
                %                     this_file=dir([folders(nF).folder filesep folders(nF).name filesep '*_' num2str(k) '.eeg']);
                %                 end
            elseif type_File==2
                this_file=dir([folders(nF).folder filesep folders(nF).name filesep '*_' num2str(k) '.bdf']);
            end
            if isempty(this_file)
                continue;
            end
            
            file_name = this_file(1).name;
            file_folder = this_file(1).folder;
            hdr=ft_read_header([file_folder filesep file_name]);
            if isempty(all_channels)
                all_channels=hdr.label(find((cellfun(@isempty,regexp(hdr.label,'EOG')) & ~cellfun(@isempty,regexp(hdr.chantype,'eeg')))));
            else
                all_channels=intersect(all_channels,hdr.label(find((cellfun(@isempty,regexp(hdr.label,'EOG')) & ~cellfun(@isempty,regexp(hdr.chantype,'eeg'))))));
            end
        end
        
        % Epoch by trial
        for k=1:numBlocks
            if type_File==1
                this_file=dir([folders(nF).folder filesep folders(nF).name filesep '*M' num2str(k) '.eeg']);
            elseif type_File==2
                this_file=dir([folders(nF).folder filesep folders(nF).name filesep '*_' num2str(k) '.bdf']);
            end
            if isempty(this_file)
                continue;
            end
            file_name = this_file(1).name;
            file_folder = this_file(1).folder;
            FileID=file_name(1:end-4);
            fprintf('... working on subject %s (%g/%g) - file %s (%g/%g)\n',SubID,nF,length(folders),FileID,k,length(files))
            
            %%% Read headers
            hdr=ft_read_header([file_folder filesep file_name]);
            
            %%% Define epochs
            cfg=[];
            cfg.trialfun            = 'LS_RMD_blockfun';
            cfg.SubID               = SubID;
            cfg.dataset             = [file_folder filesep file_name];
            cfg.trialdef.prestim    = 0;
            cfg.trialdef.poststim   = 0;
            cfg = ft_definetrial(cfg);
            
            cfg.channel        = all_channels;
            cfg.demean         = 'yes';
            cfg.lpfilter       = 'yes';        % enable high-pass filtering
            cfg.lpfilttype     = 'but';
            cfg.lpfiltord      = 4;
            cfg.lpfreq         = 35;
            cfg.hpfilter       = 'no';        % enable high-pass filtering
            %             cfg.hpfilttype     = 'but';
            %             cfg.hpfiltord      = 4;
            %             cfg.hpfreq         = 0.1;
            cfg.dftfilter      = 'yes';        % enable notch filtering to eliminate power line noise
            cfg.dftfreq        = [50 100 150]; % set up the frequencies for notch filtering
            
            cfg.reref      = 'yes';
            cfg.refchannel = 'all';
            
            dat                = ft_preprocessing(cfg); % read raw data
            
            cfgbs=[];
            cfgbs.resamplefs      = 250;
            cfgbs.detrend         = 'no';
            cfgbs.demean          = 'yes';
            dat2                  = ft_resampledata(cfgbs,dat); % read raw data
            if k==1
                data=dat2;
            else
                data.trial=[data.trial dat2.trial];
                data.time=[data.time dat2.time];
            end
        end
        
        
        %%% take out trial
        thisF=match_str(LSRMDBadTrialsChannels.SubID,SubID);
        if isempty(thisF)
            continue;
        end
        badChannels=[];
        
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
        
        
        %%% Reject bad component
        thisF=match_str(LSRMDBadComponents.SubID,SubID);
        if isempty(thisF)
            warning(sprintf('... missing %s in ICA detection\n',SubID));
            continue;
        end
        cfg = [];
        these_components=LSRMDBadComponents.EyeComponents{thisF};
        if ~isempty(these_components)
            eval(sprintf('cfg.component = [%s];',these_components)); % to be removed component(s)
        else
            cfg.component=[];
        end

        load([preproc_path filesep 'ICf_etrial_ft_f_etrial_ft_' SubID '.mat'],'comp');
        data = ft_rejectcomponent(cfg, comp, data);
        
        cfg=[];
        cfg.demean          = 'yes';
        cfg.baselinewindow  = [-1 0];
        data = ft_preprocessing(cfg,data);
        
        save([preproc_path filesep 'ICAcleaned_eblock_ft_' SubID],'data');
        
    end
    toc;
end


