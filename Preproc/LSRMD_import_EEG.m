%%
clear all;
close all;

%%
run ../LS_RMD_localdef.m
addpath((path_fieldtrip));
ft_defaults;
addpath(genpath(path_LSCPtools));

folders=dir([data_path filesep]); % DP 19/07 removing _ to incorporate older data
folders(1:2) = []; %DP 19/07 removing hidden folders - probably more elegant solution needed here!

%% loop on subjects
redo=1;
for nF=1:length(folders)
%     redo=0; %DP 19/07 adding so I can specify only older adults
    files=dir([folders(nF).folder filesep folders(nF).name filesep '*M*.eeg']); % DP 19/07 adding M to distinguish from older adults
    type_File=1;
    these_names={files.name};
    files(find(~(cellfun(@isempty,regexp(these_names,'RS')))))=[];
    SubID=folders(nF).name;
%     if ~strcmp(SubID(1),'A')
%         continue; %trying to solve fixation break errors
%     end
%     if strcmp(SubID,'A107') || strcmp(SubID,'A108') || strcmp(SubID,'A109') || strcmp(SubID,'A111')
%         warning('SKIPPING (problem with events)')
%         continue; %trying to solve fixation break errors
%     end
    if isempty(files)
        files=dir([folders(nF).folder filesep folders(nF).name filesep '*.bdf']);
        type_File=2;
        if isempty(files) %DP 19/07 adding Bryce older adult data
            files=dir([folders(nF).folder filesep folders(nF).name filesep SubID '90*.eeg']);
            type_File=3;
            
%             redo=0; %DP 19/07 adding so I can specify only older adults
            if isempty(files) %DP 19/07 adding Megan older adult data
                files=dir([folders(nF).folder filesep folders(nF).name filesep 'HN9*_*.eeg']);
                type_File=4;
%                 redo=0; %DP 19/07 adding so I can specify only older adults
                if isempty(files) %DP 28/07 adding Megan younger adult data
                    files=dir([folders(nF).folder filesep folders(nF).name filesep 'HN8*_*.eeg']);
                    type_File=5;
%                     redo=0; %DP 28/07 adding so I can specify only older adults
                end
            end
        end
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
  
numBlocks=length(files);
fprintf('Processing %s...',SubID);
    tic;
    if redo==1 || exist([preproc_path filesep 'f_etrial_ft_' SubID '.mat'])==0
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
            elseif type_File==3
                this_file=dir([folders(nF).folder filesep folders(nF).name filesep '*90' num2str(k) '.eeg']); %DP 19/07 adding Bryce older adult data
            elseif type_File==4 || type_File==5
                this_file=dir([folders(nF).folder filesep folders(nF).name filesep '*_' num2str(k) '.eeg']); %DP 28/07 adding Megan older adult + younger adult data
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
        
        for k=1:numBlocks
            if strcmp(SubID,'HN880') && k==1 % DP 29/07 - skipping block 1 - no events
                continue
            end
            if type_File==1
                this_file=dir([folders(nF).folder filesep folders(nF).name filesep '*M' num2str(k) '.eeg']);
            elseif type_File==2
                this_file=dir([folders(nF).folder filesep folders(nF).name filesep '*_' num2str(k) '.bdf']);
            elseif type_File==3
                this_file=dir([folders(nF).folder filesep folders(nF).name filesep '*90' num2str(k) '.eeg']); %DP 19/07 adding Bryce older adult data
                behavfiles=dir([folders(nF).folder filesep '..' filesep 'Behav' filesep SubID '90'  num2str(k) '.mat']);
            elseif type_File==4 || type_File==5
                this_file=dir([folders(nF).folder filesep folders(nF).name filesep '*_' num2str(k) '.eeg']); %DP 28/07 adding Megan older adult + younger adult data
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
            
            %%% Load behaviour data if relevant
            if exist('behavfiles')~=0 && ~isempty(behavfiles)
               behav_data=load([behavfiles.folder filesep behavfiles.name]);
            else
                behav_data=[];
            end
            
            %%% Define epochs
            cfg=[];
            cfg.trialfun            = 'LS_RMD_trialfun_v2';
            cfg.SubID               = SubID;
            cfg.behav               = behav_data;
            cfg.dataset             = [file_folder filesep file_name];
            cfg.trialdef.prestim    = 0.5;
            cfg.trialdef.poststim   = 0.2;
            cfg.type_File           = type_File;
            try
            cfg = ft_definetrial(cfg);
            catch
                warning('problem with trial definition');
                continue;
            end
            trl=cfg.trl;
            
            cfg.channel        = all_channels;
            cfg.demean         = 'yes';
            cfg.lpfilter       = 'yes';        % enable high-pass filtering
            cfg.lpfilttype     = 'but';
            cfg.lpfiltord      = 4;
            cfg.lpfreq         = 40;
            cfg.hpfilter       = 'yes';        % enable high-pass filtering
            cfg.hpfilttype     = 'but';
            cfg.hpfiltord      = 4;
            cfg.hpfreq         = 0.1;
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
            if k==1 || strcmp(SubID,'HN880') && k==2 % DP 29/07 - block 1 was skipped so treat k==2 as block 1
                data=dat2;
            else
                data.trial=[data.trial dat2.trial];
                data.time=[data.time dat2.time];
            end
        end
        save([preproc_path filesep 'f_etrial_ft_' SubID],'data','hdr','trl');
    end
    toc;
end