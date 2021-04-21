%%
clear all;
close all;

%%
run ../LS_RMD_localdef.m
addpath((path_fieldtrip));
ft_defaults;

folders=dir([data_path filesep '*_*']);

%% loop on subjects
for nF=1:length(folders)
    files=dir([folders(nF).folder filesep folders(nF).name filesep '*.eeg']);
    these_names={files.name};
    files(find(~(cellfun(@isempty,regexp(these_names,'_')))))=[];
    SubID=folders(nF).name;
    
    tic;
    if exist([preproc_path filesep 'f_eblock_ft_' SubID '.mat'])==0
        %%% Loop on individual block files
        all_channels=[];
        for k=1:length(files)
            file_name = files(k).name;
            file_folder = files(k).folder;
            hdr=ft_read_header([file_folder filesep file_name]);
            if k==1
                all_channels=hdr.label(find((cellfun(@isempty,regexp(hdr.label,'EOG')))));
            else
                all_channels=intersect(all_channels,hdr.label(find((cellfun(@isempty,regexp(hdr.label,'EOG'))))));
            end
        end
        
        for k=1:length(files)
            file_name = files(k).name;
            file_folder = files(k).folder;
            FileID=file_name(1:end-4);
            fprintf('... working on subject %s (%g/%g) - file %s (%g/%g)\n',SubID,nF,length(folders),FileID,k,length(files))
            
            %%% Read headers
            hdr=ft_read_header([file_folder filesep file_name]);
            
            %%% Define epochs
            cfg=[];
            cfg.trialfun            = 'LS_RMD_blockfun';
            cfg.SubID               = SubID;
            cfg.dataset             = [file_folder filesep file_name];
            cfg.trialdef.prestim    = 0.5;
            cfg.trialdef.poststim   = 0.5;
            cfg = ft_definetrial(cfg);
            
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
            if k==1
                data=dat2;
            else
                data.trial{k}=dat2.trial{1};
                data.time{k}=dat2.time{1};
            end
        end
        save([preproc_path filesep 'f_eblock_ft_' SubID],'data','hdr');
    end
    toc;
end