%%
clear all;
close all;

%%
run ../LS_RMD_localdef.m
addpath((path_fieldtrip));
ft_defaults;
addpath(genpath(path_LSCPtools));

folders=dir([data_path filesep]); % DP 06/10 removing _ to incorporate older data
folders(1:2) = []; %DP 06/10 removing hidden folders - probably more elegant solution needed here!

load(['..' filesep 'LS_RMD_Bad Trials_Channels.mat']);
load(['..' filesep 'LS_RMD_Bad_Components.mat']);


%% loop on subjects
redo=1;
for nF=1:length(folders)
    files=dir([folders(nF).folder filesep folders(nF).name filesep '*M*.eeg']); % DP 06/10 adding M to distinguish from older adults
    type_File=1;
    these_names={files.name};
    files(find(~(cellfun(@isempty,regexp(these_names,'RS')))))=[];
    SubID=folders(nF).name;
%     if ~strcmp(SubID,'HN872')
%         continue;
%     end
    if isempty(files)
        files=dir([folders(nF).folder filesep folders(nF).name filesep '*.bdf']);
        type_File=2;
        if isempty(files) %DP 06/10 adding Bryce older adult data
            files=dir([folders(nF).folder filesep folders(nF).name filesep SubID '90*.eeg']);
            type_File=3;
            if isempty(files) %DP 06/10 adding Megan older adult data
                files=dir([folders(nF).folder filesep folders(nF).name filesep 'HN9*_*.eeg']);
                type_File=4;
                if isempty(files) %DP 06/10 adding Megan younger adult data
                    files=dir([folders(nF).folder filesep folders(nF).name filesep 'HN8*_*.eeg']);
                    type_File=5;
                end
            end
        end
    end
    
      if ismember(type_File,[1 2])
        continue;
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
    
%     numBlocks=length(files);
if strcmp(SubID,'HN977')
        numBlocks=9;
    else
        numBlocks=8;
    end
    fprintf('Processing %s...',SubID);
    tic;
    if redo==1 || exist([preproc_path filesep 'ICAcleaned_etrial_ft_' SubID '.mat'])==0
        %%% Loop on individual block files
        all_channels=[];
        for k=1:numBlocks
            if type_File==1
                this_file=dir([folders(nF).folder filesep folders(nF).name filesep '*M' num2str(k) '.eeg']);
            elseif type_File==2
                this_file=dir([folders(nF).folder filesep folders(nF).name filesep '*_' num2str(k) '.bdf']);
            elseif type_File==3
                this_file=dir([folders(nF).folder filesep folders(nF).name filesep '*90' num2str(k) '.eeg']); %DP 06/10 adding Bryce older adult data
            elseif type_File==4 || type_File==5
                this_file=dir([folders(nF).folder filesep folders(nF).name filesep '*_' num2str(k) '.eeg']); %DP 06/10 adding Megan older adult + younger adult data
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
        data=[];
        all_trl=[];
        all_sigtrl=[];
        for k=1:numBlocks
            if type_File==1
                this_file=dir([folders(nF).folder filesep folders(nF).name filesep '*M' num2str(k) '.eeg']);
            elseif type_File==2
                this_file=dir([folders(nF).folder filesep folders(nF).name filesep '*_' num2str(k) '.bdf']);
            elseif type_File==3
                this_file=dir([folders(nF).folder filesep folders(nF).name filesep '*90' num2str(k) '.eeg']); %DP 06/10 adding Bryce older adult data
            elseif type_File==4 || type_File==5
                this_file=dir([folders(nF).folder filesep folders(nF).name filesep '*_' num2str(k) '.eeg']); %DP 06/10 adding Megan older adult + younger adult data
            end
            if isempty(this_file)
                continue;
            end
            if type_File==1
                behavfiles=dir([folders(nF).folder filesep 'Behav' filesep this_file.name(1:end-4) '.mat']);
            elseif type_File==2
                behavfiles=dir([folders(nF).folder filesep 'Behav' filesep this_file.name(1:end-4) '.mat']);
            elseif type_File==3
                behavfiles=dir([folders(nF).folder filesep 'Behav' filesep SubID '90'  num2str(k) '.mat']);
            elseif type_File==4 || type_File==5
                behavfiles=dir([folders(nF).folder filesep 'Behav' filesep this_file.name(1:end-4) '.mat']); %DP 28/07 adding Megan older adult + younger adult data
            end
            
            file_name = this_file(1).name;
            file_folder = this_file(1).folder;
            FileID=file_name(1:end-4);
            fprintf('... working on subject %s (%g/%g) - file %s (%g/%g)\n',SubID,nF,length(folders),FileID,k,length(files))
            
             if exist('behavfiles')~=0 && ~isempty(behavfiles)
               behav_data=load([behavfiles.folder filesep behavfiles.name]);
            else
                behav_data=[];
                warning('cannot find the behavioural data!!');
             end
 
            %%% Read headers
            hdr=ft_read_header([file_folder filesep file_name]);
            
            %%% Define epochs
            cfg=[];
            cfg.trialfun            = 'LS_RMD_trialfun_v3';
            cfg.SubID               = SubID;
            cfg.dataset             = [file_folder filesep file_name];
            cfg.trialdef.prestim    = 1;
            cfg.trialdef.poststim   = 1;
            cfg.behav               = behav_data;
            cfg.type_File           = type_File;
%             cfg = ft_definetrial(cfg);
%             cfg.trl(cfg.trl(:,2)>hdr.nSamples,:)=[];
            try
                cfg = ft_definetrial(cfg);
            catch
                cfg.trl=[]; trl=[];
                warning('problem with trial definition (no trial defined)');
                all_trl=[all_trl ; [nF k 0 nan(1,2)]];
                continue;
            end
            trl=cfg.trl;
            all_trl=[all_trl ; [nF k size(trl,1) mean(trl(:,5:6),1)]];
            if size(trl,1)<10
                warning('problem with trial definition (less than 10 trials defined)');
                continue;
            end
            
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
            if isempty(data)
                data=dat2;
                all_sigtrl=[all_sigtrl ; [k*ones(size(trl,1),1) trl]];
            else
                data.trial=[data.trial dat2.trial];
                data.time=[data.time dat2.time];
                all_sigtrl=[all_sigtrl ; [k*ones(size(trl,1),1) trl]];
%                 fprintf('!!!!!!! block %g %g trials in trl and %g in data\n',k,size(trl,1),length(dat2.trial))
            end
        end
        
        
        %%% take out trial
        thisF=match_str(LSRMDBadTrialsChannels.SubID,SubID);
        if isempty(thisF)
            continue;
        end
        if isempty(data)
            continue;
        end
        eval(['badChannels=[' LSRMDBadTrialsChannels.ExcludedChannels{thisF} '];']);
        eval(['badTrials=[' LSRMDBadTrialsChannels.ExcludedTrials{thisF} '];']);
        cfg=[];
        cfg.trials          = setdiff(1:length(data.trial),badTrials);
        data = ft_preprocessing(cfg, data);
        badTrials(badTrials>length(all_sigtrl))=[]; % Removing bad trials that are in block 9
        all_sigtrl(badTrials,:)=[];
        
        %%% Layout
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
        cfg = [];
        cfg.layout = 'acticap-64ch-standard2.mat';
        cfg.channel=newlabels;
        cfg.center      = 'yes';
        layout=ft_prepare_layout(cfg);
        
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
            cfg.badchannel     = layout.label(badChannels);
            cfg.missingchannel = [];
            cfg.neighbours     = neighbours;
            cfg.trials         = 'all';
            cfg.layout         = layout;
            cfg.channel = layout.label;
            
            data.label=newlabels;
            [data] = ft_channelrepair(cfg, data);
        end
        
        %         %%% CSD
        %         elec=[];
        %         elec.label=layout.label;
        %         elec.pos=layout.pos;
        %         cfg=[];
        %         cfg.method       = 'spline';
        %         cfg.elec         = elec;
        %         cfg.order        = 4;
        %         cfg.degree       = 14;
        %         [data] = ft_scalpcurrentdensity(cfg, data);
        
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
        
%         cfg           = [];
%         cfg.toilim    = [-0.7 1.8];
%         [data] = ft_redefinetrial(cfg, data);
        
        cfg=[];
        cfg.demean          = 'yes';
        cfg.baselinewindow  = [-0.1 0];
        data = ft_preprocessing(cfg,data);
        
        %%%%%% Re-order channels
        load(['..' filesep 'LS_RMD_Common_Electrodes_64ch.mat']);
        % fix channel names
        mylabels=data.label;
        newchanidx=[];
        newlabels=[];
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
        for  nCh=1:length(all_channels)
            newchanidx(nCh)=find(ismember(newlabels,all_channels{nCh}));
        end
        for  nTr=1:length(data.trial)
            data.trial{nTr}=data.trial{nTr}(newchanidx,:);
        end
        data.label=all_channels;
        
        save([preproc_path filesep 'ICAcleaned_etrial_ft_' SubID],'data','hdr');
        table=array2table(all_sigtrl,'VariableNames',{'Block','Stat','End','Offset','Cond','RT','DiffOnset'});
        if size(table,1)~=length(data.trial)
            warning('NOT THE SAME SIZE BETWEEN DATA AND BEHAV')
        end
        table.RT((table.RT*hdr.Fs+table.Stat)>table.End)=NaN;
        writetable(table,[preproc_path filesep 'behav_' SubID '.csv']);
        Duration=(table.End-table.Stat)/hdr.Fs;
        fprintf('... ... %g trials detected: min duration %g | max duration %g | average duration %g\n',length(Duration),min(Duration),max(Duration),mean(Duration))
        cfgerp        = [];
        cfgerp.trials = 1:length(data.trial);
        av_ERP = ft_timelockanalysis(cfgerp, data);
        
    end
    toc;
end
