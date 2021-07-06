%%
clear all;
close all;

%%
run ../LS_RMD_localdef.m
addpath((path_fieldtrip));
ft_defaults;
addpath(genpath(path_LSCPtools));

folders=dir([data_path filesep '*_*']);


%% loop on subjects
redo=1;
        all_channels=[];
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
   
    
    numBlocks=length(files);
    fprintf('Processing %s...\n',SubID);
    tic;
    if redo==1 || exist([preproc_path filesep 'ICAcleaned_etrial_ft_' SubID '.mat'])==0
        %%% Loop on individual block files
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
            
            % fix channel names
            mylabels=hdr.label(find((cellfun(@isempty,regexp(hdr.label,'EOG')) & ~cellfun(@isempty,regexp(hdr.chantype,'eeg')))));
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
            
            if isempty(all_channels)
                all_channels=newlabels;
            else
                all_channels=intersect(all_channels,newlabels);
            end
        end
    fprintf('Common denominator has %g channels...\n',length(all_channels));
    end
end