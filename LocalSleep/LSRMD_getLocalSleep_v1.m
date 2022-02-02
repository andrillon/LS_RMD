%%
clear all
close all

run ../LS_RMD_localdef.m

addpath(path_fieldtrip);
ft_defaults;
addpath(genpath(path_LSCPtools));

files=dir([preproc_path filesep 'ICAcleaned_eblock_ft_*.mat']);

%%
redo=0;
for nF=1:length(files)
    file_name = files(nF).name;
    folder_name = files(nF).folder;
    SubID=file_name(1:end-4);
    seps=findstr(SubID,'ft_');
    SubID=SubID(seps(1)+3:end);
    tic;
    fprintf('... working on %s (%g/%g)\n',SubID,nF,length(files))
    
    
    
    if redo==1 || exist([preproc_path filesep '..' filesep 'SWdetection' filesep 'allSW_' SubID '.mat'])==0
         load([folder_name filesep file_name]);
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
        
        all_Waves=[];
        for nBl=1:size(data.trial,2)
            temp_data=data.trial{nBl}(:,:);
            temp_data=temp_data-repmat(mean(temp_data(match_str(newlabels,{'TP7','TP8'}),:),1),size(temp_data,1),1);
            temp_data=temp_data-repmat(mean(temp_data,2),1,size(temp_data,2));
            
            [twa_results]=twalldetectnew_TA_v2(temp_data,data.fsample,0);
            for nE=1:size(data.label,1)
                all_Waves=[all_Waves ; [repmat([nF nBl nE],length(abs(cell2mat(twa_results.channels(nE).maxnegpkamp))),1) abs(cell2mat(twa_results.channels(nE).maxnegpkamp))'+abs(cell2mat(twa_results.channels(nE).maxpospkamp))' ...
                    cell2mat(twa_results.channels(nE).negzx)' ...
                    cell2mat(twa_results.channels(nE).poszx)' ...
                    cell2mat(twa_results.channels(nE).wvend)' ...
                    cell2mat(twa_results.channels(nE).maxnegpk)' ...
                    cell2mat(twa_results.channels(nE).maxnegpkamp)' ...
                    cell2mat(twa_results.channels(nE).maxpospk)' ...
                    cell2mat(twa_results.channels(nE).maxpospkamp)' ...
                    cell2mat(twa_results.channels(nE).mxdnslp)' ...
                    cell2mat(twa_results.channels(nE).mxupslp)' ...
                    cell2mat(twa_results.channels(nE).maxampwn)' ...
                    cell2mat(twa_results.channels(nE).minampwn)' ...
                    ]];
            end
        end
        fprintf('\n')
        save([preproc_path filesep '..' filesep 'SWdetection' filesep 'allSW_' SubID],'all_Waves')
        
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
        
        thr_Wave=[];
        slow_Waves=[];
        for nE=1:size(data.label,1)
            thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
            temp_p2p=thisE_Waves(:,paramSW.AmpCriterionIdx);
            
            if ~isempty(paramSW.fixThr)
                thr_Wave(nE)=paramSW.fixThr;
            else
                thr_Wave(nE)=prctile(thisE_Waves(:,paramSW.AmpCriterionIdx),paramSW.prticle_Thr);
            end
            slow_Waves=[slow_Waves ; thisE_Waves(temp_p2p>thr_Wave(nE),:)];
        end
        save([preproc_path filesep '..' filesep 'SWdetection' filesep 'SW_' SubID],'slow_Waves','paramSW')
    end
end
