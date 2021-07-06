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
all_nSW=[];
all_newlabels=[];
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
      fprintf('... ... %g channels\n',length(newlabels))
  all_newlabels=[all_newlabels ; newlabels];
    for nB=1:16
        nout=histc(slow_Waves(slow_Waves(:,2)==nB,3),1:64);
        all_nSW_perE=[all_nSW_perE ; [nF nB nout']];
        all_nSW(nFc,nB)=mean(nout);
    end
end
