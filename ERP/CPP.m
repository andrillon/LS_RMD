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
        
        fprintf('%2.0f-%2.0f\n',0,0)

        % crop trials
        cfg=[];
        cfg.toilim =[0 1.0];
        data=ft_redefinetrial(cfg,data);
        
        % remove weird extra field
        if isfield(data,'elec')
            data=rmfield(data,'elec')
        end
        
        % estimate the ERP
        tlck(nF) = ft_timelockanalysis([],data);
        
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

for pps=1:nFc
    PzERP(pps,:)=tlck(pps).avg(strcmp(tlck(pps).label,'Pz'),:);
end

% young_CPP=nanmean(PzERP(all_agegroup==0,:));
% old_CPP=nanmean(PzERP(all_agegroup==1,:));

t=tlck(1).time;

figure; 
[~,hp(1)]=simpleTplot(t,PzERP(all_agegroup==0,:),0,'r',[0 0.05 0.0001 1000],'-',0.1,1,0,1,1); 
[~,hp(2)]=simpleTplot(t,PzERP(all_agegroup==1,:),0,'b',[0 0.05 0.0001 1000],'-',0.1,1,0,1,1);
legend(hp,'Younger','Older','Location','northeast'); xlim([0 1]); hp(1).LineWidth=3; hp(2).LineWidth=3; xticks([0 .25 .50 .75 1.00]);
xlabel('Post-target (s)'); set(gcf,'Color','white'); title('CPP (Pz)');