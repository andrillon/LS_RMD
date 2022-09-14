%% PC of recordings
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
pathFit = 'G:\Local Sleep Analysis\LSRMD\LS_RMD\PCA\FittingResults\';

% The percentage of variance to be explained by the extracted PCs
pVar = 99;

nFc=0; nFo=0; nFy=0;

%% Load files and calculate group averages
for nF=1:length(files)
    file_name = files(nF).name;
    folder_name = files(nF).folder;
    SubID=file_name(1:end-4);
    seps=findstr(SubID,'ft_');
    SubID=SubID(seps(1)+3:end);
    tic;
    fprintf('... working on %s (%g/%g)\n',SubID,nF,length(files))
    
    load([folder_name filesep file_name]);
    
    nFc=nFc+1;
    
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
        cfg.toilim =[-0.5 1];
        data=ft_redefinetrial(cfg,data);
        
        temp_data=data.trial;
        temp_data_cat=cat(3,temp_data{:});
        temp_data=mean(temp_data_cat,3);
             
        if length(SubID)==4 && SubID(1)=='A' % OLD (MONASH) - UP & DOWN 90%COH
            thisgroup=2;
            thisagegroup=NaN;
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
        all_agegroup(nFc)=thisagegroup; % STORE HERE THE INFO ABOUT THE AGE GROUP
        
        if thisgroup==4
            nFy=nFy+1;
            young90_subjdata{nFy}=temp_data;
        elseif thisgroup==5
            nFo=nFo+1;
            old90_subjdata{nFo}=temp_data;
        else
            continue;
        end
end
young90_subjdata=cat(3,young90_subjdata{:});
young90_subjdata=young90_subjdata(:,127:end,:);
young90_mean=mean(young90_subjdata(:,:,:),3);
young90_mean=young90_mean';
old90_subjdata=cat(3,old90_subjdata{:});
old90_subjdata=old90_subjdata(:,127:end,:);
old90_mean=mean(old90_subjdata(:,:,:),3);
old90_mean=old90_mean';

%% Run PCA
[young90_coeff,young90_score,young90_latent,young90_tsquared,young90_explained] = pca(young90_mean);
[old90_coeff,old90_score,old90_latent,old90_tsquared,old90_explained] = pca(old90_mean);

%% Find the components which explain more than pVar% of the variance 
% First for younger adults
young90i = 1;
young90a = young90_explained(young90i);
while young90a <= pVar
    young90a = young90a + young90_explained(young90i+1);
    young90i = young90i+1;
end
young90_nDim = young90i;

% Then older adults
old90i = 1;
old90a = old90_explained(old90i);
while old90a <= pVar
    old90a = old90a +  old90_explained(old90i+1);
    old90i = old90i+1;
end
old90_nDim = old90i;

%% Weight each individual's data according to the selected components  
for idx = 1:size(young90_subjdata,3) % First with younger group 
    young90_allsubj{idx} = young90_subjdata(:,:,idx);
    young90_subjdata_mean = mean(young90_allsubj{idx},2); % Find mean of data (columns)
    young90_subjdata_rep = repmat(young90_subjdata_mean,1,length(young90_allsubj{idx})); %Replicate mean vector to matrix for subtraction
    young90_subjdata_norm{idx} = young90_allsubj{idx} - young90_subjdata_rep; % Normalize data to zero mean y subtraction
    young90_subjScores{idx} = young90_subjdata_norm{idx}'*young90_coeff(:,1:young90_nDim); % Manually calculate scores using PCA coeff and normalized data
    young90_allSubjScores (idx, :,:) =  young90_subjScores{idx};
end
young90_meanSubjScores = squeeze(mean(young90_allSubjScores,1));

% Then with older group
for idx = 1:size(old90_subjdata,3) % Then older
    old90_allsubj{idx} = old90_subjdata(:,:,idx);
    old90_subjdata_mean = mean(old90_allsubj{idx},2); % Find mean of data (columns)
    old90_subjdata_rep = repmat(old90_subjdata_mean,1,length(old90_allsubj{idx})); %Replicate mean vector to matrix for subtraction
    old90_subjdata_norm{idx} = old90_allsubj{idx} - old90_subjdata_rep; % Normalize data to zero mean y subtraction
    old90_subjScores{idx} = old90_subjdata_norm{idx}'*old90_coeff(:,1:old90_nDim); % Manually calculate scores using PCA coeff and normalized data
    old90_allSubjScores (idx, :,:) =  old90_subjScores{idx};
end
old90_meanSubjScores = squeeze(mean(old90_allSubjScores,1));

save([pathFit 'PCA_clusters']);

%% scalp maps : Plot the components which explain more than pVar% of the variance - Younger Adults
% load([pathFit 'PCA_clusters']);
% 
% close all;
% F = figure;
% set(F,'color','w');
% % s = title(['PCs explaining more than ' , num2str(pVar),'% of variance']);
% % s.FontSize = 28;
% % gs = get(s,'position');
% % set(s,'position',[gs(1) gs(2)-0.2 gs(3)])
% for k = 1:young90_nDim
%     subplot(1,young90_nDim,k)
%     grandAverage = [];
%     grandAverage.label = newlabels;
%     grandAverage.time = [k];
%     grandAverage.dimord = 'chan_time';
%     grandAverage.avg = young90_coeff(:,k);
%     cfg = [];
% %   cfg.zlim =[-0.3 0.3];
%     cfg.parameter = 'avg';
%     cfg.markersize = 2;
%     cfg.comment = 'no';
%     cfg.layout = 'biosemi64.lay';
%     cfg.colormap = 'hot';
%     ft_topoplotER(cfg, grandAverage);
%     t = title([num2str(young90_explained(k)),'%'], 'FontSize', 14);
% end
% 
% %% Time series : Plot the components which explain more than pVar% of the variance - Younger Adults
% for n = 1:young90_nDim
%     subplot(1,young90_nDim,n)
%     p = plot (1:length(young90_meanSubjScores), young90_meanSubjScores(:,n),'r');
%     box off
%     p.LineWidth = 2;
%     ylim([-30 30]);
%     % yticks(-2:1:2);
%     % xticks(-100:100:400);
%     xlabel('Time (ms)');
%     ylabel('Amplitude (\muV)');
% %     v = vline(0, '--k');
% %     v.LineWidth =  1;
%     set(gcf,'color','w')
%     set(gca,'FontSize',26,'FontWeight','bold','linewidth',2)
%     set(gca,'TickDir','out')
% end

%% scalp maps : Plot the components which explain more than pVar% of the variance - Older Adults
load([pathFit 'PCA_clusters']);

close all;
F = figure;
set(F,'color','w');
% s = title(['PCs explaining more than ' , num2str(pVar),'% of variance']);
% s.FontSize = 28;
% gs = get(s,'position');
% set(s,'position',[gs(1) gs(2)-0.2 gs(3)])
for k = 1:old90_nDim
    subplot(1,old90_nDim,k)
    grandAverage = [];
    grandAverage.label = newlabels;
    grandAverage.time = [k];
    grandAverage.dimord = 'chan_time';
    grandAverage.avg = old90_coeff(:,k);
    cfg = [];
%   cfg.zlim =[-0.3 0.3];
    cfg.parameter = 'avg';
    cfg.markersize = 2;
    cfg.comment = 'no';
    cfg.layout = 'biosemi64.lay';
    cfg.colormap = 'hot';
    ft_topoplotER(cfg, grandAverage);
    t = title([num2str(old90_explained(k)),'%'], 'FontSize', 14);
end

%% Time series : Plot the components which explain more than pVar% of the variance - Younger Adults
for n = 1:old90_nDim
    subplot(1,old90_nDim,n)
    p = plot (1:length(old90_meanSubjScores), old90_meanSubjScores(:,n),'r');
    box off
    p.LineWidth = 2;
    ylim([-30 30]);
    % yticks(-2:1:2);
    % xticks(-100:100:400);
    xlabel('Time (ms)');
    ylabel('Amplitude (\muV)');
%     v = vline(0, '--k');
%     v.LineWidth =  1;
    set(gcf,'color','w')
    set(gca,'FontSize',26,'FontWeight','bold','linewidth',2)
    set(gca,'TickDir','out')
end