%%
close all;
clear all;

%%
run ../LS_RMD_localdef.m
% adding relevant toolboxes to the path
addpath((path_fieldtrip));
ft_defaults;

%% choose and load subject
List_Subj=dir([preproc_path filesep 'ICf_etrial_ft_*.mat']);
ListNames={List_Subj.name};
pick=listdlg('ListString',ListNames);
load([preproc_path filesep ListNames{pick}])

%% Layout
mylabels=data.label;
for nCh=1:length(mylabels)
    findspace=findstr(mylabels{nCh},' ');
    if ismember(mylabels{nCh}(1),{'1','2','3','4','5','6','7','8','9'})
        newlabels{nCh}=mylabels{nCh}(findspace+1:end);
    else
        newlabels{nCh}=mylabels{nCh}(1:findspace-1);
    end
end

cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.channel=newlabels;
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);

comp.topolabel=newlabels;
data.label=newlabels;
%%
figure;
cfg = [];
cfg.component = 1:30;       % specify the component(s) that should be plotted
cfg.layout    = layout; % specify the layout file that should be used for plotting
cfg.comment   = 'no';
%     cfg.marker = 'labels';
ft_topoplotIC(cfg, comp)
set(gcf,'Position',[1           1        1871         984]);

figure;
cfg = [];
cfg.component = 31:length(comp.label);       % specify the component(s) that should be plotted
cfg.layout    = layout; % specify the layout file that should be used for plotting
cfg.comment   = 'no';
%     cfg.marker = 'labels';
ft_topoplotIC(cfg, comp)
set(gcf,'Position',[1           1        1871         984]);


%% reject trials
pickComponents=input('Select component you want to plot: ');

%%
cfg          = [];
cfg.channel  = pickComponents; % components to be plotted
cfg.viewmode = 'component';
cfg.layout   = layout; % specify the layout file that should be used for plotting
cfg.allowoverlap='true';
cfg.continuous='no';
ft_databrowser(cfg, comp);
fprintf('... working on %s\n',ListNames{pick})

