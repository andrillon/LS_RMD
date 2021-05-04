%%
close all;
clear all;

%%
run ../LS_RMD_localdef.m
% adding relevant toolboxes to the path
addpath((path_fieldtrip));
ft_defaults;

%% choose and load subject
List_Subj=dir([preproc_path filesep 'f_etrial_ft_*.mat']);
ListNames={List_Subj.name};
pick=listdlg('ListString',ListNames);
load([preproc_path filesep ListNames{pick}])

%% reject trials
cfg          = [];
cfg.method   = 'summary';
cfg.alim     = 5e-5;
data2        = ft_rejectvisual(cfg,data);

%% display data
cfg=[];
cfg.continuous='no';
cfg.allowoverlap='true';
cfg.viewmode='vertical';
cfg = ft_databrowser(cfg, data);

