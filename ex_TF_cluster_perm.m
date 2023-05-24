clusteralpha=0.05;
montecarloalpha=0.05;
nperm=1000;
tailFlag=0;

dat=cat(1,tf_groupA,tf_groupB); %tf_groupA: subj * freq * time
design=[ones(size(tf_groupA,1),1) ; 2*ones(size(tf_groupB,1),1)];

cfg                     = [];
cfg.neighbours          = [];
cfg.dimord              = 'rpt_freq_times';
cfg.freq                = freqP;
cfg.time                = timeP;
cfg.channel             = {'eeg'};
cfg.dim                 = [1 nFreq nTime];

cfg.numrandomization    = nperm;                   % number of randomizations, can be 'all'
cfg.correctm            = 'cluster';              % apply multiple-comparison correction, 'no', 'max', cluster', 'bonferoni', 'holms', 'fdr' (default = 'no')
cfg.alpha               = montecarloalpha;                   % critical value for rejecting the null-hypothesis per tail (default = 0.05)
cfg.clusteralpha        = clusteralpha;                   % critical value for rejecting the null-hypothesis per tail (default = 0.05)
cfg.tail                = tailFlag;                      % -1, 1 or 0 (default = 0)
cfg.correcttail         = 'no';                % correct p-values or alpha-values when doing a two-sided test, 'alpha','prob' or 'no' (default = 'no')
cfg.ivar                = 1;                      % number or list with indices, independent variable(s)
cfg.uvar                = 2;                      % number or list with indices, unit variable(s)
cfg.wvar                = [];                     % number or list with indices, within-cell variable(s)
cfg.cvar                = [];                     % number or list with indices, control variable(s)
cfg.feedback            = 'textbar';              % 'gui', 'text', 'textbar' or 'no' (default = 'text')
cfg.randomseed          = 'yes';                  % 'yes', 'no' or a number (default = 'yes')

cfg.clusterstatistic    = 'maxsum';               % no, max, maxsize, maxsum, wcm init: maxsum
cfg.clusterthreshold    = 'nonparametric_common'; % parametric, nonparametric_individual, nonparametric_common
cfg.clusteralpha        = clusteralpha;                   %
cfg.clustercritval      = [] ;                    %
cfg.clustertail         = cfg.tail;               %

cfg.statistic           = 'depsamplesT';
cfg.avgoverchan         = 'no';
[stat, cfg]             = ft_statistics_montecarlo(cfg, dat, design);

%%% ploting result
TF_diff=squeeze(mean(tf_groupA,1)-mean(tf_groupB,1));
% mytval(~stat.mask)=0;
h=simpleTFplot(TF_diff,freqP,timeP,0,0); %colormap(cmap); %ylim([freqs(1) 8])
caxis([-1 1]*max(max(abs(TF_diff))));
newmask = 1*double(reshape(stat.mask, [length(timeP) length(freqP)]))';
newmask(newmask(:) == 0) = 0.5;
alpha(h, newmask); hold on;
if ~isempty(newmask==1)
    contour(timeP, freqP,newmask,0.5,'Color','k','LineWidth',3)
end
