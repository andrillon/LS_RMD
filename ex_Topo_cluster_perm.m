clusteralpha=0.05;
montecarloalpha=0.05;
nperm=1000;
tailFlag=0;
layout=[];

dat=cat(1,tf_groupA,tf_groupB); %tf_groupA: subj * freq * time
design=[ones(size(tf_groupA,1),1) ; 2*ones(size(tf_groupB,1),1)];


cfg_neighb=[];
cfg_neighb.method = 'template';
cfg_neighb.layout=layout;
cfg_neighb.channel=layout.label;
neighbours = ft_prepare_neighbours(cfg_neighb);

[SW_clus]=get_clusterperm_lme_lsneurom(SW_est,clusteralpha,montecarloalpha,nperm,neighbours);