% RTmean_perG_perB=zeros(8,3);
clearvars Y50_RTmean Y90_RTmean O90_RTmean
for nB=1:8
    Y50_RTmean{nB}=nanmean(all_nSW_perE_perB(all_nSW_perE_perB(:,2)==nB&(all_nSW_perE_perB(:,3)==1|all_nSW_perE_perB(:,3)==3),63),1);
    Y90_RTmean{nB}=nanmean(all_nSW_perE_perB(all_nSW_perE_perB(:,2)==nB&(all_nSW_perE_perB(:,3)==4),63),1);
    O90_RTmean{nB}=nanmean(all_nSW_perE_perB(all_nSW_perE_perB(:,2)==nB&(all_nSW_perE_perB(:,3)==5),63),1);
end;
Y50=cell2table(Y50_RTmean, 'VariableNames',{'Block1','Block2','Block3','Block4','Block5','Block6','Block7','Block8'})
Y90=cell2table(Y90_RTmean, 'VariableNames',{'Block1','Block2','Block3','Block4','Block5','Block6','Block7','Block8'})
O90=cell2table(O90_RTmean, 'VariableNames',{'Block1','Block2','Block3','Block4','Block5','Block6','Block7','Block8'})
T=[Y50; Y90]
T1=[T;O90]
T1.Properties.RowNames={'Young50%','Young90%','Old90%'}
mspec = ['Blocks~GroupIDs'];
rm = fitrm(T1,mspec,'WithinDesign');
ranovatbl = ranova(rm)

tab = table(GroupIDs,RT_all(:,1),RT_all(:,2),'VariableNames',{'GroupIDs','Left','Right'});
% idx = tab.GroupIDs ~= 3; tab = tab(idx,:); % eliminate a group - optional
tab.GroupIDs = nominal(tab.GroupIDs);
mspec = ['Left-Right~GroupIDs'];
within = table([1:2]','VariableNames',{'Sides'});
rm = fitrm(tab,mspec,'WithinDesign',within);
ranovatbl = ranova(rm) % for the within subjects effect
anovatbl = anova(rm) % for the between subjects effect