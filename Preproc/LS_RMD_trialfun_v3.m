function trl = LS_RMD_trialfun_v2(cfg)
trl=[];

% this function requires the following fields to be specified
% cfg.dataset
% cfg.trialdef.eventtype
% cfg.trialdef.eventvalue
% cfg.trialdef.prestim
% cfg.trialdef.poststim

hdr=ft_read_header(cfg.dataset);
evt=ft_read_event(cfg.dataset);
evt_samples=[evt.sample];
% Duration range: 3.06 to 7.29s
% 101: Coherent Motion
% 102: Coherent Motion
% 103: Coherent Motion
% 104: Coherent Motion
% 105: Coherent Motion
% 106: Coherent Motion
% Fixation break trigger: 28 (get rid of trials with a 28 between the S1**
% and S5 or S12)
try
    strformat=1;
    evt_values={evt.value};
    evt_values(find(cellfun(@isempty,evt_values)))={'void'};
    possible_values={'S101','S102','S103','S104','S105','S106','S107','S108','S109','S110',...
        'S111','S112'};
    stim_idx=find(ismember(evt_values,possible_values));
    
    offset_idx=find(ismember(evt_values,'S  5'));
    all_breaks=[];
    if ~isempty(find_trials(evt_values,{'S 28'}))
        %%% discard trials just before a 28
        fixbreaks_idx=find_trials(evt_values,{'S 28'});
        for brks=length(fixbreaks_idx):-1:1
            thisbrk_idx=fixbreaks_idx(brks);
            find_stim_idx=stim_idx(find(stim_idx<thisbrk_idx));
            if ~isempty(find_stim_idx)
                find_stim_idx=find_stim_idx(end);
                if ~isempty(find(ismember(evt_values(find_stim_idx:thisbrk_idx),{'S  4','S  5','S 12'})))
                else
                    all_breaks=[all_breaks find(stim_idx==find_stim_idx)];
                    stim_idx(stim_idx==find_stim_idx)=[];
                end
            end
        end
    end
catch
    strformat=0;
    evt_values={evt.value};
    evt_values(find(cellfun(@isempty,evt_values)))={NaN};
    evt_values=cell2mat(evt_values);
    possible_values=101:112;
    stim_idx=find(ismember(evt_values,possible_values));
    all_breaks=[];
    if ~isempty(find(ismember(evt_values,28)))
        %%% discard trials just before a 28
        fixbreaks_idx=find(ismember(evt_values,28));
        for brks=length(fixbreaks_idx):-1:1
            thisbrk_idx=fixbreaks_idx(brks);
            find_stim_idx=stim_idx(find(stim_idx<thisbrk_idx));
            if ~isempty(find_stim_idx)
                find_stim_idx=find_stim_idx(end);
                if ~isempty(find(ismember(evt_values(find_stim_idx:thisbrk_idx),[4 5 12])))
                else
                    all_breaks=[all_breaks find(stim_idx==find_stim_idx)];
                    stim_idx(stim_idx==find_stim_idx)=[];
                end
            end
        end
    end
end

%%% realign with behavioural data
if ~isempty(cfg.behav)
    behav_data=cfg.behav;
    realigned_trials=[];
    clean_realigned_trials=[];
    PTBtrig=cfg.behav.PTBtrig(cfg.behav.PTBtrig~=5)-100;
    if max(all_breaks)>length(PTBtrig)
        %                     warning(sprintf('PTB times do not match condition (%g vs %g)',length(PTBtrig),length(behav_data.trialCond)));
        %             return;
        all_breaks(all_breaks>length(PTBtrig))=[];
    end
    PTBtrig(all_breaks)=[];
    if length(PTBtrig)==length(behav_data.trialCond) && sum(PTBtrig==behav_data.trialCond)==length(behav_data.trialCond)
        PTBtimes=cfg.behav.PTBtrigT(cfg.behav.PTBtrig~=5)-100;
        PTBtimes(all_breaks)=[];
    elseif length(behav_data.trialCond)-1==length(PTBtrig)
        trialCond=behav_data.trialCond;
        sumcheck=[];
        for m=1:length(trialCond)
            sumcheck(m)=sum(PTBtrig==trialCond(setdiff(1:length(trialCond),m)))==length(trialCond(setdiff(1:length(trialCond),m)));
        end
        suppressed_trial=find(sumcheck);
        if length(suppressed_trial)==1
            behav_data.trialCond(suppressed_trial)=[];
        end
        PTBtimes=cfg.behav.PTBtrigT(cfg.behav.PTBtrig~=5)-100;
        PTBtimes(all_breaks)=[];
    elseif length(behav_data.trialCond)==length(PTBtrig)-1     && PTBtrig(1)==PTBtrig(end)
        PTBtrig(end)=[];
        if sum(PTBtrig==behav_data.trialCond)==length(behav_data.trialCond)
            PTBtimes=cfg.behav.PTBtrigT(cfg.behav.PTBtrig~=5)-100;
            PTBtimes(all_breaks)=[];
            PTBtimes(end)=[];
        else
            warning(sprintf('PTB times do not match condition (%g vs %g)',length(PTBtrig),length(behav_data.trialCond)));
            return;
        end
    else
        warning(sprintf('PTB times do not match condition (%g vs %g)',length(PTBtrig),length(behav_data.trialCond)));
        return;
    end
    EVTtimes=evt_samples(stim_idx);
    
    step=1;
    start=1;
    getrid=zeros(1,length(stim_idx));
    if strformat
        thiscond=evt_values(stim_idx(start));
        thiscond=str2num(thiscond{1}(3:end));
    else
        thiscond=evt_values(stim_idx(start))-100;
    end
    while thiscond~=behav_data.trialCond(step) && start<length(stim_idx)
        getrid(start)=1;
        start=start+1;
        if strformat
            thiscond=evt_values(stim_idx(start));
            thiscond=str2num(thiscond{1}(3:end));
        else
            thiscond=evt_values(stim_idx(start))-100;
        end
    end
    timeOri_EVT=evt_samples(stim_idx(start));
    timeOri_PTB=PTBtimes(step);
    realigned_trials=[step stim_idx(start) 0 0 0 thiscond];
    
    EVTtimes=(EVTtimes-timeOri_EVT)/hdr.Fs;
    PTBtimes=PTBtimes-timeOri_PTB;
    
    for k=2:length(PTBtimes)
        thiscond=behav_data.trialCond(k);
        if strformat
            these_stimidx=find(ismember(evt_values(stim_idx),['S' num2str(thiscond+100)]));
        else
            these_stimidx=find(ismember(evt_values(stim_idx),(thiscond+100)));
        end
        if ~isempty(these_stimidx)
            [closesttime, closestindex]=findclosest(EVTtimes(these_stimidx),PTBtimes(k));
            realigned_trials(k,:)=[k stim_idx(these_stimidx(closestindex)) PTBtimes(k) closesttime PTBtimes(k)-closesttime thiscond];
        else
            realigned_trials(k,:)=[k NaN PTBtimes(k) NaN +Inf thiscond];
        end
    end
    
    clean_realigned_trials=realigned_trials(abs(realigned_trials(:,5))<0.1,:);
    
    stim_idx=clean_realigned_trials(:,2);
    
    fprintf('... ... realigned %g out of %g trials with a mean absolute jitter of %1.5f s\n',length(stim_idx),size(realigned_trials,1),mean(abs(clean_realigned_trials(:,5))));
end


%%% find trial condition
% For the first cohort of younger adults (for reference, these are denoted as type_File=1:2 in the import_EEG script) the following applies for coherent motion triggers:
% S 101:S 106.
% Left target trials: S 101, S 103, S 105.
% Right target trials: S 102, S 104, S 106.
% For the second cohort of younger adults and the older adults (type_File=3:5) the following are the coherent motion triggers:
% S 101:S 112.
% Left target trials: S 101, S 102, S 105, S 106, S 109, S 110.
% Right target trials: S 103, S 104, S 107, S 108, S 111, S 112.
if ismember(cfg.type_File,3:5)
    if strformat
        left_targets={'S101','S102','S105','S106','S109','S110'};
        right_targets={'S103','S104','S107','S108','S111','S112'};
    else
        left_targets=[101 102 105 106 109 110];
        right_targets=[103 104 107 108 111 112];
    end
elseif ismember(cfg.type_File,1:2)
    if strformat
        left_targets={'S101','S103','S105'};
        right_targets={'S102','S104','S106'};
    else
        left_targets=[101 103 105];
        right_targets=[102 104 106];
    end
else
    left_targets=[];
    right_targets=[];
    warning('unexpected number of evt values');
end
stim_cond=nan(1,length(stim_idx));
stim_cond(find(ismember(evt_values(stim_idx),left_targets)))=1;
stim_cond(find(ismember(evt_values(stim_idx),right_targets)))=2;

%%% find subjects' respones
if strformat
    code_resp='S 12';
else
    code_resp=12;
end
resp_idx=find(ismember(evt_values,code_resp));
resp_sample=evt_samples(resp_idx);
stim_sample=evt_samples(stim_idx);
stim_rt=nan(1,length(stim_idx));
for rsp=1:length(resp_sample)
    thisresp_sample=resp_sample(rsp);
    find_stim_idx=stim_idx(find(stim_sample<thisresp_sample));
    if ~isempty(find_stim_idx)
        find_stim_idx=find_stim_idx(end);
        this_rt=(thisresp_sample-stim_sample(stim_idx==find_stim_idx))/hdr.Fs;
        stim_rt(stim_idx==find_stim_idx)=this_rt;
        %     else
        %         this_rt=nan;
    end
end

for k=1:length(stim_idx)-1
    begsample     = evt_samples(stim_idx(k)) - cfg.trialdef.prestim*hdr.Fs;
    temp_offset=offset_idx(offset_idx>stim_idx(k));
    if isempty(temp_offset)
        continue;
    end
    temp_offset=temp_offset(1);
    endsample     = evt_samples(temp_offset) + cfg.trialdef.poststim*hdr.Fs - 1;
    offset        = -cfg.trialdef.prestim*hdr.Fs;
    trl = [trl ; [round([begsample endsample offset]) stim_cond(k) stim_rt(k)]];
end
k=length(stim_idx);
begsample     = evt_samples(stim_idx(k)) - cfg.trialdef.prestim*hdr.Fs;
endsample     = min(hdr.nSamples,begsample+min(trl(:,2)-trl(:,1)));
offset        = -cfg.trialdef.prestim*hdr.Fs;
trl = [trl ; [round([begsample endsample offset]) stim_cond(k) stim_rt(k)]];

trl=[trl abs(clean_realigned_trials(:,5))];