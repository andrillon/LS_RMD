function trl = LS_RMD_trialfun(cfg)

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
    if ~isempty(find_trials(evt_values,{'S 28'}))
        %%% discard trials just before a 28
        fixbreaks_idx=find_trials(evt_values,{'S 28'});
        for brks=length(fixbreaks_idx):-1:1
            thisbrk_idx=fixbreaks_idx(brks);
            find_stim_idx=stim_idx(find(stim_idx<thisbrk_idx));
            if ~isempty(find_stim_idx)
                find_stim_idx=find_stim_idx(end);
                stim_idx(stim_idx==find_stim_idx)=[];
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
    if ~isempty(find(ismember(evt_values,28)))
        %%% discard trials just before a 28
        fixbreaks_idx=find(ismember(evt_values,28));
        for brks=length(fixbreaks_idx):-1:1
            thisbrk_idx=fixbreaks_idx(brks);
            find_stim_idx=stim_idx(find(stim_idx<thisbrk_idx));
            if ~isempty(find_stim_idx)
                find_stim_idx=find_stim_idx(end);
            stim_idx(stim_idx==find_stim_idx)=[];
               end
     end
    end
end

%%% realign with behavioural data
if ~isempty(cfg.behav)
    behav_data=cfg.behav;
    step=1;
    start=1;
    getrid=zeros(1,length(stim_idx));
    if strformat
        thiscond=evt_values(stim_idx(start));
        thiscond=str2num(thiscond{1}(3:end));
    else
        thiscond=evt_values(stim_idx(start))-100;
    end
    while thiscond~=behav_data.trialCond(step)
        getrid(start)=1;
        start=start+1;
        if strformat
            thiscond=evt_values(stim_idx(start));
            thiscond=str2num(thiscond{1}(3:end));
        else
            thiscond=evt_values(stim_idx(start))-100;
        end
    end
    step=step+1;
    
    k=start+1;
    while k<length(stim_idx)
        if strformat
            thiscond=evt_values(stim_idx(k));
            thiscond=str2num(thiscond{1}(3:end));
        else
            thiscond=evt_values(stim_idx(k))-100;
        end
        while thiscond~=behav_data.trialCond(step) && k<length(stim_idx)
            getrid(k)=1;
            k=k+1;
            if strformat
            thiscond=evt_values(stim_idx(k));
            thiscond=str2num(thiscond{1}(3:end));
        else
            thiscond=evt_values(stim_idx(k))-100;
        end
        end
        k=k+1;
        step=step+1;
    end
    stim_idx(getrid==1)=[];
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

trl=[];
for k=1:length(stim_idx)-1
    begsample     = evt_samples(stim_idx(k)) - cfg.trialdef.prestim*hdr.Fs;
    endsample     = evt_samples(stim_idx(k+1)) + cfg.trialdef.poststim*hdr.Fs - 1;
    offset        = -cfg.trialdef.prestim*hdr.Fs;
    trl = [trl ; [round([begsample endsample offset]) stim_cond(k) stim_rt(k)]];
end
