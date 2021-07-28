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
    evt_values={evt.value};
    evt_values(find(cellfun(@isempty,evt_values)))={'void'};
    stim_idx=find_trials(evt_values,{'S10.'});
    if ~isempty(find_trials(evt_values,{'28'}))
        pause;
    end
catch
    evt_values={evt.value};
    evt_values(find(cellfun(@isempty,evt_values)))={NaN};
    evt_values=cell2mat(evt_values);
    stim_idx=find(ismember(evt_values,101:106));
    if ~isempty(find(ismember(evt_values,28)))
        pause;
    end
end
trl=[];
for k=1:length(stim_idx)-1
    begsample     = evt_samples(stim_idx(k)) - cfg.trialdef.prestim*hdr.Fs;
    endsample     = evt_samples(stim_idx(k+1)) + cfg.trialdef.poststim*hdr.Fs - 1;
    offset        = -cfg.trialdef.prestim*hdr.Fs;
    trl = [trl ; [round([begsample endsample offset])]];
end
