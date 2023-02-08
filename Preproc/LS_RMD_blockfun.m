function trl = LS_RMD_blockfun(cfg)

% this function requires the following fields to be specified
% cfg.dataset
% cfg.trialdef.eventtype
% cfg.trialdef.eventvalue
% cfg.trialdef.prestim
% cfg.trialdef.poststim

hdr=ft_read_header(cfg.dataset);
evt=ft_read_event(cfg.dataset);
evt_samples=[evt.sample];
try
    evt_values={evt.value};
    evt_values(find(cellfun(@isempty,evt_values)))={'void'};
    evt_values(find(cellfun(@isnumeric,evt_values)))={'void'};
    stim_idx=find_trials(evt_values,{'S  5'});
catch
    evt_values={evt.value};
    evt_values(find(cellfun(@isempty,evt_values)))={NaN};
    evt_values(find(cellfun(@isnumeric,evt_values)))={NaN};
    evt_values=cell2mat(evt_values);
    stim_idx=find(evt_values==5);
end
if isempty(stim_idx)
    return;
end
begsample     = evt_samples(stim_idx(1)) - cfg.trialdef.prestim*hdr.Fs;
if evt_values{end}(1)=='S'
    endsample     = evt_samples(end) + cfg.trialdef.poststim*hdr.Fs - 1;
else
    endsample     = evt_samples(stim_idx(end)) + cfg.trialdef.poststim*hdr.Fs - 1;
end
offset        = -cfg.trialdef.prestim*hdr.Fs;
trl = [round([begsample endsample offset])];

