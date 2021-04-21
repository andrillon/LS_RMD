function trl = LS_RMD_blockfun(cfg)

% this function requires the following fields to be specified
% cfg.dataset
% cfg.trialdef.eventtype
% cfg.trialdef.eventvalue
% cfg.trialdef.prestim
% cfg.trialdef.poststim

hdr=ft_read_header(cfg.dataset);
evt=ft_read_event(cfg.dataset);
evt_values={evt.value};
evt_samples=[evt.sample];
evt_values(find(cellfun(@isempty,evt_values)))={'void'};
stim_idx=find_trials(evt_values,{'S'});

begsample     = evt_samples(stim_idx(1)) - cfg.trialdef.prestim*hdr.Fs;
endsample     = evt_samples(stim_idx(end)) + cfg.trialdef.poststim*hdr.Fs - 1;
offset        = -cfg.trialdef.prestim*hdr.Fs;
trl = [round([begsample endsample offset])];

