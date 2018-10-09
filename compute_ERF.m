% Computes ERFs & covariance matrices for by-condition data, combined cue & combined target

function [erf_clean, erf_cue_combined, erf_target_combined] = compute_ERF (trials_clean)

    % run the #define section
    global conds_cue; global conds_target; global eventnames;
    common();
    
    % We kept bad trials as NaN, so need to exclude them now
    for j = 1:length(eventnames)
        %length(trials_clean.(eventnames{j}).trial)
        good_trials_idx = find(~isnan(cell2mat(cellfun(@(isgood)isgood(1),trials_clean.(eventnames{j}).trial,'uni',0)))); % Just need to evaluate the first element as all samples in bad trial are NaN
        cfg             = [];
        cfg.trials      = good_trials_idx;
        trials_clean.(eventnames{j}) = ft_redefinetrial(cfg, trials_clean.(eventnames{j}));
    end

    % An alternative method to remove NaN trials - doesn't seem to work!
    %{
    for j = 1:length(eventnames)
        cfg         = [];
        cfg.nanmean = 'yes';
        trials_clean.(eventnames{j}) = ft_selectdata(cfg, trials_clean.(eventnames{j})); % Do this because we kept bad trials as NaN
    end
    %}

    % Cue Window
    % Compute erf & cov matrix on the combined data (cue window in all 4 conditions appended together)
    cfg = [];
    cue_combined = ft_appenddata(cfg, trials_clean.cuechstay, trials_clean.cuechswitch, trials_clean.cueenstay, trials_clean.cueenswitch);
    cfg.covariance       = 'yes';
    cfg.covariancewindow = [0 0.65];
    erf_cue_combined = ft_timelockanalysis(cfg, cue_combined);
    % Compute erf & cov matrix for each condition 
    for j = conds_cue
        erf_clean.(eventnames{j}) = ft_timelockanalysis(cfg, trials_clean.(eventnames{j}));
    end

    % Target Window
    % Compute erf & cov matrix on the combined data (target window in all 4 conditions appended together)
    cfg = [];
    target_combined = ft_appenddata(cfg, trials_clean.targetchstay, trials_clean.targetchswitch, trials_clean.targetenstay, trials_clean.targetenswitch);
    cfg.covariance       = 'yes';
    cfg.covariancewindow = [0 0.5]; % do not include any period after vocal response onset
    erf_target_combined = ft_timelockanalysis(cfg, target_combined);
    % Compute erf & cov matrix for each condition 
    for j = conds_target
        erf_clean.(eventnames{j}) = ft_timelockanalysis(cfg, trials_clean.(eventnames{j}));
    end
    
    % Compute erf for 'Response' events (just for completeness' sake)
    cfg = [];
    erf_clean.(eventnames{end}) = ft_timelockanalysis(cfg, trials_clean.(eventnames{end}));

end