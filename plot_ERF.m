% Plots ERF for each condition. 
% Generates and plots GFP.

% can use this to regen all plots from saved erf results:
%load([ResultsFolder SubjectID '_erf.mat']); % select which subject to load
%load([ResultsFolder 'lay.mat']); 

% if you want to compare before & after artefact removal (erf vs. erf_clean),
% set plot_uncleaned to 'true' & supply the 'erf' input arg;
% if only plotting cleaned erf & GFP, can specify the 'erf' arg as []

% if you only want the GFP, set both plot_uncleaned & plot_multiplot to false

function plot_ERF (erf, erf_clean, lay, plot_uncleaned, plot_multiplot)

    % run the #define section
    global conds_cue; global conds_target; global eventnames;
    common();


    % required configurations b4 any calls to ft_multiplotER()
    cfg              = [];
    cfg.showlabels   = 'yes';
    cfg.fontsize     = 6;
    cfg.layout       = lay;
    cfg.baseline     = [-0.2 0]; % makes no diff if we've already done baseline correction earlier
    cfg.baselinetype = 'absolute';

    % all pairwise comparisons btwn raw data (erf) & after artefact removal (erf_clean)
    % just to see how good the artefact removal is (atm we reject components 1:5, we can adjust this for more/less removal)
    if (plot_uncleaned)
        for j = 1:length(eventcodes)
            figure;
            ft_multiplotER(cfg, erf.(eventnames{j}), erf_clean.(eventnames{j}));
        end
    end


    % to compare the 4 conds
    % (requires the list of 'cfg' assignments above, if running in console)
    if (plot_multiplot)
        figure('Name','ft_multiplotER: erf_clean.cuechstay, erf_clean.cuechswitch, erf_clean.cueenstay, erf_clean.cueenswitch');
        cfg.xlim = [-0.2 0.75];
        ft_multiplotER(cfg, erf_clean.cuechstay, erf_clean.cuechswitch, erf_clean.cueenstay, erf_clean.cueenswitch);
        legend(eventnames(conds_cue));

        figure('Name','ft_multiplotER: erf_clean.targetchstay, erf_clean.targetchswitch, erf_clean.targetenstay, erf_clean.targetenswitch');
        cfg.xlim = [-0.2 0.75];
        ft_multiplotER(cfg, erf_clean.targetchstay, erf_clean.targetchswitch, erf_clean.targetenstay, erf_clean.targetenswitch);
        legend(eventnames(conds_target));
    end


    %% Calc global averages across all sensors (GFP = global field power)
    cfg        = [];
    cfg.method = 'power';
    %cfg.channel = alldata.label([81:89]); % trigger spike is present in all of these channels: 81-89
    for j = 1:length(eventnames)
        erf_clean_GFP.(eventnames{j}) = ft_globalmeanfield(cfg, erf_clean.(eventnames{j}));
    end

    % plot GFP for cue-locked 
    figure('Name','GFP_cue'); hold on
    %for j = 1:length(eventcodes)
    for j = conds_cue
        %{
        % epoch was [-1 1], we only want to plot [-0.2 0.75]
        % so here we customise the time interval to plot
        start_time = find(erf_clean_GFP.(eventnames{j}).time >= -0.2);
        start_time = start_time(1);
        end_time = find(erf_clean_GFP.(eventnames{j}).time <= 0.75);
        end_time = end_time(end);
        %}
        plot(erf_clean_GFP.(eventnames{j}).time, erf_clean_GFP.(eventnames{j}).avg);
        xlim([-0.2 0.75]); % epoch was [-1 1], we only want to plot [-0.2 0.75]
    end
    legend(eventnames(conds_cue));

    % plot GFP for target-locked 
    figure('Name','GFP_target'); hold on
    for j = conds_target
        plot(erf_clean_GFP.(eventnames{j}).time, erf_clean_GFP.(eventnames{j}).avg);
        xlim([-0.2 0.75]); % epoch was [-1 1], we only want to plot [-0.2 0.75]
    end
    legend(eventnames(conds_target));
end
