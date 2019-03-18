function MakeFigures_ROI_timeseries

    % run the #define section to obtain values for global vars
    global ResultsFolder_ROI; 
    global PLOT_SHADE; global colours; % for plotting shaded boundary on each time course
    common();

    temp = load([ResultsFolder_ROI 'GA.mat']);
    GA = temp.GA;
    temp = load([ResultsFolder_ROI 'GA_individuals.mat']);
    GA_indi = temp.GA_indi;


    %% Effect 1: cue_interaction_LIFG_315-345ms, p = 0.0300
    % [L2 positive sw$, L1 negative sw$; neither was sig in posthoc t-test]
    ROI_name = 'LIFG';
    start_time = 0.315;
    end_time = 0.345;

    
    % collapse into 2 conditions, based on the cross-subject grand avg
    % (these will be the 2 lines in the plot, each line shows the time course for 1 collapsed cond)
    en_switchcost = GA.(ROI_name).cueenstay;
    en_switchcost.avg = GA.(ROI_name).cueenswitch.avg - GA.(ROI_name).cueenstay.avg;
    ch_switchcost = GA.(ROI_name).cuechstay;
    ch_switchcost.avg = GA.(ROI_name).cuechswitch.avg - GA.(ROI_name).cuechstay.avg;

    % baseline correction
    cfg = [];
    cfg.baseline = [-0.1 0];
    en_switchcost = ft_timelockbaseline(cfg, en_switchcost); 
    ch_switchcost = ft_timelockbaseline(cfg, ch_switchcost); 


    % if plotting shaded boundary, need to do the following:
    % collapse into 2 conditions, based on each individual subject's ROI time course
    % (these will be used to calc the shaded boundary around each line in the plot)
    if ~strcmp(PLOT_SHADE, 'no')
        en_indi = GA_indi.(ROI_name).cueenstay;
        en_indi.individual = GA_indi.(ROI_name).cueenswitch.individual - GA_indi.(ROI_name).cueenstay.individual;
        ch_indi = GA_indi.(ROI_name).cuechstay;
        ch_indi.individual = GA_indi.(ROI_name).cuechswitch.individual + GA_indi.(ROI_name).cuechstay.individual;
    end
    

    % plot
    figure('Name', 'cue_interaction_LIFG_315-345ms'); hold on;

    if strcmp(PLOT_SHADE, 'no') % do not plot shaded boundary, only plot the lines                 
        plot(en_switchcost.time, en_switchcost.avg); % 'LineWidth',3
        plot(ch_switchcost.time, ch_switchcost.avg); % 'LineWidth',3
    else
        margin_en = calc_margin(en_indi.individual, PLOT_SHADE);
        margin_ch = calc_margin(ch_indi.individual, PLOT_SHADE);
        
        % plot time courses with shaded boundary
        boundedline(en_switchcost.time, en_switchcost.avg, margin_en(:), 'alpha', 'transparency',0.15, colours(1));
        boundedline(ch_switchcost.time, ch_switchcost.avg, margin_ch(:), 'alpha', 'transparency',0.15, colours(2));        
    end
    
    xlim([-0.1 0.75]);
    xlabel('Seconds');
    ylabel('Ampere per square metre');
    set(gca, 'LineWidth',1.5, 'FontSize',22); % set axes properties
    box on; % draw a border around the figure

    % create shaded region indicating effect duration
    ylimits = ylim; ylow = ylimits(1); yhigh = ylimits(2);
    x = [start_time end_time end_time start_time]; % specify x,y coordinates of the 4 corners
    y = [ylow ylow yhigh yhigh];
    patch(x,y,'black', 'FaceAlpha',0.15) % draw the shade (FaceAlpha is transparency)
    ylim(ylimits); % ensure ylim doesn't get expanded

    % specify the legend manually (otherwise it will include
    % each shaded patch as an item too). For some reason,
    % the order of the lines are reversed when you grab them
    lines = findall(gcf, 'Type','line');
    legend([lines(2) lines(1)], {'English (L2) switch cost', 'Mandarin (L1) switch cost'}, 'Location','northwest', 'FontSize',30);
    
    set(lines, 'Linewidth',3); % line thickness
    hold off;


    %% Effect 2: target_lang_RIFG_200-235ms_p=0.046
    % [L2 > L1]
    ROI_name = 'RIFG';
    start_time = 0.200;
    end_time = 0.235;


    % collapse into 2 conditions, based on the cross-subject grand avg
    % (these will be the 2 lines in the plot, each line shows the time course for 1 collapsed cond)
    en = GA.(ROI_name).targetenstay;
    en.avg = (GA.(ROI_name).targetenstay.avg + GA.(ROI_name).targetenswitch.avg) / 2;
    ch = GA.(ROI_name).targetchstay;
    ch.avg = (GA.(ROI_name).targetchstay.avg + GA.(ROI_name).targetchswitch.avg) / 2;

    % baseline correction
    cfg = [];
    cfg.baseline = [-0.1 0];
    en = ft_timelockbaseline(cfg, en); 
    ch = ft_timelockbaseline(cfg, ch); 


    % if plotting shaded boundary, need to do the following:
    % collapse into 2 conditions, based on each individual subject's ROI time course
    % (these will be used to calc the shaded boundary around each line in the plot)
    if ~strcmp(PLOT_SHADE, 'no')
        en_indi = GA_indi.(ROI_name).targetenstay;
        en_indi.individual = (GA_indi.(ROI_name).targetenstay.individual + GA_indi.(ROI_name).targetenswitch.individual) / 2;
        ch_indi = GA_indi.(ROI_name).targetchstay;
        ch_indi.individual = (GA_indi.(ROI_name).targetchstay.individual + GA_indi.(ROI_name).targetchswitch.individual) / 2;
    end

    
    % plot
    figure('Name', 'target_lang_RIFG_200-235ms'); hold on;

    if strcmp(PLOT_SHADE, 'no') % do not plot shaded boundary, only plot the lines                 
        plot(en.time, en.avg); % 'LineWidth',3
        plot(ch.time, ch.avg); % 'LineWidth',3
    else
        margin_en = calc_margin(en_indi.individual, PLOT_SHADE);
        margin_ch = calc_margin(ch_indi.individual, PLOT_SHADE);
        
        % plot time courses with shaded boundary
        boundedline(en.time, en.avg, margin_en(:), 'alpha', 'transparency',0.15, colours(1));
        boundedline(ch.time, ch.avg, margin_ch(:), 'alpha', 'transparency',0.15, colours(2));        
    end

    xlim([-0.1 0.75]);
    xlabel('Seconds');
    ylabel('Ampere per square metre');
    set(gca, 'LineWidth',1.5, 'FontSize',22); % set axes properties
    box on; % draw a border around the figure

    % create shaded region indicating effect duration
    ylimits = ylim; ylow = ylimits(1); yhigh = ylimits(2);
    x = [start_time end_time end_time start_time]; % specify x,y coordinates of the 4 corners
    y = [ylow ylow yhigh yhigh];
    patch(x,y,'black', 'FaceAlpha',0.15) % draw the shade (FaceAlpha is transparency)
    ylim(ylimits); % ensure ylim doesn't get expanded

    % specify the legend manually (otherwise it will include
    % each shaded patch as an item too). For some reason,
    % the order of the lines are reversed when you grab them
    lines = findall(gcf, 'Type','line');
    legend([lines(2) lines(1)], {'English (L2)', 'Mandarin (L1)'}, 'Location','northwest', 'FontSize',30);
    
    set(lines, 'Linewidth',3); % line thickness
    hold off;


    %% Marginal effects
    %{
    % cue_lang

    ROI_name = 'RSMA';

    en = GA.(ROI_name).cueenstay;
    en.avg = (GA.(ROI_name).cueenstay.avg + GA.(ROI_name).cueenswitch.avg) / 2;
    ch = GA.(ROI_name).cuechstay;
    ch.avg = (GA.(ROI_name).cuechstay.avg + GA.(ROI_name).cuechswitch.avg) / 2;

    % baseline correction
    cfg = [];
    cfg.baseline = [-0.1 0];
    en = ft_timelockbaseline(cfg, en); 
    ch = ft_timelockbaseline(cfg, ch); 

    figure('Name', 'cue_lang_RSMA_445-465ms_p=0.0500'); hold on;
    plot(en.time, en.avg, 'LineWidth',3);
    plot(ch.time, ch.avg, 'LineWidth',3);
    xlim([-0.1 1]);
    xlabel('Seconds');
    ylabel('Ampere per square metre');
    set(gca, 'LineWidth',1.5, 'FontSize',22); % set axes properties
    box on; % draw a border around the figure

    % create shaded regions indicating effect duration

    % cluster 1
    %{
    start_time = -0.160;
    end_time = -0.135;

    ylimits = ylim; ylow = ylimits(1); yhigh = ylimits(2);
    x = [start_time end_time end_time start_time]; % specify x,y coordinates of the 4 corners
    y = [ylow ylow yhigh yhigh];
    patch(x,y,'black', 'FaceAlpha',0.15) % draw the shade (FaceAlpha is transparency)
    ylim(ylimits); % ensure ylim doesn't get expanded
    %}
    % cluster 2

    start_time = 0.445;
    end_time = 0.465;

    ylimits = ylim; ylow = ylimits(1); yhigh = ylimits(2);
    x = [start_time end_time end_time start_time]; % specify x,y coordinates of the 4 corners
    y = [ylow ylow yhigh yhigh];
    patch(x,y,'black', 'FaceAlpha',0.15) % draw the shade (FaceAlpha is transparency)
    ylim(ylimits); % ensure ylim doesn't get expanded


    legend({'English (L2)', 'Mandarin (L1)'}, 'Location','northeast', 'FontSize',30);
    hold off;


    % cue_ttype

    ROI_name = 'LACC';

    st = GA.(ROI_name).cuechstay;
    st.avg = (GA.(ROI_name).cuechstay.avg + GA.(ROI_name).cueenstay.avg) / 2;
    sw = GA.(ROI_name).cuechswitch;
    sw.avg = (GA.(ROI_name).cuechswitch.avg + GA.(ROI_name).cueenswitch.avg) / 2;

    % baseline correction
    cfg = [];
    cfg.baseline = [-0.1 0];
    st = ft_timelockbaseline(cfg, st); 
    sw = ft_timelockbaseline(cfg, sw); 

    figure('Name', 'cue_ttype_LACC_235-265ms_p=0.0980_and_320-355ms_p=0.0790'); hold on;
    plot(st.time, st.avg, 'LineWidth',3);
    plot(sw.time, sw.avg, 'LineWidth',3);
    xlim([-0.1 1]);
    xlabel('Seconds');
    ylabel('Ampere per square metre');
    set(gca, 'LineWidth',1.5, 'FontSize',22); % set axes properties
    box on; % draw a border around the figure

    % create shaded regions indicating effect duration

    % cluster 1

    start_time = 0.235;
    end_time = 0.265;

    ylimits = ylim; ylow = ylimits(1); yhigh = ylimits(2);
    x = [start_time end_time end_time start_time]; % specify x,y coordinates of the 4 corners
    y = [ylow ylow yhigh yhigh];
    patch(x,y,'black', 'FaceAlpha',0.15) % draw the shade (FaceAlpha is transparency)
    ylim(ylimits); % ensure ylim doesn't get expanded

    % cluster 2

    start_time = 0.320;
    end_time = 0.355;

    ylimits = ylim; ylow = ylimits(1); yhigh = ylimits(2);
    x = [start_time end_time end_time start_time]; % specify x,y coordinates of the 4 corners
    y = [ylow ylow yhigh yhigh];
    patch(x,y,'black', 'FaceAlpha',0.15) % draw the shade (FaceAlpha is transparency)
    ylim(ylimits); % ensure ylim doesn't get expanded


    legend({'Stay', 'Switch'}, 'Location','northeast', 'FontSize',30);
    hold off;
    %}

   
%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SUBFUNCTIONS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end % main function end

    

%% BU for old effects

% not acceptable definition of dlPFC:
%{
%% Effect 1: cue_ttype_LdlPFC_245-275ms_p=0.038
ROI_name = 'LdlPFC';
start_time = 0.245;
end_time = 0.275;

st = GA.(ROI_name).cuechstay;
st.avg = (GA.(ROI_name).cuechstay.avg + GA.(ROI_name).cueenstay.avg) / 2;
sw = GA.(ROI_name).cuechswitch;
sw.avg = (GA.(ROI_name).cuechswitch.avg + GA.(ROI_name).cueenswitch.avg) / 2;

% baseline correction
cfg = [];
cfg.baseline = [-0.1 0];
st = ft_timelockbaseline(cfg, st); 
sw = ft_timelockbaseline(cfg, sw); 

figure('Name', 'cue_ttype_LdlPFC_245-275ms'); hold on;
plot(st.time, st.avg, 'LineWidth',3);
plot(sw.time, sw.avg, 'LineWidth',3);
xlim([-0.1 0.75]);
xlabel('Seconds');
ylabel('Ampere per square metre');
set(gca, 'LineWidth',1.5, 'FontSize',22); % set axes properties
box on; % draw a border around the figure

% create shaded region indicating effect duration
ylimits = ylim; ylow = ylimits(1); yhigh = ylimits(2);
x = [start_time end_time end_time start_time]; % specify x,y coordinates of the 4 corners
y = [ylow ylow yhigh yhigh];
patch(x,y,'black', 'FaceAlpha',0.15) % draw the shade (FaceAlpha is transparency)
ylim(ylimits); % ensure ylim doesn't get expanded

legend({'Stay', 'Switch'}, 'Location','northwest', 'FontSize',30);
hold off;
%}

% incorrect coreg from MEMES1/HCP due to lack of facial points in my hsp:
%{
%% Effect 1: cue_ttype_LdlPFC_420-450ms_p=0.04
ROI_name = 'LdlPFC';
start_time = 0.420;
end_time = 0.450;

st = GA.(ROI_name).cuechstay;
st.avg = (GA.(ROI_name).cuechstay.avg + GA.(ROI_name).cueenstay.avg) / 2;
sw = GA.(ROI_name).cuechswitch;
sw.avg = (GA.(ROI_name).cuechswitch.avg + GA.(ROI_name).cueenswitch.avg) / 2;

% baseline correction
cfg = [];
cfg.baseline = [-0.1 0];
st = ft_timelockbaseline(cfg, st); 
sw = ft_timelockbaseline(cfg, sw); 

figure('Name', 'cue_ttype_LdlPFC_420-450ms'); hold on;
plot(st.time, st.avg, 'LineWidth',3);
plot(sw.time, sw.avg, 'LineWidth',3);
xlim([-0.1 0.75]);
xlabel('seconds');
set(gca, 'LineWidth',1.5, 'FontSize',22); % set axes properties
box on; % draw a border around the figure

% create shaded region indicating effect duration
ylimits = ylim; ylow = ylimits(1); yhigh = ylimits(2);
x = [start_time end_time end_time start_time]; % specify x,y coordinates of the 4 corners
y = [ylow ylow yhigh yhigh];
patch(x,y,'black', 'FaceAlpha',0.15) % draw the shade (FaceAlpha is transparency)
ylim(ylimits); % ensure ylim doesn't get expanded

legend({'Stay', 'Switch'}, 'Location','northwest', 'FontSize',22);
hold off;


%% Effect 2: cue_lang_RACC_705-745ms_p=0.03
ROI_name = 'RACC';
start_time = 0.705;
end_time = 0.745;

en = GA.(ROI_name).cueenstay;
en.avg = (GA.(ROI_name).cueenstay.avg + GA.(ROI_name).cueenswitch.avg) / 2;
ch = GA.(ROI_name).cuechstay;
ch.avg = (GA.(ROI_name).cuechstay.avg + GA.(ROI_name).cuechswitch.avg) / 2;

% baseline correction
cfg = [];
cfg.baseline = [-0.1 0];
en = ft_timelockbaseline(cfg, en); 
ch = ft_timelockbaseline(cfg, ch); 

figure('Name', 'cue_lang_RACC_705-745ms'); hold on;
plot(en.time, en.avg, 'LineWidth',3);
plot(ch.time, ch.avg, 'LineWidth',3);
xlim([-0.1 1]);
xlabel('seconds');
set(gca, 'LineWidth',1.5, 'FontSize',22); % set axes properties
box on; % draw a border around the figure

% create shaded region indicating effect duration
ylimits = ylim; ylow = ylimits(1); yhigh = ylimits(2);
x = [start_time end_time end_time start_time]; % specify x,y coordinates of the 4 corners
y = [ylow ylow yhigh yhigh];
patch(x,y,'black', 'FaceAlpha',0.15) % draw the shade (FaceAlpha is transparency)
ylim(ylimits); % ensure ylim doesn't get expanded

xticks([-0.2 0 0.2 0.4 0.6 0.8 1])
legend({'L2', 'L1'}, 'Location','northwest' ,'FontSize',22);
hold off;


%% Effect 2 (alt plot in target window): target_lang_RACC_-45_to_-5ms
% NOT feasible cos the stats on target window (when incl. the 200ms pre-target) did not detect this effect!
% So we can only report it as part of the cue window


%}