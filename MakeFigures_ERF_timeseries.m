% run the #define section to obtain values for global vars
global ResultsFolder; 
common();

load([ResultsFolder 'GA_erf_allConditions.mat']);


%% Effect 1: cue_ttype_425-550ms_p=0.01
start_time = 0.425;
end_time = 0.550;

% collapse into 2 conds (stay vs switch)
st = GA_erf.cuechstay;
st.avg = (GA_erf.cuechstay.avg + GA_erf.cueenstay.avg) / 2;
sw = GA_erf.cuechswitch;
sw.avg = (GA_erf.cuechswitch.avg + GA_erf.cueenswitch.avg) / 2;

% baseline correction
cfg = [];
cfg.baseline = [-0.1 0];
st = ft_timelockbaseline(cfg, st); 
sw = ft_timelockbaseline(cfg, sw); 

figure('Name', 'Average ERF of significant channels: cue_ttype_425-550ms'); hold on;
cfg        = [];
cfg.channel = {'AG017', 'AG018', 'AG019', 'AG022', 'AG023', 'AG025', 'AG029', 'AG063', 'AG064', 'AG143'}; % 10 sig channels in cluster
ft_singleplotER(cfg, st, sw);

%plot(st.time, st.avg, 'LineWidth',3);
%plot(sw.time, sw.avg, 'LineWidth',3);
xlim([-0.1 0.75]);
xlabel('seconds');
set(gca, 'LineWidth',1.5, 'FontSize',22); % set axes properties
legend({'Stay', 'Switch'}, 'Location','northwest', 'FontSize',22);
box on; % draw a border around the figure

% create shaded region indicating effect duration
ylimits = ylim; ylow = ylimits(1); yhigh = ylimits(2);
x = [start_time end_time end_time start_time]; % specify x,y coordinates of the 4 corners
y = [ylow ylow yhigh yhigh];
patch(x,y,'black', 'FaceAlpha',0.15) % draw the shade (FaceAlpha is transparency)
ylim(ylimits); % ensure ylim doesn't get expanded
hold off;
