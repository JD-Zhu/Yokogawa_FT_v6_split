% run the #define section to obtain values for global vars
global ResultsFolder_ROI; 
common();

load([ResultsFolder_ROI 'GA.mat']);


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
cfg.baseline = [-0.2 0];
st = ft_timelockbaseline(cfg, st); 
sw = ft_timelockbaseline(cfg, sw); 

figure('Name', 'cue_ttype_LdlPFC_420-450ms'); hold on;
plot(st.time, st.avg, 'LineWidth',3);
plot(sw.time, sw.avg, 'LineWidth',3);
xlim([-0.2 0.75]);
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
cfg.baseline = [-0.2 0];
en = ft_timelockbaseline(cfg, en); 
ch = ft_timelockbaseline(cfg, ch); 

figure('Name', 'cue_lang_RACC_705-745ms'); hold on;
plot(en.time, en.avg, 'LineWidth',3);
plot(ch.time, ch.avg, 'LineWidth',3);
xlim([-0.2 0.75]);
xlabel('seconds');
set(gca, 'LineWidth',1.5, 'FontSize',22); % set axes properties
legend({'L2', 'L1'}, 'Location','northwest' ,'FontSize',22);
box on; % draw a border around the figure

% create shaded region indicating effect duration
ylimits = ylim; ylow = ylimits(1); yhigh = ylimits(2);
x = [start_time end_time end_time start_time]; % specify x,y coordinates of the 4 corners
y = [ylow ylow yhigh yhigh];
patch(x,y,'black', 'FaceAlpha',0.15) % draw the shade (FaceAlpha is transparency)
ylim(ylimits); % ensure ylim doesn't get expanded
hold off;


%% Effect 2 (alt plot in target window): target_lang_RACC_-45_to_-5ms
% NOT feasible cos the stats on target window (when incl. the 200ms pre-target) did not detect this effect!
% So we can only report it as part of the cue window
%{
ROI_name = 'RACC';
start_time = -0.045;
end_time = -0.005;

en = GA.(ROI_name).targetenstay;
en.avg = (GA.(ROI_name).targetenstay.avg + GA.(ROI_name).targetenswitch.avg) / 2;
ch = GA.(ROI_name).targetchstay;
ch.avg = (GA.(ROI_name).targetchstay.avg + GA.(ROI_name).targetchswitch.avg) / 2;

figure('Name', 'target_lang_RACC_-45_to_-5ms'); hold on;
plot(en.time, en.avg);
plot(ch.time, ch.avg);
xlim([-0.2 0.75]);
legend({'English  (L2)', 'Mandarin (L1)'}, 'Location','northwest' ,'FontSize',16);
line([start_time start_time], ylim, 'Color','black'); % plot a vertical line at start_time
line([end_time end_time], ylim, 'Color','black'); % plot a vertical line at end_time
hold off;
%}