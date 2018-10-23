% run the #define section to obtain values for global vars
global ResultsFolder_ROI; 
common();

load([ResultsFolder_ROI 'GA.mat']);


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
cfg.baseline = [-0.2 0];
st = ft_timelockbaseline(cfg, st); 
sw = ft_timelockbaseline(cfg, sw); 

figure('Name', 'cue_ttype_LdlPFC_245-275ms'); hold on;
plot(st.time, st.avg, 'LineWidth',3);
plot(sw.time, sw.avg, 'LineWidth',3);
xlim([-0.2 0.75]);
xlabel('seconds');
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


%% Effect 2: target_lang_RIFG_190-235ms_p=0.019
ROI_name = 'RIFG';
start_time = 0.19;
end_time = 0.235;

en = GA.(ROI_name).targetenstay;
en.avg = (GA.(ROI_name).targetenstay.avg + GA.(ROI_name).targetenswitch.avg) / 2;
ch = GA.(ROI_name).targetchstay;
ch.avg = (GA.(ROI_name).targetchstay.avg + GA.(ROI_name).targetchswitch.avg) / 2;

% baseline correction
cfg = [];
cfg.baseline = [-0.2 0];
en = ft_timelockbaseline(cfg, en); 
ch = ft_timelockbaseline(cfg, ch); 

figure('Name', 'target_lang_RIFG_190-235ms'); hold on;
plot(en.time, en.avg, 'LineWidth',3);
plot(ch.time, ch.avg, 'LineWidth',3);
xlim([-0.2 0.75]);
xlabel('seconds');
set(gca, 'LineWidth',1.5, 'FontSize',22); % set axes properties
box on; % draw a border around the figure

% create shaded region indicating effect duration
ylimits = ylim; ylow = ylimits(1); yhigh = ylimits(2);
x = [start_time end_time end_time start_time]; % specify x,y coordinates of the 4 corners
y = [ylow ylow yhigh yhigh];
patch(x,y,'black', 'FaceAlpha',0.15) % draw the shade (FaceAlpha is transparency)
ylim(ylimits); % ensure ylim doesn't get expanded

legend({'English  (L2)', 'Mandarin (L1)'}, 'Location','northwest', 'FontSize',30);
hold off;


%% below are for old effects 
% (incorrect coreg from MEMES1/HCP due to lack of facial points in my hsp)

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
cfg.baseline = [-0.2 0];
st = ft_timelockbaseline(cfg, st); 
sw = ft_timelockbaseline(cfg, sw); 

figure('Name', 'cue_ttype_LdlPFC_420-450ms'); hold on;
plot(st.time, st.avg, 'LineWidth',3);
plot(sw.time, sw.avg, 'LineWidth',3);
xlim([-0.2 0.75]);
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
cfg.baseline = [-0.2 0];
en = ft_timelockbaseline(cfg, en); 
ch = ft_timelockbaseline(cfg, ch); 

figure('Name', 'cue_lang_RACC_705-745ms'); hold on;
plot(en.time, en.avg, 'LineWidth',3);
plot(ch.time, ch.avg, 'LineWidth',3);
xlim([-0.2 1]);
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