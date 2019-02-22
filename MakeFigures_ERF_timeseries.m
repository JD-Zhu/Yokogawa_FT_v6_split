% run the #define section to obtain values for global vars
global ResultsFolder; 
common();

% load the results
load([ResultsFolder 'GA_erf_allConditions.mat']);
load([ResultsFolder 'stats.mat']);


%% Target-locked: collapse GA into 2 conds 
% (Eng vs Chn)
en = GA_erf.targetenstay;
en.avg = (GA_erf.targetenstay.avg + GA_erf.targetenswitch.avg) / 2;
ch = GA_erf.targetchstay;
ch.avg = (GA_erf.targetchstay.avg + GA_erf.targetchswitch.avg) / 2;

% (stay vs switch)
st = GA_erf.targetchstay;
st.avg = (GA_erf.targetchstay.avg + GA_erf.targetenstay.avg) / 2;
sw = GA_erf.targetchswitch;
sw.avg = (GA_erf.targetchswitch.avg + GA_erf.targetenswitch.avg) / 2;

% baseline correction
cfg = [];
cfg.baseline = [-0.1 0];
en = ft_timelockbaseline(cfg, en); 
ch = ft_timelockbaseline(cfg, ch);
st = ft_timelockbaseline(cfg, st); 
sw = ft_timelockbaseline(cfg, sw);


%% Effect 1: target_lang, 360-515ms, p = 0.02
stat = target_lang;
start_time = 0.360;
end_time = 0.515;

% produce a list of all channels that were sig 
% at one or more time points during the effect interval
sig_channels = [];

% each cycle checks one channel
for i = 1:size(stat.mask, 1)
    if find(stat.mask(i,:)) % this channel was sig at some point
        sig_channels = [sig_channels i]; % so add it to the list
    end
end

% plot the ERF (averaged over all sig channels)
figure('Name', 'Average ERF of significant channels: target_lang_360-515ms'); hold on;
cfg        = [];
cfg.channel = sig_channels; % only include the sig channels
cfg.linewidth = 3;
ft_singleplotER(cfg, en, ch);

%plot(en.time, en.avg, 'LineWidth',3);
%plot(ch.time, ch.avg, 'LineWidth',3);
xlim([-0.1 0.75]);
xlabel('Seconds');
ylabel('Tesla');
set(gca, 'LineWidth',1.5, 'FontSize',22); % set axes properties
legend({'English (L2)', 'Mandarin (L1)'}, 'Location','northwest', 'FontSize',30);
box on; % draw a border around the figure

% create shaded region indicating effect duration
ylimits = ylim; ylow = ylimits(1); yhigh = ylimits(2);
x = [start_time end_time end_time start_time]; % specify x,y coordinates of the 4 corners
y = [ylow ylow yhigh yhigh];
patch(x,y,'black', 'FaceAlpha',0.15) % draw the shade (FaceAlpha is transparency)
ylim(ylimits); % ensure ylim doesn't get expanded
hold off;


%% Effect 2: target_ttype, 355-465ms, p = 0.07 (marginal)
stat = target_ttype;
start_time = 0.355;
end_time = 0.465;

% Extra step for marginal effects:
% Because the "stats" output we saved had alpha=0.05, so marginal effect 
% is not showing up in the .mask field. Here we recreate the .mask field 
% to show marginal effect.
% This approach for creating the .mask field has been verified to be correct.
idx = find(stat.prob < 0.1); % find indices where p < 0.1
stat.mask(idx) = 1; % fill these positions with 1s


% produce a list of all channels that were sig 
% at one or more time points during the effect interval
sig_channels = [];

% each cycle checks one channel
for i = 1:size(stat.mask, 1)
    if find(stat.mask(i,:)) % this channel was sig at some point
        sig_channels = [sig_channels i]; % so add it to the list
    end
end

% plot the ERF (averaged over all sig channels)
figure('Name', 'Average ERF of significant channels: target_ttype_355-465ms'); hold on;
cfg        = [];
cfg.channel = sig_channels; % only include the sig channels
cfg.linewidth = 3;
ft_singleplotER(cfg, st, sw);

%plot(st.time, st.avg, 'LineWidth',3);
%plot(sw.time, sw.avg, 'LineWidth',3);
xlim([-0.1 0.75]);
xlabel('Seconds');
ylabel('Tesla');
set(gca, 'LineWidth',1.5, 'FontSize',22); % set axes properties
legend({'Stay', 'Switch'}, 'Location','northwest', 'FontSize',30);
box on; % draw a border around the figure

% create shaded region indicating effect duration
ylimits = ylim; ylow = ylimits(1); yhigh = ylimits(2);
x = [start_time end_time end_time start_time]; % specify x,y coordinates of the 4 corners
y = [ylow ylow yhigh yhigh];
patch(x,y,'black', 'FaceAlpha',0.15) % draw the shade (FaceAlpha is transparency)
ylim(ylimits); % ensure ylim doesn't get expanded
hold off;
