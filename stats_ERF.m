close all;
clear all;

% run the #define section
global conds_cue; global conds_target; global eventnames;
global ResultsFolder; % all subjects' erf data are stored here
common();

% remove the 'response' event type, leaving us with 8 actual event types
eventnames_8 = eventnames(1:8);


% initialise allSubjects_erf (each field holds all subjects' erf in that condition)
allSubjects_erf.cuechstay = {};
allSubjects_erf.cuechswitch = {};
allSubjects_erf.cueenstay = {};
allSubjects_erf.cueenswitch = {};

allSubjects_erf.targetchstay = {};
allSubjects_erf.targetchswitch = {};
allSubjects_erf.targetenstay = {};
allSubjects_erf.targetenswitch = {};



%% Read data

% find all .mat files in ResultsFolder
files = dir([ResultsFolder '*_erf.mat']);

% each cycle reads in one '.mat' file (ie. one subject's erf data)
for i = 1:length(files)
    filename = [ResultsFolder files(i).name];
    load(filename);
        
    for j = 1:length(eventnames_8) % 4 conditions in cue & 4 conditions in target (total 8)
        allSubjects_erf.(eventnames_8{j}) = [allSubjects_erf.(eventnames_8{j}) erf_clean.(eventnames_8{j})];
    end
end


%% Descriptives
% http://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock#within-subjects_experiments

fprintf('\n= COMPUTING & PLOTTING CROSS-SUBJECT AVERAGES =\n');

% CALCULATE the grand average (across all subjects) for each condition
cfg = [];
cfg.channel   = {'all', '-AG083', '-AG087', '-AG088', '-AG082', '-AG084', '-AG086'}; % remove sensors suffering from high-freq noise & trigger leaks
cfg.latency   = 'all';
cfg.parameter = 'avg';
for j = 1:length(eventnames_8)
    GA_erf.(eventnames_8{j}) = ft_timelockgrandaverage(cfg, allSubjects_erf.(eventnames_8{j}){:});  
    % "{:}" means to use data from all elements of the variable
end

%save([ResultsFolder 'GA_erf_allConditions.mat'], 'GA_erf');

% plot
cfg              = [];
cfg.showlabels   = 'yes';
cfg.fontsize     = 6;
cfg.layout       = lay;
cfg.baseline     = [-0.2 0];
cfg.baselinetype = 'absolute';

% cue
figure('Name','ft_multiplotER: GA_erf.cuechstay, GA_erf.cuechsw, GA_erf.cueenstay, GA_erf.cueensw');
ft_multiplotER(cfg, GA_erf.cuechstay, GA_erf.cuechswitch, GA_erf.cueenstay, GA_erf.cueenswitch);
legend(eventnames_8(conds_cue));

% target
figure('Name','ft_multiplotER: GA_erf.targetchstay, GA_erf.targetchsw, GA_erf.targetenstay, GA_erf.targetensw');
ft_multiplotER(cfg, GA_erf.targetchstay, GA_erf.targetchswitch, GA_erf.targetenstay, GA_erf.targetenswitch);
legend(eventnames_8(conds_target));



% CALCULATE global averages across all sensors (i.e. GFP = global field potentials)
cfg        = [];
cfg.method = 'power';
%cfg.channel = {'AG017', 'AG018', 'AG019', 'AG022', 'AG023', 'AG025', 'AG029', 'AG063', 'AG064', 'AG143'}; % 10 sig channels in cluster
for j = 1:length(eventnames_8)
    GA_erf_GFP.(eventnames_8{j}) = ft_globalmeanfield(cfg, GA_erf.(eventnames_8{j}));
end

% plot GFP for cue-locked 
figure('Name','GFP_cue'); hold on
for j = conds_cue
    plot(GA_erf_GFP.(eventnames_8{j}).time, GA_erf_GFP.(eventnames_8{j}).avg);
end
legend(eventnames_8(conds_cue));

% plot GFP for target-locked 
figure('Name','GFP_target'); hold on
for j = conds_target
    plot(GA_erf_GFP.(eventnames_8{j}).time, GA_erf_GFP.(eventnames_8{j}).avg);
end
legend(eventnames_8(conds_target));


% average across all 4 conds (for selecting windows for peaks)
averageAcrossConds_cue = GA_erf_GFP.cuechstay;
averageAcrossConds_cue.avg = (GA_erf_GFP.cuechstay.avg + GA_erf_GFP.cuechswitch.avg + GA_erf_GFP.cueenstay.avg + GA_erf_GFP.cueenswitch.avg) / 4;
averageAcrossConds_target = GA_erf_GFP.targetchstay;
averageAcrossConds_target.avg = (GA_erf_GFP.targetchstay.avg + GA_erf_GFP.targetchswitch.avg + GA_erf_GFP.targetenstay.avg + GA_erf_GFP.targetenswitch.avg) / 4;

figure('Name','Average across all conds - cue window'); 
plot(averageAcrossConds_cue.time, averageAcrossConds_cue.avg); xlim([-0.2 0.75]);
figure('Name','Average across all conds - target window'); 
plot(averageAcrossConds_target.time, averageAcrossConds_target.avg); xlim([-0.2 0.75]);


%% Statistical analysis

fprintf('\n= STATS: CLUSTER-BASED PERMUTATION TESTS =\n');

cfg = [];
cfg.channel = {'all', '-AG083', '-AG087', '-AG088', '-AG082', '-AG084', '-AG086'}; % {'MEG'};
%cfg.latency = [0 0.75]; % time interval over which the experimental 
                     % conditions must be compared (in seconds)
latency_cue = [0.385 0.585];%[0.4 0.6];%[0.425 0.55]; % time window for cue-locked effect
latency_target = [0.2 0.3];%[0.22 0.32];%[0.25 0.3]; % time window for target-locked effect 
                            % tried [0.2 0.4], not sig

load([ResultsFolder 'neighbours.mat']); % this is the sensor layout - it's the same for all subjects (even same across experiments). So just prepare once & save, then load here
cfg.neighbours = neighbours;  % same as defined for the between-trials experiment

cfg.avgovertime = 'yes'; % if yes, this will average over the entire time window chosen in cfg.latency 
                        % (useful when you want to look at a particular component, e.g. to look at M100,
                        % cfg.latency = [0.08 0.12]; cfg.avgovertime = 'yes'; )
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT'; %cfg.statistic = 'ft_statfun_indepsamplesT'; OR 'ft_statfun_depsamplesFmultivariate';
cfg.correctm = 'cluster'; %'no'; % its common in MEG studies to run uncorrected at alpha = 0.001
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 3; % minimum number of neighbourhood channels required to be significant 
                   % in order to form a cluster 
                   % (default: 0, ie. each single channel can be considered a cluster).
                   % 4 or 5 is a good choice; 2 is too few coz it's even below
                   % the resolution of the sensor layout(??)

cfg.tail = 0;
cfg.clustertail = 0; % 2 tailed test
cfg.alpha = 0.05;
cfg.correcttail = 'prob'; % correct for 2-tailedness
cfg.numrandomization = 1000; % Rule of thumb: use 500, and double this number if it turns out 
    % that the p-value differs from the critical alpha-level (0.05 or 0.01) by less than 0.02

numSubjects = length(files);
within_design_2x2 = zeros(2,2*numSubjects);
within_design_2x2(1,:) = repmat(1:numSubjects,1,2);
within_design_2x2(2,1:numSubjects)        = 1;
within_design_2x2(2,numSubjects+1:2*numSubjects) = 2;

cfg.design = within_design_2x2;
cfg.uvar  = 1; % row of design matrix that contains unit variable (in this case: subjects)
cfg.ivar  = 2; % row of design matrix that contains independent variable (i.e. the conditions)

% Interaction (i.e. calc sw$ in each lang, then test the 2 sw$)
% http://www.fieldtriptoolbox.org/faq/how_can_i_test_an_interaction_effect_using_cluster-based_permutation_tests
for i = 1:numSubjects
    allSubj_cue_ch_switchCost{i} = allSubjects_erf.cuechstay{i}; % difference btwn chstay & chswitch (calc'd on next line)
    allSubj_cue_ch_switchCost{i}.avg = allSubjects_erf.cuechswitch{i}.avg - allSubjects_erf.cuechstay{i}.avg; % chswitch - chstay
    allSubj_cue_en_switchCost{i} = allSubjects_erf.cueenstay{i};
    allSubj_cue_en_switchCost{i}.avg = allSubjects_erf.cueenswitch{i}.avg - allSubjects_erf.cueenstay{i}.avg; % enswitch - enstay
    allSubj_target_ch_switchCost{i} = allSubjects_erf.targetchstay{i};
    allSubj_target_ch_switchCost{i}.avg = allSubjects_erf.targetchswitch{i}.avg - allSubjects_erf.targetchstay{i}.avg; % chswitch - chstay
    allSubj_target_en_switchCost{i} = allSubjects_erf.targetenstay{i};
    allSubj_target_en_switchCost{i}.avg = allSubjects_erf.targetenswitch{i}.avg - allSubjects_erf.targetenstay{i}.avg; % enswitch - enstay
end
fprintf('\nCUE window -> Testing lang x ttype interaction:\n');
cfg.latency = latency_cue; % time interval over which the experimental 
[cue_interaction] = ft_timelockstatistics(cfg, allSubj_cue_ch_switchCost{:}, allSubj_cue_en_switchCost{:});
fprintf('\nTARGET window -> Testing lang x ttype interaction:\n');
cfg.latency = latency_target; % time interval over which the experimental 
[target_interaction] = ft_timelockstatistics(cfg, allSubj_target_ch_switchCost{:}, allSubj_target_en_switchCost{:});

% Main effect of lang (collapse across stay-switch)
for i = 1:numSubjects
    allSubj_cue_ch{i} = allSubjects_erf.cuechstay{i};
    allSubj_cue_ch{i}.avg = (allSubjects_erf.cuechstay{i}.avg + allSubjects_erf.cuechswitch{i}.avg) / 2; % cue_ch_all.avg = (cuechstay.avg + cuechsw.avg) / 2
    allSubj_cue_en{i} = allSubjects_erf.cueenstay{i};
    allSubj_cue_en{i}.avg = (allSubjects_erf.cueenstay{i}.avg + allSubjects_erf.cueenswitch{i}.avg) / 2;
    allSubj_target_ch{i} = allSubjects_erf.targetchstay{i};
    allSubj_target_ch{i}.avg = (allSubjects_erf.targetchstay{i}.avg + allSubjects_erf.targetchswitch{i}.avg) / 2; % target_ch_all.avg = (targetchstay.avg + targetchsw.avg) / 2
    allSubj_target_en{i} = allSubjects_erf.targetenstay{i};
    allSubj_target_en{i}.avg = (allSubjects_erf.targetenstay{i}.avg + allSubjects_erf.targetenswitch{i}.avg) / 2;    
end
fprintf('\nCUE window -> Main effect of lang:\n');
cfg.latency = latency_cue; % time interval over which the experimental 
[cue_lang] = ft_timelockstatistics(cfg, allSubj_cue_ch{:}, allSubj_cue_en{:});
fprintf('\nTARGET window -> Main effect of lang:\n');
cfg.latency = latency_target; % time interval over which the experimental 
[target_lang] = ft_timelockstatistics(cfg, allSubj_target_ch{:}, allSubj_target_en{:});

% Main effect of switch (collapse across langs)
for i = 1:numSubjects
    allSubj_cue_stay{i} = allSubjects_erf.cuechstay{i};
    allSubj_cue_stay{i}.avg = (allSubjects_erf.cuechstay{i}.avg + allSubjects_erf.cueenstay{i}.avg) / 2;
    allSubj_cue_switch{i} = allSubjects_erf.cuechswitch{i};
    allSubj_cue_switch{i}.avg = (allSubjects_erf.cuechswitch{i}.avg + allSubjects_erf.cueenswitch{i}.avg) / 2;
    allSubj_target_stay{i} = allSubjects_erf.targetchstay{i};
    allSubj_target_stay{i}.avg = (allSubjects_erf.targetchstay{i}.avg + allSubjects_erf.targetenstay{i}.avg) / 2;
    allSubj_target_switch{i} = allSubjects_erf.targetchswitch{i};
    allSubj_target_switch{i}.avg = (allSubjects_erf.targetchswitch{i}.avg + allSubjects_erf.targetenswitch{i}.avg) / 2;
end
fprintf('\nCUE window -> Main effect of ttype:\n');
cfg.latency = latency_cue; % time interval over which the experimental 
[cue_ttype] = ft_timelockstatistics(cfg, allSubj_cue_stay{:}, allSubj_cue_switch{:});
fprintf('\nTARGET window -> Main effect of ttype:\n');
cfg.latency = latency_target; % time interval over which the experimental 
[target_ttype] = ft_timelockstatistics(cfg, allSubj_target_stay{:}, allSubj_target_switch{:});

%save([ResultsFolder 'stats.mat'], 'cue_interaction', 'cue_lang', 'cue_ttype', 'target_interaction', 'target_lang', 'target_ttype');


%% Plotting: use ft_clusterplot & ft_topoplot

load([ResultsFolder 'stats.mat']);
load([ResultsFolder 'lay.mat']);
load([ResultsFolder 'GA_erf_allConditions.mat']); % only required if using ft_topoplot

% select which comparison to plot
stat = cue_ttype; % here we plot the only effect that seems to survive correction (at minnbchan = 0)
                  % to explore where (both in terms of time & location) the effect might have possibly
                  % occurred
                  % [TODO] then we can define more precise time window &
                  % set avgovertime = 'yes', which should give us more
                  % sensitivity, and allow us to increase the minnbchan to
                  % a reasonable number: 2 (ft tutorial) or 4 (Paul)

%% ft_clusterplot (based on t-values)
cfg = [];
%cfg.zlim = [-5 5]; % set scaling (range of t-values) (usually using automatic is ok) 
cfg.highlightcolorpos = [1 1 1]; % white for pos clusters
cfg.highlightcolorneg = [255/255 192/255 203/255]; % pink for neg clusters
cfg.alpha = 0.05;
%cfg.colorbar = 'yes'; % shows the scaling
cfg.layout = lay;
ft_clusterplot(cfg, stat);


%% ft_topoplot (based on actual erf amplitude) 

% first, define the 2 conds to be compared (this time using cross-subject averages, i.e. GA)
% here we look at main effect of ttype in cue window, so we collapse across langs
GA_cue_stay = GA_erf.cuechstay;
GA_cue_stay.avg = (GA_erf.cuechstay.avg + GA_erf.cueenstay.avg) / 2;
GA_cue_switch = GA_erf.cuechswitch;
GA_cue_switch.avg = (GA_erf.cuechswitch.avg + GA_erf.cueenswitch.avg) / 2;

% then, calc the diff btwn the 2 conds
cfg  = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
GA_cue_stayvsswitch = ft_math(cfg, GA_cue_stay, GA_cue_switch);


% define parameters for plotting
start_time = stat.cfg.latency(1); % get the time window specified earlier in stat analysis
end_time = stat.cfg.latency(end);
timestep = 0.05; %(end_time - start_time) / 15; % length of time interval you want in each subplot (in seconds); 
                                                % alt: specify how many subplots you want (e.g. 15)
sampling_rate = 200; % we downsampled to 200Hz
sample_count = length(stat.time); % number of samples in MEG data (in the ERF time window)
j = [start_time : timestep : end_time];   % define the time interval (in seconds) for each subplot
m = [1 : timestep*sampling_rate : sample_count];  % corresponding sample indices in MEG data

% ensure stat.cfg.alpha (the alpha level we specified earlier in ft_timelockstatistics) still exists
if ~isfield(stat.cfg,'alpha'); stat.cfg.alpha = 0.05; end; % if not, set it (new version corrects for 2-tailedness, so no need to use 0.025)

%{
if (length(stat.posclusters) == 0) % if no clusters were found at all, code below will throw error
    % so create a fake one (just to allow code below to run w/o error)
    stat.posclusters(1).prob = 1; 
    stat.posclusters(1).clusterstat = -9; 
    stat.posclusters(1).stddev = 0; 
    stat.posclusters(1).cirange = 0;
end
if (length(stat.negclusters) == 0) % do the same for neg clusters
    stat.negclusters(1).prob = 1; 
    stat.negclusters(1).clusterstat = -9; 
    stat.negclusters(1).stddev = 0; 
    stat.negclusters(1).cirange = 0;
end
%}

% get all p-values associated with the clusters
pos_cluster_pvals = [stat.posclusters(:).prob];
neg_cluster_pvals = [stat.negclusters(:).prob];
% find which clusters are significant, outputting their indices as held in stat.posclusters
pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
neg_signif_clust = find(neg_cluster_pvals < stat.cfg.alpha);
% make a boolean matrix of which (channel,time)-pairs are part of a significant cluster
pos = ismember(stat.posclusterslabelmat, pos_signif_clust); % I think stat.mask is simply combining pos & neg 
neg = ismember(stat.negclusterslabelmat, neg_signif_clust); % (i.e. stat.mask == pos | neg)

% Ensure the channels have the same order in the grand average and in the statistical output
% This might not be the case, because ft_math might shuffle the order  
[i1,i2] = match_str(GA_cue_stayvsswitch.label, stat.label);
% i1 holds a list of channel numbers in the grand averages
% i2 holds a list of channel numbers in the stat output

figure;  
for k = 1:length(j)-1; % create one subplot for each time interval
     subplot(3,5,k); % 3 * 5 = 15 subplots 
     
     cfg = [];   
     cfg.xlim=[j(k) j(k+1)];   
     %cfg.zlim = [-5e-14 5e-14];  % set scaling (usually using automatic is ok) 
     pos_int = zeros(numel(GA_cue_stayvsswitch.label),1); % initialise the arrays with 0s
     neg_int = zeros(numel(GA_cue_stayvsswitch.label),1);
     pos_int(i1) = all(pos(i2, m(k):m(k+1)), 2); % if a channel maintains significance thruout this time interval, then
     neg_int(i1) = all(neg(i2, m(k):m(k+1)), 2); % we set this channel to 1 (to be highlighted)
     % not sure why it has to "maintain significance"; here I try with only requiring sig for half of time pts in this interval
     a = neg(i2, m(k):m(k+1));
     neg_int(i1) = sum(a, 2) > size(a, 2) / 2;
     
     sig_channels = find(pos_int | neg_int); % get indices of all significant channels
     if length(sig_channels) ~= 0 % if any sig channels found, report which channels these are
         fprintf(['In time interval [' num2str(cfg.xlim) '], these channels were significant:\n']);
         stat.label(sig_channels)
     end
     cfg.highlight = 'on';
     cfg.highlightchannel = sig_channels; % highlight these channels on topoplot
     cfg.highlightcolor = [255/255 192/255 203/255]; % pink colour

     cfg.comment = ['time = [' num2str(cfg.xlim) ']   ' strjoin(stat.label(sig_channels))]; % display time interval & names of sig channels
     %cfg.comment = 'auto'; % display date, xlim (time interval), zlim (amplitude range)
     cfg.commentpos = 'title';   
     %cfg.colorbar = 'yes'; % shows the scaling
     cfg.layout = lay;
     ft_topoplotER(cfg, GA_cue_stayvsswitch);
end  

%% To plot the actual effect (i.e. average ERF of sig channels)
% alt: simply go to the multiplotER generated earlier, select the sig channels & plot

% effect in cue window
cfg        = [];
cfg.channel = stat.label(find(cue_ttype.mask)); % autoly retrieve sig channels (only works with cfg.avgovertime = 'yes')
%{'AG017', 'AG018', 'AG019', 'AG022', 'AG023', 'AG025', 'AG029', 'AG063', 'AG064', 'AG143'}; % 10 sig channels in cluster

figure('Name','Average ERF of significant channels - cue window');
ft_singleplotER(cfg, GA_erf.cuechstay, GA_erf.cuechswitch, GA_erf.cueenstay, GA_erf.cueenswitch);
legend(eventnames_8(conds_cue));

% effect in target window
cfg        = [];
cfg.channel = stat.label(find(target_lang.mask)); % autoly retrieve sig channels (only works with cfg.avgovertime = 'yes')

figure('Name','Average ERF of significant channels - target window');
ft_singleplotER(cfg, GA_erf.targetchstay, GA_erf.targetchswitch, GA_erf.targetenstay, GA_erf.targetenswitch);
legend(eventnames_8(conds_target));


% Ref code: copied directly from FT_compare_conditions.m
% see also: http://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock#the_format_of_the_output
%{

neg_cluster_pvals=[];
pos_cluster_pvals=[];

figure;% % plot negative
cfg = [];
cfg.comment = 'no';
cfg.layout = layout;

cfg.xlim=latency;
cfg.highlight = 'off';

subplot(2,4,1)
eval(['ft_topoplotER(cfg, GM_meg_',cond1,');']);
subplot(2,4,2)
eval(['ft_topoplotER(cfg, GM_meg_',cond2,');']);
subplot(2,4,5)
eval(['ft_topoplotER(cfg, GM_meg_',cond1,'_planar);']);
subplot(2,4,6)
eval(['ft_topoplotER(cfg, GM_meg_',cond2,'_planar);']);

if isfield(eval([cond1,'_vs_',cond2,'_stat']),'posclusters') && ~isempty(eval([cond1,'_vs_',cond2,'_stat.posclusters']))
    eval(['pos_cluster_pvals = [',cond1,'_vs_',cond2,'_stat.posclusters(:).prob];']);
    pos_signif_clust = find(pos_cluster_pvals < 0.05);
    eval(['pos = ismember(',cond1,'_vs_',cond2,'_stat.posclusterslabelmat, pos_signif_clust);']);
    eval(['poscluster_p=([',cond1,'_vs_',cond2,'_stat.posclusters.prob])']);
    cfg.highlight = 'on';
    cfg.highlightchannel = find(pos);
    subplot(3,4,4)
    eval(['plot(GM_meg_',cond1,'.time,mean(GM_meg_',cond1,'.avg(logical(',cond1,'_vs_',cond2,'_stat.posclusterslabelmat),:)))'])
    xlim([-0.1 0.5])
    xlabel('Time (s)')
    ylabel('Amplitude (fT)')
    hold on
    eval(['plot(GM_meg_',cond2,'.time,mean(GM_meg_',cond2,'.avg(logical(',cond1,'_vs_',cond2,'_stat.posclusterslabelmat),:)))'])
else
    sprintf('no significant pos clusters')
end
subplot(2,4,3)
eval(['ft_topoplotER(cfg, GM_meg_',cond1,'_vs_GM_meg_',cond2,');']);
cfg.highlight = 'off';
subplot(2,4,7)
eval(['ft_topoplotER(cfg, GM_meg_planar_',cond1,'_vs_GM_meg_planar_',cond2,');']);
set(gcf, 'Color', 'w');
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure. 

figure;% % plot positive
cfg = [];
cfg.comment = 'no';
cfg.layout = layout;

cfg.xlim=latency;
cfg.highlight = 'off';

subplot(2,4,1)
eval(['ft_topoplotER(cfg, GM_meg_',cond1,');']);
subplot(2,4,2)
eval(['ft_topoplotER(cfg, GM_meg_',cond2,');']);
subplot(2,4,5)
eval(['ft_topoplotER(cfg, GM_meg_',cond1,'_planar);']);
subplot(2,4,6)
eval(['ft_topoplotER(cfg, GM_meg_',cond2,'_planar);']);

if isfield(eval([cond1,'_vs_',cond2,'_stat']),'negclusters') && ~isempty(eval([cond1,'_vs_',cond2,'_stat.negclusters']))
    eval(['neg_cluster_pvals = [',cond1,'_vs_',cond2,'_stat.negclusters(:).prob];']);
    neg_signif_clust = find(neg_cluster_pvals < 0.05);
    eval(['neg = ismember(',cond1,'_vs_',cond2,'_stat.negclusterslabelmat, neg_signif_clust);']);
    eval(['negcluster_p=([',cond1,'_vs_',cond2,'_stat.negclusters.prob])']);
    cfg.highlight = 'on';
    cfg.highlightchannel = find(neg);
    subplot(3,4,4)
    eval(['plot(GM_meg_',cond1,'.time,mean(GM_meg_',cond1,'.avg(logical(',cond1,'_vs_',cond2,'_stat.negclusterslabelmat),:)))'])
    xlim([-0.1 0.5])
    xlabel('Time (s)')
    ylabel('Amplitude (fT)')
    hold on
    eval(['plot(GM_meg_',cond2,'.time,mean(GM_meg_',cond2,'.avg(logical(',cond1,'_vs_',cond2,'_stat.negclusterslabelmat),:)))'])
else
    sprintf('no significant neg clusters')
end
subplot(2,4,3)
eval(['ft_topoplotER(cfg, GM_meg_',cond1,'_vs_GM_meg_',cond2,');']);
cfg.highlight = 'off';
subplot(2,4,7)
eval(['ft_topoplotER(cfg, GM_meg_planar_',cond1,'_vs_GM_meg_planar_',cond2,');']);
set(gcf, 'Color', 'w');
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure. 
end

%}
