%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% stats_ROI.m
%
% grand average & statistical analysis on reconstructed ROI activities
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% run the #define section
global conds_cue; global conds_target; global eventnames;
global ResultsFolder_ROI; % all subjects' ROI data are stored here
common();

% remove the 'response' event type, leaving us with 8 actual event types
eventnames_8 = eventnames(1:8);


%% Read data

% find all .mat files in ResultsFolder_ROI
files = dir([ResultsFolder_ROI '*_ROI.mat']);

% each cycle reads in one '.mat' file (ie. one subject's ROI results)
for i = 1:length(files)
    filename = [ResultsFolder_ROI files(i).name];
    load(filename);
    allSubjects_ROIs_bySubjects(i) = ROI_activity;
end

% get a list of all the ROI labels
ROIs_label = fieldnames(allSubjects_ROIs_bySubjects(1));

% reformat allSubjects_ROIs: Subject|ROI|condition -> ROI|condition|Subjects
for k = 1:length(ROIs_label)
    ROI_name = ROIs_label{k};
    allSubjects_ROIs.(ROI_name) = allSubjects_reformat(allSubjects_ROIs_bySubjects, ROI_name, eventnames_8);
%{
    for j = 1:length(eventnames_8) % 4 conditions in cue & 4 conditions in target (total 8)
       % if exist('allSubjects_ROI') % if var already exists, append to it
        %if ~isempty(['allSubjects_ROI.' ROIs_label{k} '.' (eventnames_8{j})]) % if var already exists, append to it
        if isfield(allSubjects_ROI, ROIs_label{k}) % if var already exists, append to it
            allSubjects_ROI.(ROIs_label{k}).(eventnames_8{j}) = [allSubjects_ROI.(ROIs_label{k}).(eventnames_8{j}) ROI_activity.(ROIs_label{k}).(eventnames_8{j})];
        else % if first time, simply assign to the first item
            allSubjects_ROI.(ROIs_label{k}).(eventnames_8{j}) = ROI_activity.(ROIs_label{k}).(eventnames_8{j});
        end
    end
%}
end


%% Grand average across subjects

fprintf('\n= COMPUTING & PLOTTING CROSS-SUBJECT AVERAGES =\n');

% each cycle processes one ROI
for k = 1:length(ROIs_label)
    ROI_name = ROIs_label{k};
    fprintf(['\nROI: ' ROI_name '\n']);

    % compute grand average (across subjects) for each condition
    cfg = [];
    cfg.channel   = {'all'}; % there is only one channel (i.e. the virtual sensor for this ROI)
    cfg.latency   = 'all';
    cfg.parameter = 'avg';
    for j = 1:length(eventnames_8)
        GA_erf.(eventnames_8{j}) = ft_timelockgrandaverage(cfg, allSubjects_ROIs.(ROI_name).(eventnames_8{j}){:});  
        % "{:}" means to use data from all elements of the variable
    end

    GA.(ROI_name) = GA_erf; % store it in the correct field
    
    % Plot the GAs
    
    % cue-locked
    figure('Name', ['Cue window: GA in ' ROI_name]); hold on
    for j = conds_cue
        plot(GA_erf.(eventnames_8{j}).time, GA_erf.(eventnames_8{j}).avg);
    end
    legend(eventnames_8(conds_cue));

    % target-locked 
    figure('Name', ['Target window: GA in ' ROI_name]); hold on
    for j = conds_target
        plot(GA_erf.(eventnames_8{j}).time, GA_erf.(eventnames_8{j}).avg);
    end
    legend(eventnames_8(conds_target));
end

save([ResultsFolder_ROI 'GA.mat'], 'GA');


%% Statistical analysis (to identify time interval of each effect, i.e. temporal clusters)

fprintf('\n= STATS: CLUSTER-BASED PERMUTATION TESTS =\n');

% each cycle processes one ROI
for k = 1:length(ROIs_label)
    
    ROI_name = ROIs_label{k};
    fprintf(['\nROI: ' ROI_name '\n']);

    data = allSubjects_ROIs.(ROI_name); % data for the current ROI
    
    % set some config for the statistical test
    cfg = [];
    cfg.channel = {'all'}; % there is only one channel (i.e. the virtual sensor for this ROI)
    cfg.avgoverchan = 'yes'; % this is necessary (or else FT will ask for cfg.neighbours)
    
    cfg.latency = [0 0.75]; % time interval over which the experimental 
                         % conditions must be compared (in seconds)
    %latency_cue = [0.385 0.585];%[0.4 0.6];%[0.425 0.55]; % time window for cue-locked effect
    %latency_target = [0.22 0.32];%[0.2 0.3];%[0.25 0.3]; % time window for target-locked effect 
                                % tried [0.2 0.4], not sig
    cfg.avgovertime = 'no'; % if yes, this will average over the entire time window chosen in cfg.latency 
                            % (useful when you want to look at a particular component, e.g. to look at M100,
                            % cfg.latency = [0.08 0.12]; cfg.avgovertime = 'yes'; )

    %load([ResultsFolder_ROI 'neighbours.mat']); % this is the sensor layout - it's the same for all subjects (even same across experiments). So just prepare once & save, then load here
    %cfg.neighbours = neighbours;  % same as defined for the between-trials experiment

    cfg.method = 'montecarlo';
    cfg.statistic = 'depsamplesT'; %cfg.statistic = 'ft_statfun_indepsamplesT'; OR 'ft_statfun_depsamplesFmultivariate';
    cfg.correctm = 'cluster';
    cfg.clusteralpha = 0.05;
    cfg.clusterstatistic = 'maxsum';
    %cfg.minnbchan = 3; % minimum number of neighbourhood channels required to be significant 
                       % in order to form a cluster 
                       % (default: 0, ie. each single channel can be considered a cluster).
                       % 4 or 5 is a good choice; 2 is too few coz it's even below
                       % the resolution of the sensor layout(??)

    cfg.tail = 0;
    cfg.clustertail = 0; % 2 tailed test
    cfg.alpha = 0.05;
    cfg.correcttail = 'prob'; % correct for 2-tailedness
    cfg.numrandomization = 500; % Rule of thumb: use 500, and double this number if it turns out 
        % that the p-value differs from the critical alpha-level (0.05 or 0.01) by less than 0.02

    numSubjects = length(files);
    within_design_2x2 = zeros(2, 2*numSubjects);
    within_design_2x2(1, :) = repmat(1:numSubjects, 1, 2);
    within_design_2x2(2, 1:numSubjects) = 1;
    within_design_2x2(2, numSubjects+1:2*numSubjects) = 2;

    cfg.design = within_design_2x2;
    cfg.uvar  = 1; % row of design matrix that contains unit variable (in this case: subjects)
    cfg.ivar  = 2; % row of design matrix that contains independent variable (i.e. the conditions)

    % Run the statistical tests
    
    % Interaction (i.e. calc sw$ in each lang, then submit the 2 sw$ for comparison)
    fprintf('\nCUE window -> Testing lang x ttype interaction:\n');
    [timelock1, timelock2] = combine_conds_for_T_test('interaction', data.cuechstay, data.cuechswitch, data.cueenstay, data.cueenswitch);
    %cfg.latency = latency_cue; % time interval over which the experimental 
    [cue_interaction.(ROI_name)] = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:});
    fprintf('\nTARGET window -> Testing lang x ttype interaction:\n');
    [timelock1, timelock2] = combine_conds_for_T_test('interaction', data.targetchstay, data.targetchswitch, data.targetenstay, data.targetenswitch); %'2-1 vs 4-3');
    %cfg.latency = latency_target; % time interval over which the experimental 
    [target_interaction.(ROI_name)] = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:});
    
    % Main effect of lang (collapse across stay-switch)
    fprintf('\nCUE window -> Main effect of lang:\n');
    [timelock1, timelock2] = combine_conds_for_T_test('main_12vs34', data.cuechstay, data.cuechswitch, data.cueenstay, data.cueenswitch);
    %cfg.latency = latency_cue; % time interval over which the experimental 
    [cue_lang.(ROI_name)] = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:});
    fprintf('\nTARGET window -> Main effect of lang:\n');
    [timelock1, timelock2] = combine_conds_for_T_test('main_12vs34', data.targetchstay, data.targetchswitch, data.targetenstay, data.targetenswitch); %'2-1 vs 4-3');
    %cfg.latency = latency_target; % time interval over which the experimental 
    [target_lang.(ROI_name)] = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:});

    % Main effect of switch (collapse across langs)
    fprintf('\nCUE window -> Main effect of ttype:\n');
    [timelock1, timelock2] = combine_conds_for_T_test('main_13vs24', data.cuechstay, data.cuechswitch, data.cueenstay, data.cueenswitch);
    %cfg.latency = latency_cue; % time interval over which the experimental 
    [cue_ttype.(ROI_name)] = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:});
    fprintf('\nTARGET window -> Main effect of ttype:\n');
    [timelock1, timelock2] = combine_conds_for_T_test('main_13vs24', data.targetchstay, data.targetchswitch, data.targetenstay, data.targetenswitch); %'2-1 vs 4-3');
    %cfg.latency = latency_target; % time interval over which the experimental 
    [target_ttype.(ROI_name)] = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:});

end

save([ResultsFolder_ROI 'stats.mat'], 'cue_interaction', 'cue_lang', 'cue_ttype', 'target_interaction', 'target_lang', 'target_ttype');


%% Find the effects & plot them

% Automatically check all the stats output & read out the time interval
% of each effect (from the stat.mask field)

stats = load([ResultsFolder_ROI 'stats.mat']);
load([ResultsFolder_ROI 'GA.mat']);

% loop thru all 6 stats output (cue/target lang/ttype/interxn) and loop thru all ROIs in each,
% check whether the .mask field has any non-zero entries 
% (these are the effects & they are already cluster-corrected, so doesn't need to be consecutive 1s)
stats_names = fieldnames(stats);
for i = 1:length(stats_names) % each cycle handles one effect (e.g. cue_lang)
    stat_name = stats_names{i};
    ROIs_names = fieldnames(stats.(stat_name)); % get the list of ROI names
    for k = 1:length(ROIs_names) % each cycle handles one ROI
        ROI_name = ROIs_names{k};
        % if the .mask contains any non-zero entries, that's an effect
        effect = find(stats.(stat_name).(ROI_name).mask); 
        if ~isempty(effect) % if there is an effect, we print it out
            %time_points = sprintf(' %d', effect);
            %fprintf('%s has an effect in %s, at these time points:%s.\n', ROI_name, stat_name, time_points);            
            start_time = stats.(stat_name).(ROI_name).time(effect(1));
            end_time = stats.(stat_name).(ROI_name).time(effect(end));
            fprintf('\n%s has an effect in %s, between %.f~%.f ms.\n', ROI_name, stat_name, start_time*1000, end_time*1000); % convert units to ms

            % plot the effect period, overlaid onto the GA plot for this ROI
            if strcmp(stat_name(1:3), 'cue') % this effect occurs in cue window
                figure('Name', [stat_name ' in ' ROI_name]); hold on
                for j = conds_cue
                    plot(GA.(ROI_name).(eventnames_8{j}).time, GA.(ROI_name).(eventnames_8{j}).avg);
                    xlim([-0.2 0.75]); 
                end
                line([start_time start_time], ylim, 'Color','black'); % plot a vertical line at start_time
                line([end_time end_time], ylim, 'Color','black'); % plot a vertical line at end_time
                % make a colour patch for the time interval of the effect
                % (this keeps occupying the front layer, blocking the GA plot)
                %x = [start_time end_time end_time start_time]; % shade between 2 values on x-axis
                %y = [min(ylim)*[1 1] max(ylim)*[1 1]]; % fill up throughout y-axis
                %patch(x,y,'white'); % choose colour
                legend(eventnames_8(conds_cue));
            elseif strcmp(stat_name(1:6), 'target') % this effect occurs in target window
                figure('Name', [stat_name ' in ' ROI_name]); hold on
                for j = conds_target
                    plot(GA.(ROI_name).(eventnames_8{j}).time, GA.(ROI_name).(eventnames_8{j}).avg);
                    xlim([-0.2 0.75]); 
                end
                line([start_time start_time], ylim, 'Color','black'); % plot a vertical line at start_time
                line([end_time end_time], ylim, 'Color','black'); % plot a vertical line at end_time
                legend(eventnames_8(conds_target));                
            else % should never be here
                fprintf('Error: an effect is found, but its not in either cue nor target window.\n');
            end
        else % output a msg even if there's no effect, just so we know the script ran correctly
            %fprintf('%s: No effect in %s\n', stat_name, ROI_name);
        end
    end
end


%% TODO

% Rename the "results_ROI" folder,
% then re-run source_v1 using Centroid method (Robert) & run this script again


%%
% below is old stuff from stats_ERF, irrelevant here
%{
%% Plotting: use ft_clusterplot & ft_topoplot

load([ResultsFolder_ROI 'stats.mat']);
load([ResultsFolder_ROI 'lay.mat']);
load([ResultsFolder_ROI 'GA.mat']); % only required if using ft_topoplot

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
%}