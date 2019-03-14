%{
Plot topography showing location of the sensors which formed the cluster.

ft_clusterplot is convenient for early visualisation of significant channel 
locations. It is a wrapper around ft_topoplot, and automatically extracts
the spatial and temporal extent of the cluster from the stat output. It 
produces a series of subplots, spanning the time window in which the effect
was found (1 subplot per sample).

If you do not wish to include all the subplots in your paper, you can produce 
a plot averaged over time points (i.e. this script). You need to make the 
decision on how to include the spatial extent (i.e. electrodes) of your
cluster. You could show all electrodes that were part of the cluster at any
one time point, or you could plot the electrode size as a function of how
often (i.e. at how many time points) the electrode was part of the cluster.

Here we plot the topography averaged over the interval of the effect (425-550ms), 
and highlight the channels which are significant at the "middle" time sample (485ms).

See here for more information:
http://www.fieldtriptoolbox.org/faq/how_not_to_interpret_results_from_a_cluster-based_permutation_test
%}

%% run the #define section to obtain values for global vars
global ResultsFolder; 
common();


%% load the results
load([ResultsFolder 'stats.mat']);
load([ResultsFolder 'lay.mat']);
load([ResultsFolder 'GA_erf_allConditions.mat']);
load([ResultsFolder 'neighbours.mat']);

% load nice colourmap
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colours = colormap(flipud(brewermap(64, 'RdBu')));

% select which effect to plot
stat = target_lang;
start_time = 0.360;
end_time = 0.515;


%% plot topography for axial gradiometers (based on t-values)

cfg                   = [];
cfg.layout            = lay;
cfg.colormap          = colours;
cfg.colorbar          = 'yes'; % shows the scales
%cfg.colorbar         = 'EastOutside';
cfg.zlim              = [-4 4];%'maxabs'; % set the scaling

cfg.xlim = [start_time end_time]; % duration of the effect (as reported by ft_clusterplot)
                          % topography will be averaged over this interval

cfg.highlight         = 'on'; % highlight significant channels
cfg.highlightsymbol = '.';
cfg.highlightsize = 26;
cfg.highlightcolor = 'w';


% which channels to highlight? 

% here we highlight all channels that were sig 
% at one or more time points during the effect interval
sig_channels = [];
% each cycle checks one channel
for i = 1:size(stat.mask, 1)
    if find(stat.mask(i,:)) % this channel was sig at some point
        sig_channels = [sig_channels i]; % so add it to the list
    end
end
cfg.highlightchannel  = sig_channels;

% Alternatively, we highlight the channels that are sig at the middle time point
%{
time_point = (start_time + end_time) / 2; 
time_point = time_point * 1000; % convert unit to ms
%time_point = 440; % set manually if needed
time_index = round(time_point / 5 + 1); % index of the time point you want, round it up/down to an integer
                                        % note: this is only correct if the epoch starts from time 0;
                                        % if it starts from -100ms, then need to adjust index accordingly
time_index = 109; % set manually if needed
cfg.highlightchannel  = stat.label(ismember(stat.negclusterslabelmat(:,time_index),1)); % find all the channels showing '1' at that time point
%}


cfg.style             = 'straight';
cfg.comment           = 'no';
cfg.gridscale         = 512;
cfg.marker            = 'off'; % show the location of all channels?

cfg.parameter = 'stat'; % what field (in the stats output) to plot, e.g. selecting 'stat' will plot t-values

ft_topoplotER(cfg, stat);
set(gca,'fontsize',22); % set colorbar fontsize
%title(gca, 'T-values');


%% transform into planar gradient (based on actual erf amplitude) 

% first, define the 2 conds to be compared:
% here we look at main effect of lang in target window
GA_en = GA_erf.targetenstay;
GA_en.avg = (GA_erf.targetenstay.avg + GA_erf.targetenswitch.avg) / 2;
GA_ch = GA_erf.targetchstay;
GA_ch.avg = (GA_erf.targetchstay.avg + GA_erf.targetchswitch.avg) / 2;

% then, calc the diff btwn the 2 conds
cfg  = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
GA_diff = ft_math(cfg, GA_en, GA_ch);

% fill in the "grad" field so that ft_megplanar will work
GA_diff.grad = stat.grad;


% Compute the planar gradient at each sensor location
% in both the horizontal and the vertical direction 
cfg                 = [];
cfg.feedback        = 'yes';
cfg.method          = 'template';
cfg.neighbours      = neighbours;

cfg.planarmethod    = 'sincos';
planar              = ft_megplanar(cfg, GA_diff); % plug in any timelock struct

% Combine the horizontal and vertical components 
% of the planar gradient using Pythagoras rule
planarComb = ft_combineplanar([], planar);

% Plot the planar gradient
cfg.zlim = 'maxabs';
cfg.xlim = [start_time end_time]; % duration of the effect (as reported by ft_clusterplot)
                          % topography will be averaged over this interval
cfg.layout = lay;
cfg.colormap = colours;
cfg.colorbar = 'yes';

ft_topoplotER(cfg, planarComb);
set(gca,'fontsize',22); % set colorbar fontsize
