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

% select which effect to plot
stat = cue_ttype;


%% plot topography

% load nice colourmap
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colours = colormap(flipud(brewermap(64, 'RdBu')));

% call topoplot
cfg                   = [];
cfg.layout            = lay;
cfg.colormap          = colours;
cfg.colorbar          = 'yes'; % shows the scales
%cfg.colorbar         = 'EastOutside';
cfg.zlim              = [-4 4];%'maxabs'; % set the scaling

cfg.xlim = [0.425 0.550]; % duration of the effect (as reported by ft_clusterplot)
                          % topography will be averaged over this interval

cfg.highlight         = 'on'; % highlight significant channels
cfg.highlightsymbol = '.';
cfg.highlightsize = 26;
cfg.highlightcolor = 'w';
miny = 98; % index of the time point you want (we selected 485ms - middle of the effect interval)
cfg.highlightchannel  = stat.label(ismember(stat.negclusterslabelmat(:,miny),1)); % find all the channels showing '1' at that time point

cfg.style             = 'straight';
cfg.comment           = 'no';
cfg.gridscale         = 512;
cfg.marker            = 'off'; % show the location of all channels?

cfg.parameter = 'stat'; % what field (in the stats output) to plot, e.g. selecting 'stat' will plot t-values

figure; ft_topoplotER(cfg, stat);
set(gca,'fontsize',22); % set colorbar fontsize
%title(gca, 'T-values');
