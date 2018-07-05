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
