%
% preprocessing & time-domain analysis for vep data
%
%
%%
%close all
%clear all

% addpath to access custom functions in all subfolders
addpath(genpath(pwd));


% select SubjectFolder here
global SubjectID;
SubjectID = 'vep_QZ'; % select SubjectFolder here
SubjectFolder = ['..\\..\\vep\\' SubjectID '\\']; % relative path
ResultsFolder = '..\\..\\vep_results_ERF\\'; % no need to change this, all subjects' erf data will be stored here

eventcodes = {{'chnfemale'},{'17'};{'chnmale'},{'21'};{'engfemale'},{'25'};{'engmale'},{'29'}}; %...          % face cues
              %{'chnfemaleITI'},{'18'};{'chnmaleITI'},{'22'};{'engfemaleITI'},{'26'};{'engmaleITI'},{'30'}}; % fixation cross (ITI) after the face cue

eventnames = eventcodes(:,1); % extract a list of all event names
eventnames = [eventnames{:}]; % convert into strings


%% preprocessing

% find all .con files
files = dir([SubjectFolder '*.con']);
% skip "_test.con" file
for i = 1:length(files)
    if (strcmp(files(i).name(end-8:end), '_test.con') == 1)
        files(i) = [];
    end
end

% each cycle processes one '.con' file
for i = 1:length(files)
    rawfile = [SubjectFolder files(i).name];
    
    % ft_definetrial: defines the segments of data that will be read in by FT_PREPROCESSING
    cfg                      = [];
    cfg.headerfile           = rawfile;
    cfg.datafile             = rawfile;
    cfg.trialdef.triallength = Inf;
    cfg.trialdef.ntrials     = 1; % read in all data as a single segment
    
    cfg = ft_definetrial(cfg);
    
    % ft_preprocessing: reads in MEG data
    cfg.continuous = 'yes';
    cfg.bpfilter   = 'yes';
    cfg.bpfreq     = [0.5 30]; % bandpass filter
    
    alldata = ft_preprocessing(cfg);
    
    % Select channels 1-160 (i.e. MEG data)
    cfg         = [];
    cfg.channel = alldata.label(1:160);
    
    alldata = ft_selectdata(cfg, alldata);
    
    % deal with 50Hz line noise
    cfg          = [];
    cfg.bsfilter = 'yes';
    cfg.bsfreq   = [49.5 50.5];
    
    alldata = ft_preprocessing(cfg, alldata);
    
    % SAVE preprocessed data
    %save([SubjectFolder 'preprocessed_data_B' num2str(i) '.mat'], 'alldata');
    
    %     Create layout file for later + save
    cfg      = [];
    cfg.grad = alldata.grad; % struct containing gradiometer definition
    lay      = ft_prepare_layout(cfg, alldata); % creates a 2-D layout of the channel locations
    %save([SubjectFolder 'lay'], 'lay');
    
    %     Define trials using custom trialfun
    cfg                   = [];
    cfg.dataset           = rawfile;
    cfg.continuous        = 'yes';
    cfg.trialfun          = 'trig_fun_160_basic_v2';
    cfg.trialdef.prestim  = 0.5;      % pre-stimulus interval
    cfg.trialdef.poststim = 0.5;        % post-stimulus interval
    
    %trialinfo             = ft_definetrial(cfg);
    %eval(['trialinfo_b',num2str(i),' = trialinfo']);
    trialinfo_b(i) = ft_definetrial(cfg);
    alldata = ft_redefinetrial(trialinfo_b(i), alldata);
    
    cfg         = [];
    cfg.demean  = 'yes';
    cfg.detrend = 'yes';
    
    eval(['block_',num2str(i),' = ft_preprocessing(cfg, alldata);']);
    %block(i) = ft_preprocessing(cfg, alldata);

end
% preprocessing complete!

% combine data from all blocks into one dataset
%for i = 1:length(files) - 1
%    blocks{i} = ['block_',num2str(i),','];
%end
%blocks{i+1} = ['block_',num2str(i+1)];
%ie. blocks = {'block_1', 'block_2', etc};
blocks = {'block_1'}; % fix this up - there is only 1 block

cfg = [];
eval(['all_blocks = ft_appenddata(cfg,',cell2mat(blocks),');']);
%ie. all_blocks = ft_appenddata(cfg, block_1, block_2);

% now "all_blocks" contains all the data from the whole exp

%{
% prepare neighbours & save for use in stat analysis
cfg_neighb        = [];
cfg_neighb.method = 'distance'; % or 'triangulation'        
neighbours        = ft_prepare_neighbours(cfg_neighb, all_blocks);
save([ResultsFolder 'neighbours.mat'], 'neighbours');
%}

%%

% Chn stay: cue - channel 17, target - channel 18
% Chn sw:	cue - channel 19, target - channel 20
% Eng stay: cue - channel 21, target - channel 22
% Eng sw:	cue - channel 23, target - channel 24


% Create a list of trial numbers for each event type (as defined in the eventcodes{} array)

trialsgone = 0;

% each cycle is one "block" (i.e. one '.con' file)
for i = 1:length(trialinfo_b)
    for j = 1:length(eventcodes)
        %eval([cell2mat(eventcodes{j,1}),'= find(strcmp({trialinfo_b',num2str(i),'.event.value},','''',cell2mat(eventcodes{j,2}),'''','));']);
        events.(eventnames{j}) = find(strcmp({trialinfo_b(i).event.value}, eventcodes{j,2})); % 9 fields representing the 9 types of events
                                                                                            % each field contains a list of all events belonging to this type
                                                                                            % (by matching event code)
                                                                                            % NB. this is a tmp var, it gets overwritten in each cycle
    end
    if i == 1 % first block only
        for j = 1:length(eventcodes)
            %eval([cell2mat(eventcodes{j,1}),'_b=',cell2mat(eventcodes{j,1}),'''',';']);
            events_b.(eventnames{j}) = events.(eventnames{j})'; % save the lists to a perm var, also transpose each list
        end
    else % all other blocks
        %eval(['trialsinblock=length(trialinfo_b',num2str(i-1),'.event);']);
        trialsinblock = length(trialinfo_b(i-1).event); % how many "trials" (i.e. events) were identified in previous block
        trialsgone = trialsgone + trialsinblock; % add this number to the total number of "past" trials
        
        for j = 1:length(eventcodes)
            %eval([cell2mat(eventcodes{j,1}),'_b=[',cell2mat(eventcodes{j,1}),'_b;',cell2mat(eventcodes{j,1}),'''','+trialsgone];'])
            events_b.(eventnames{j}) = [events_b.(eventnames{j}); events.(eventnames{j})' + trialsgone]; % continue to append to the perm lists (stored in "events_b")
                                                                                                         % in the end, each perm list will contain all events of that type from all blocks
        end
    end
end

%{
% ft_redefine the "response" trials & calc erf
cfg        = [];
cfg.trials = events_b.response; % list of "response" events
response   = ft_redefinetrial(cfg, all_blocks);

cfg          = [];
response_erf = ft_timelockanalysis(cfg, response);

%Run ICA on the "response" trials
disp('About to run ICA using the SVD method')
cfg           = [];
cfg.method    = 'svd';
response_comp = ft_componentanalysis(cfg, response_erf);

%Change the colourmap
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64, 'RdBu'))) % change the colormap

%save([SubjectFolder 'response_comp.mat'],'response_comp','-v7.3')

%Display Components - change layout as needed
cfg          = [];
cfg.viewmode = 'component';
cfg.layout   = lay;
ft_databrowser(cfg, response_comp)

%**unmix here...is that necessary?***
cfg              = [];
cfg.component    = 1:5; % reject top 5 components
all_blocks_clean = ft_rejectcomponent(cfg, response_comp, all_blocks); % reject these comps from all trials
%}

% Reject Outlier Trials
% Display visual trial summary to reject deviant trials.
cfg              = [];
cfg.method       = 'summary';
cfg.keepchannel  = 'no';
cfg.keeptrial    = 'nan';
all_blocks = ft_rejectvisual(cfg, all_blocks);


% ft_redefine all other event types (to get ready for analysis)
for j = 1:length(eventcodes)
    cfg = [];
    %eval(['cfg.trials =', cell2mat(eventcodes{j,1}),'_b;']);
    %eval([cell2mat(eventcodes{j,1}),'= ft_redefinetrial(cfg,all_blocks);']);
    cfg.trials = events_b.(eventnames{j});
    trials.(eventnames{j}) = ft_redefinetrial(cfg, all_blocks);
end

% do the same for the cleaned data
for j = 1:length(eventcodes)
    cfg = [];
    %eval(['cfg.trials =', cell2mat(eventcodes{j,1}),'_b;']);
    %eval([cell2mat(eventcodes{j,1}),'_clean = ft_redefinetrial(cfg,all_blocks_clean);']);
    cfg.trials = events_b.(eventnames{j});
    trials_clean.(eventnames{j}) = ft_redefinetrial(cfg, all_blocks);
end

%% time-domain analysis (i.e. compute erf)
for j = 1:length(eventcodes)
    cfg         = [];
    cfg.nanmean = 'yes';
    %eval([cell2mat(eventcodes{j,1}),'_erf = ft_selectdata(cfg,',cell2mat(eventcodes{j,1}),');']);
    %eval([cell2mat(eventcodes{j,1}),'_erf = ft_timelockanalysis(cfg,',cell2mat(eventcodes{j,1}),'_erf);']);
    erf.(eventnames{j}) = ft_selectdata(cfg, trials.(eventnames{j})); % Do this because we kept bad trials as NaN
    erf.(eventnames{j}) = ft_timelockanalysis(cfg, erf.(eventnames{j})); % Do this to create average field and timelock struct
end
% the timelockanalysis output is not exactly in the format we want,
% therefore we compute it manually below. The call to ft_timelockanalysis
% above is just to create the output structure

% compute erf manually (by calc'ing averages) & fill in some missing struct?
for j = 1:length(eventcodes)
    %eval(['tmp=',cell2mat(eventcodes{j,1}),'.trial;']);
    tmp = trials.(eventnames{j}).trial;
    tmp_mean = nanmean(reshape([tmp{:}],size(tmp{1},1),size(tmp{1},2),size(tmp,2)),3); % calc averages
    %eval([cell2mat(eventcodes{j,1}),'_erf = ',cell2mat(eventcodes{j,1}),';']);
    %eval([cell2mat(eventcodes{j,1}),'_erf.avg = tmp_mean;']);
    %eval([cell2mat(eventcodes{j,1}),'_erf.time =',cell2mat(eventcodes{j,1}),'_erf.time{1};']);
    %eval([cell2mat(eventcodes{j,1}),'_erf = rmfield(',cell2mat(eventcodes{j,1}),'_erf,''trial''',');']);
    %eval([cell2mat(eventcodes{j,1}),'_erf.dimord = ''chan_time''',';']);
    erf.(eventnames{j}) = trials.(eventnames{j}); % copy the trial struct from "trials"
    erf.(eventnames{j}).avg = tmp_mean; % fill in the average we just calc'ed
    erf.(eventnames{j}).time = erf.(eventnames{j}).time{1};
    erf.(eventnames{j}) = rmfield(erf.(eventnames{j}), 'trial');
    erf.(eventnames{j}).dimord = 'chan_time';
end

% do the same for the cleaned data
for j = 1:length(eventcodes)
    %eval(['tmp=',cell2mat(eventcodes{j,1}),'_clean.trial;']);
    tmp = trials_clean.(eventnames{j}).trial;
    tmp_mean = nanmean(reshape([tmp{:}],size(tmp{1},1),size(tmp{1},2),size(tmp,2)),3);
    %eval([cell2mat(eventcodes{j,1}),'_erf_clean = ',cell2mat(eventcodes{j,1}),'_clean;']);
    %eval([cell2mat(eventcodes{j,1}),'_erf_clean.avg = tmp_mean;']);
    %eval([cell2mat(eventcodes{j,1}),'_erf_clean.time =',cell2mat(eventcodes{j,1}),'_erf_clean.time{1};']);
    %eval([cell2mat(eventcodes{j,1}),'_erf_clean = rmfield(',cell2mat(eventcodes{j,1}),'_erf_clean,''trial''',');']);
    %eval([cell2mat(eventcodes{j,1}),'_erf_clean.dimord = ''chan_time''',';']);
    erf_clean.(eventnames{j}) = trials_clean.(eventnames{j});
    erf_clean.(eventnames{j}).avg = tmp_mean;
    erf_clean.(eventnames{j}).time = erf_clean.(eventnames{j}).time{1};
    erf_clean.(eventnames{j}) = rmfield(erf_clean.(eventnames{j}), 'trial');
    erf_clean.(eventnames{j}).dimord = 'chan_time';
end
    
% baseline correction
cfg = [];
cfg.baseline = [-0.1 0];
erf.(eventnames{j}) = ft_timelockbaseline(cfg, erf.(eventnames{j})); % uncleaned data
erf_clean.(eventnames{j}) = ft_timelockbaseline(cfg, erf_clean.(eventnames{j})); % cleaned data

% SAVE all relevant variables in the workspace
%save([SubjectFolder 'data.mat'], 'eventcodes', 'lay', 'erf', 'erf_clean');
save([ResultsFolder SubjectID '_erf.mat'], 'SubjectFolder', 'eventcodes', 'eventnames', 'lay', 'erf', 'erf_clean');


%% FOR PLOTTING (can use this to regen all plots from saved variables)

% load workspace vars
%load([SubjectFolder 'data.mat']);
ResultsFolder = '..\\vep_results_all\\'; % should match the path selected at the top
%load([OutputFolder 'vep_Phillip' '_erf.mat']); % select which subject to load

% required configurations b4 any calls to ft_multiplotER()
cfg              = [];
cfg.showlabels   = 'yes';
cfg.fontsize     = 6;
cfg.layout       = lay;
cfg.baseline     = [-0.1 0]; % makes no diff if we've already done baseline correction earlier
cfg.baselinetype = 'absolute';

%{
% all pairwise comparisons btwn raw data (_erf) & after artefact removal (_erf_clean)
% just to see how good the artefact removal is (atm we reject components 1:5, we can adjust this for more/less removal)
for j = 1:length(eventcodes)
    figure;
    %eval(['ft_multiplotER(cfg,',cell2mat(eventcodes{j,1}),'_erf,',cell2mat(eventcodes{j,1}),'_erf_clean);']);
    ft_multiplotER(cfg, erf.(eventnames{j}), erf_clean.(eventnames{j}));
end
%}

% to compare the 4 conds
% (requires the list of 'cfg' assignments above, if running in console)
figure('Name','ft_multiplotER: erf_clean.chnfemale, erf_clean.chnmale, erf_clean.engfemale, erf_clean.engmale');
ft_multiplotER(cfg, erf_clean.chnfemale, erf_clean.chnmale, erf_clean.engfemale, erf_clean.engmale);
legend([eventcodes{1:4, 1}]);

%{
figure('Name','ft_multiplotER: erf_clean.chnfemaleITI, erf_clean.chnmaleITI, erf_clean.engfemaleITI, erf_clean.engmaleITI');
ft_multiplotER(cfg, erf_clean.chnfemaleITI, erf_clean.chnmaleITI, erf_clean.engfemaleITI, erf_clean.engmaleITI);
legend([eventcodes{5:8, 1}]);
%}

% Calc global averages across all sensors (GFP = global field potentials)
cfg        = [];
cfg.method = 'power';
for j = 1:length(eventcodes)
    %eval([cell2mat(eventcodes{j,1}),'_erf_clean_GFP = ft_globalmeanfield(cfg,',cell2mat(eventcodes{j,1}),'_erf_clean);']);
    erf_clean_GFP.(eventnames{j}) = ft_globalmeanfield(cfg, erf_clean.(eventnames{j}));
end

% plot global averages for cue-locked 
cfg = []; 
figure('Name','cue_GFP'); hold on
%for j = 1:length(eventcodes)
for j = 1:4
    %eval(['plot(',cell2mat(eventcodes{j,1}),'_erf_clean_GFP.time,',cell2mat(eventcodes{j,1}),'_erf_clean_GFP.avg)'])
    plot(erf_clean_GFP.(eventnames{j}).time, erf_clean_GFP.(eventnames{j}).avg);
end
legend([eventcodes{1:4, 1}])
%{
% plot global averages for ITI (fixation cross) 
figure('Name','fixation_GFP'); hold on
for j = 5:8
    %eval(['plot(',cell2mat(eventcodes{j,1}),'_erf_clean_GFP.time,',cell2mat(eventcodes{j,1}),'_erf_clean_GFP.avg)'])
    plot(erf_clean_GFP.(eventnames{j}).time, erf_clean_GFP.(eventnames{j}).avg);
end
legend([eventcodes{5:8, 1}])
%}














% Reject Trials
% Display visual trial summary to reject deviant trials.
% You need to load the mag + grad separately due to different scales
%
% cfg             = [];
% cfg.method      = 'summary';
% cfg.keepchannel = 'yes';
% data            = ft_rejectvisual(cfg, data);
%
%
% % Display Data
% Displaying the (raw) preprocessed MEG data
%
% cfg          = [];
% cfg.viewmode = 'vertical';
% ft_databrowser(cfg,data)
%
% Load the summary again so you can manually remove any deviant trials
% cfg             = [];
% cfg.method      = 'summary';
% cfg.keepchannel = 'yes';
% data            = ft_rejectvisual(cfg, data);
%
% data_clean_noICA = data
% save data_clean_noICA data_clean_noICA
% clear data_clean_noICA
% close all
%
% %% !!! ICA !!!
% % Downsample Data
%
% data_orig = data; %save the original CLEAN data for later use
% cfg = [];
% cfg.resamplefs = 150; %downsample frequency
% cfg.detrend = 'no';
% disp('Downsampling data');
% data = ft_resampledata(cfg, data_orig);
%
% % Run ICA
% disp('About to run ICA using the Runica method')
% cfg            = [];
% cfg.method     = 'fastica';
% comp           = ft_componentanalysis(cfg, data);
% save('comp.mat','comp','-v7.3')
%
% % Display Components - change layout as needed
% cfg = [];
% cfg.viewmode = 'component';
% ft_databrowser(cfg, comp)
%
% %% Remove components from original data
% %% Decompose the original data as it was prior to downsampling
% diary on;
% disp('Decomposing the original data as it was prior to downsampling...');
% cfg           = [];
% cfg.unmixing  = comp.unmixing;
% cfg.topolabel = comp.topolabel;
% comp_orig     = ft_componentanalysis(cfg, data_orig);
%
% %% The original data can now be reconstructed, excluding specified components
% % This asks the user to specify the components to be removed
% disp('Enter components in the form [1 2 3]')
% comp2remove = input('Which components would you like to remove?\n');
% cfg           = [];
% cfg.component = [comp2remove]; %these are the components to be removed
% data_clean    = ft_rejectcomponent(cfg, comp_orig,data_orig);
%
% %% Save the clean data
% disp('Saving data_clean...');
% save('data_clean','data_clean','-v7.3')
% diary off