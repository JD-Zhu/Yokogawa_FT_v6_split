%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% stats_ROI_TFCE.m
%
% statistical analysis on reconstructed ROI activities using the TFCE method:
% https://github.com/Mensen/ept_TFCE-matlab
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

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
end

% reformat again into eeglab format (which TFCE accepts):
% in each ROI|condition, have a .Summary field containing the "subject x channel x time" matrix
for k = 1:length(ROIs_label)
    ROI_name = ROIs_label{k};
    
    for j = 1:length(eventnames_8) % loop thru each condition, to create the 3d matrix for this cond
        data_for_this_cond = allSubjects_ROIs.(ROI_name).(eventnames_8{j});
        subj_chan_time = []; % initialise the 3d matrix for this condition
        
        for subject = 1:length(data_for_this_cond) % loop thru all subjects
            % because the format requires a "channel" dimension, fake that by making 2 copies of the only channel (so we now have 2 channels)
            chan_time = vertcat(data_for_this_cond{subject}.avg, data_for_this_cond{subject}.avg);
            % add this subject's "chan x time" matrix to the 3d matrix
            subj_chan_time = cat(3, subj_chan_time, chan_time); % concatenate along the 3rd dimension
        end
        % subj_chan_time is now the 3d matrix containing all subjects ("subject" being the 3rd dimension)
        % change the order of matrix dimensions to: subj x chan x time
        subj_chan_time = permute(subj_chan_time, [3 1 2]); 
        
        % store the 3d matrix in new variable (eeglab format), under the correct ROI & condition name
        allSubjects_ROIs_eeglab.(ROI_name).(eventnames_8{j}).Summary = subj_chan_time; 
        allSubjects_ROIs_eeglab.(ROI_name).(eventnames_8{j}).label = ROI_name;
        allSubjects_ROIs_eeglab.(ROI_name).(eventnames_8{j}).dimord = 'subj_chan_time';
    end
end


%% Statistical analysis using TFCE method

fprintf('\n= STATS: Threshold-free cluster enhancement (TFCE method) =\n');

% each cycle processes one ROI
for k = 1:length(ROIs_label)
    
    ROI_name = ROIs_label{k};
    fprintf(['\nROI: ' ROI_name '\n']);

    data = allSubjects_ROIs_eeglab.(ROI_name); % data for the current ROI
    
    % run TFCE
    Results = ept_TFCE(data.cuechstay.Summary, data.cuechswitch.Summary, ... %TODO: create real comparisons for main effects & interaction
        [], ...
        'type', 'd', ...
        'flag_ft', true, ...
        'flag_tfce', true, ... % set this to 'true' to use the TFCE method
        'nPerm', 1000, ...
        'rSample', 200, ...
        'saveName', [ResultsFolder_ROI 'TFCE_temp\\ept_' ROI_name '.mat']); % set a location to temporarily store the output. we don't need to save it, but if you don't set a location, it will litter arond your current directory

    %load([ResultsFolder_ROI 'TFCE_temp\\ept_' ROI_name '.mat']);
    effects = find(Results.P_Values < 0.05)
    if ~isempty(effects)
        % store into a list of effects found, and output at the end?
    end

end

%{
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
    %{
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
    %}
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
    cfg.numrandomization = 1000; % Rule of thumb: use 500, and double this number if it turns out 
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
fprintf('\nThe following effects were detected:\n');

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
            fprintf('%s has an effect in %s, between %.f~%.f ms.\n', ROI_name, stat_name, start_time*1000, end_time*1000); % convert units to ms

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

%}