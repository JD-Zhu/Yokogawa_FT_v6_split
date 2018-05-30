%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% stats_ERF_TFCE.m
%
% statistical analysis on ERF data using the TFCE method:
% https://github.com/Mensen/ept_TFCE-matlab
%
%
% HOW TO REPORT THE CLUSTERS:
% "Analysis showed that Condition A had significantly higher amplitudes 
% compared to Condition B for 56 unique channels in the frontal region 
% for the time range from 280 - 320 ms (peak channel: Fz; T = 6.024, p = 0.002)."
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stats_ERF_TFCE()

    % run the #define section
    global conds_cue; global conds_target; global eventnames;
    global ResultsFolder; % all subjects' erf data are stored here
    global filename_suffix; % erf results file suffix
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
    files = dir([ResultsFolder '*_erf' filename_suffix '.mat']);

    % each cycle reads in one '.mat' file (ie. one subject's erf data)
    for i = 1:length(files)
        filename = [ResultsFolder files(i).name];
        temp = load(filename);
        erf_clean = temp.erf_clean;

        for j = 1:length(eventnames_8) % 4 conditions in cue & 4 conditions in target (total 8)
            allSubjects_erf.(eventnames_8{j}) = [allSubjects_erf.(eventnames_8{j}) erf_clean.(eventnames_8{j})];
        end
    end
    
    % if haven't calculated GA before, do that now & save
    if exist([ResultsFolder 'GA_erf_allConditions.mat'], 'file') ~= 2
        % CALCULATE the grand average (across all subjects) for each condition
        cfg = [];
        cfg.channel   = {'all'};
        cfg.latency   = 'all';
        cfg.parameter = 'avg';
        for j = 1:length(eventnames_8)
            GA_erf.(eventnames_8{j}) = ft_timelockgrandaverage(cfg, allSubjects_erf.(eventnames_8{j}){:});  
            % "{:}" means to use data from all elements of the variable
        end

        save([ResultsFolder 'GA_erf_allConditions.mat'], 'GA_erf');
    end

    
    % reformat data into eeglab format (which TFCE accepts):
    % each condition contains a 3d ("subject x channel x time") matrix
    for j = 1:length(eventnames_8) % loop thru each condition, convert it to eeglab format
        allSubjects_erf_eeglab.(eventnames_8{j}) = convert_FT_to_eeglab(allSubjects_erf.(eventnames_8{j})); 
    end
    
    
    %% select the time interval for analysis (e.g. 0~750ms),
    % by cropping this portion from the original timecourse (-1~1s).
    time_field = allSubjects_erf.cuechstay{1,1}.time; % get the original timecourse from any condition
    start_sample = find(time_field<0.0001, 1, 'last'); % find the index of the sample at 0ms (don't use the exact number '0', coz the actual time recorded is sth like -5.4845e-14 second
    end_sample = find(time_field>0.7499, 1, 'first'); % find the index of the sample at 750ms
    
    % crop the time_field
    time_field = time_field(start_sample:end_sample);
    
    % crop the data accordingly    
    data = allSubjects_erf_eeglab;
    for j = 1:length(eventnames_8)
        data.(eventnames{j}) = data.(eventnames{j})(:,:,start_sample:end_sample);
    end
    
    
    %% Statistical analysis using TFCE method

    fprintf('\n= STATS: Threshold-free cluster enhancement (TFCE method) =\n');

    % read the channel locations
    %addpath(genpath('C:\Users\43606024\Documents\MATLAB\eeglab14_1_1b\'));
    %chanlocs = readlocs([ResultsFolder 'chanlocs_XYZ.txt'], 'filetype','custom', 'format',{'X','Y','Z'});
    chanlocs = []; load([ResultsFolder 'chanlocs.mat']);
    

    % Run the statistical tests (interaction & 2 main effects)

    % Interaction (i.e. calc sw$ in each lang, then submit the 2 sw$ for comparison)
    fprintf('\nCUE window -> Testing lang x ttype interaction:\n');
    [timelock1, timelock2] = combine_conds_for_T_test('eeglab', 'interaction', data.cuechstay, data.cuechswitch, data.cueenstay, data.cueenswitch);
    [cue_interaction] = myWrapper_ept_TFCE(timelock1, timelock2, chanlocs);
    fprintf('\nTARGET window -> Testing lang x ttype interaction:\n');
    [timelock1, timelock2] = combine_conds_for_T_test('eeglab', 'interaction', data.targetchstay, data.targetchswitch, data.targetenstay, data.targetenswitch); %'2-1 vs 4-3');
    [target_interaction] = myWrapper_ept_TFCE(timelock1, timelock2, chanlocs);

    % Main effect of lang (collapse across stay-switch)
    fprintf('\nCUE window -> Main effect of lang:\n');
    [timelock1, timelock2] = combine_conds_for_T_test('eeglab', 'main_12vs34', data.cuechstay, data.cuechswitch, data.cueenstay, data.cueenswitch);
    [cue_lang] = myWrapper_ept_TFCE(timelock1, timelock2, chanlocs);
    fprintf('\nTARGET window -> Main effect of lang:\n');
    [timelock1, timelock2] = combine_conds_for_T_test('eeglab', 'main_12vs34', data.targetchstay, data.targetchswitch, data.targetenstay, data.targetenswitch); %'2-1 vs 4-3');
    [target_lang] = myWrapper_ept_TFCE(timelock1, timelock2, chanlocs);

    % Main effect of switch (collapse across langs)
    fprintf('\nCUE window -> Main effect of ttype:\n');
    [timelock1, timelock2] = combine_conds_for_T_test('eeglab', 'main_13vs24', data.cuechstay, data.cuechswitch, data.cueenstay, data.cueenswitch);
    [cue_ttype] = myWrapper_ept_TFCE(timelock1, timelock2, chanlocs);
    fprintf('\nTARGET window -> Main effect of ttype:\n');
    [timelock1, timelock2] = combine_conds_for_T_test('eeglab', 'main_13vs24', data.targetchstay, data.targetchswitch, data.targetenstay, data.targetenswitch); %'2-1 vs 4-3');
    [target_ttype] = myWrapper_ept_TFCE(timelock1, timelock2, chanlocs);

    save([ResultsFolder 'stats_TFCE.mat'], 'cue_interaction', 'cue_lang', 'cue_ttype', ...
                                           'target_interaction', 'target_lang', 'target_ttype', ...
                                           'time_field');
            
    % Alternatively, use the built-in ANOVA function (seems to be not as sensitive)
    %{
    % put data into 2x2 cell array for ANOVA
    cue = {data.cuechstay, data.cuechswitch; data.cueenstay, data.cueenswitch};
    ept_TFCE_ANOVA(cue, chanlocs);    
    target = {data.targetchstay, data.targetchswitch; data.targetenstay, data.targetenswitch};
    ept_TFCE_ANOVA(target, chanlocs);
    
    % show the results
    ept_ResultViewer(); % GUI for visualisation (you need to load the saved output file using the GUI)

    % to calc the cluster locations, we need to tweak the ANOVA output structure into one that contains a single effect
    factor = A; % select which main effect / interaction you want: A, B, AB
    Results_single.P_Values = Results.P_Values.(factor);
    Results_single.Obs = Results.Obs.(factor);
    Results_single.TFCE_Obs = Results.TFCE_Obs.(factor);
    Results_single.maxTFCE = Results.maxTFCE.(factor);
    [cluster_results] = ept_calculateClusters(Results_single, Info.Electrodes.ChannelNeighbours, 0.05);
    %}


    %% Find the effects & plot them
    
    % Can run this section alone to output a description & plot for each effect.
    % Requires the following files to be placed in the ResultsFolder:
    %   stats_TFCE.mat, GA_erf_allConditions.mat, channel_neighbours.mat, time_field.mat (if not contained in stats_TFCE.mat)

    % load the stats output
    stats = load([ResultsFolder 'stats_TFCE.mat']);
    % get time_field: used for reading out the actual timings
    if isfield(stats, 'time_field') % in new version, time_field is saved inside the stats output
        time_field = stats.time_field;
        stats = rmfield(stats, 'time_field'); % remove this field now
    else % in old version, read the time_field saved in results folder
        temp = load([ResultsFolder 'time_field.mat']); time_field = temp.time_field;
    end
    % load GA for plotting
    temp = load([ResultsFolder 'GA_erf_allConditions.mat']); GA = temp.GA_erf;
    % prepare for calling ept_calculateClusters
    channel_neighbours = []; load([ResultsFolder 'channel_neighbours.mat']);
    threshold = 0.05; % set your alpha here
 
    fprintf('\nThe following effects were detected:\n');

    % Automatically go thru all 6 stats output (cue/target lang/ttype/interxn),
    % check if any clusters were found for each effect
    stats_names = fieldnames(stats);
    for i = 1:length(stats_names) % each cycle handles one effect (e.g. cue_lang)
        stat_name = stats_names{i};
        fprintf('\n[%s]: ', stat_name);
        
        % produce info for each cluster (eg. which channels & time points formed each cluster)
        [cluster_results] = ept_calculateClusters(stats.(stat_name), channel_neighbours, threshold);
        
        % each cycle handles one cluster
        for index = 1:length(cluster_results)
            % get the channels & time points which formed this cluster
            [channels, samples] = find(cluster_results(index).cluster_locations); 
            channels = unique(channels); % remove the repeated entries
            samples = unique(samples);
            start_time = time_field(samples(1)); % read out the actual timing (in ms)
            end_time = time_field(samples(end)); %NOTE: we are assuming the effect is continuous here (which is prob true in most cases). But really should check this!! (which is why we output the samples / time points below)

            % output a description of this effect to console
            fprintf('\n  Cluster %d: channels%s at samples%s (%.f~%.f ms). Peak channel %d, peak sample %d, T = %f, p = %f\n', ...
                index, sprintf(' %d', channels), sprintf(' %d', samples), start_time*1000, end_time*1000, ...
                cluster_results(index).channel_peak, cluster_results(index).sample_peak, cluster_results(index).max_t_value, cluster_results(index).p_value_peak);

            % plot this cluster, overlaid onto GA
            figure('Name', sprintf('%s: cluster %d', stat_name, index)); hold on;
            
            cfg        = [];
            cfg.channel = channels;
            if strcmp(stat_name(1:3), 'cue') % this effect occurs in cue window
                %ft_singleplotER(cfg, GA.cuechstay, GA.cuechswitch, GA.cueenstay, GA.cueenswitch);
                %legend(eventnames_8(conds_cue));
                
                % combine conditions to show main effect
                stay = GA.cuechstay;
                sw = GA.cuechswitch;
                stay.avg = (GA.cuechstay.avg + GA.cueenstay.avg) / 2; % average across stay & switch
                sw.avg = (GA.cuechswitch.avg + GA.cueenswitch.avg) / 2; % average across stay & switch
                ft_singleplotER(cfg, stay, sw);
                legend({'stay', 'switch'});    
                xlim([-0.2 0.75]);
            elseif strcmp(stat_name(1:6), 'target') % this effect occurs in target window
                %ft_singleplotER(cfg, GA.targetchstay, GA.targetchswitch, GA.targetenstay, GA.targetenswitch);
                %legend(eventnames_8(conds_target));                
                
                % combine conditions to show main effect
                ch = GA.targetchstay;
                en = GA.targetenstay;
                ch.avg = (GA.targetchstay.avg + GA.targetchswitch.avg) / 2; % average across stay & switch
                en.avg = (GA.targetenstay.avg + GA.targetenswitch.avg) / 2; % average across stay & switch
                ft_singleplotER(cfg, ch, en);
                legend({'ch', 'en'});                
                xlim([-0.2 0.75]); 
            else % should never be here
                fprintf('Error: an effect is found, but its not in either cue nor target window.\n');
            end
            
            line([start_time start_time], ylim, 'Color','black'); % plot a vertical line at start_time
            line([end_time end_time], ylim, 'Color','black'); % plot a vertical line at end_time
            hold off;
        end

    end

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SUBFUNCTIONS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % wrapper function for calling ept_TFCE(), so that settings only need to be changed in one place
    function Results = myWrapper_ept_TFCE(data1, data2, chanlocs)
        
        % the input format to TFCE requires a "channel" dimension, so we check this
        % https://github.com/Mensen/ept_TFCE-matlab/issues/21
        if (size(data1, 2) == 1) % if there is only 1 channel (e.g. ROI result),
                                 % we fake a 2nd channel by making a copy of the 1st channel
            data1 = repmat(data1, [1,2,1]);
            data2 = repmat(data2, [1,2,1]);
            
            % we also need to treat this data as time-frequency data (so that TFCE won't look for a channel locations file)
            flag_ft = true;
            chanlocs = [];
            
        else % for normal (multi channel) data
            flag_ft = false;
        end       
        
        Results = ept_TFCE(data1, data2, ...
            chanlocs, ...
            'type', 'd', ...
            'flag_ft', flag_ft, ...
            'flag_tfce', true, ...
            'nPerm', 2000, ...
            'rSample', 200, ...
            'flag_save', false);
            %'saveName', [ResultsFolder 'TFCE_temp\\ept_' ROI_name '.mat']); % set a location to temporarily store the output. we don't need to save it, but if you don't set a location, it will litter arond your current directory
    end
end
