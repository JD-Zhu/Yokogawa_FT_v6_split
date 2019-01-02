%{
Version log:

v3: latest stable working version with all the features
v4: most EVAL statements are replaced by using dynamic field names rather
than dynamic variable names (script should run faster)
v5: added baseline correction after calc'ing erf -> later removed, in order to preserve cov matrix;
    added error-trial rejection, based on "errorsheet" from beh data checking;
    added downsampling (for saving);
    added computation of cov matrix in timelockanalysis (to enable creation of spatial filters in sourceanalysis)
v6: split whole script into 3 main stages of processing & each stage now loops thru all subjects;
    abstracted certain sections into its own function;
    created "common.m" - equivalent to a .h file which can be included/executed at the top of each script
%}

%%
close all
clear all % disable this line if u want breakpoints to work


% = Settings =
% Please adjust as required:

CHANNEL_REPAIR = false; % repair bad/rejected channels?
CALC_UNCLEANED_ERF = false; % calculate uncleaned erf? (for quality check of response-component rejection)

REMOVE_TRIGGER_ARTEFACT_ON_INDI_EPOCHS = false; % remove trigger-leak artefact? (spike around 55ms before cue onset & target onset)
REMOVE_TRIGGER_ARTEFACT_ON_AVG_ERF = false;
REMOVE_TRIGGER_LEAK_CHANNELS = false; % remove all channels affected by trigger leak (81-89)?


%%
% run the #define section
global DataFolder; global ResultsFolder; global filename_suffix; 
global eventnames;
common();

% check the Settings, and modify stuff accordingly
%if (channelrepair == true)
%    ResultsFolder = [ResultsFolder(1:end-2) '_channelrepair\\']; % modify ResultsFolder location
%end    

% set filenames for saving the output from each stage (so that we don't have to rerun the whole thing from beginning every time)
S1_output_filename = 'S1_preprocessed_data.mat'; % Stage 1 output (stored inside each Subject folder)
S2_output_filename = ['S2_after_visual_rejection' filename_suffix '.mat']; % Stage 2 output (stored inside each Subject folder)
S3_output_filename = ['_erf' filename_suffix '.mat']; % ERF output (stored in ResultsFolder for all subjects)

% enable access to 'SubjectID' from inside "trig_fun_160_...", so that 
% correct code_delay value can be set for each subject (30ms for first 5 subjects, 100ms for the others)
global SubjectID; 

% find all subject folders containing raw MEG recording
SubjectIDs = dir([DataFolder 'M*']);
SubjectIDs = {SubjectIDs.name}; % extract the names into a cell array
%SubjectIDs = {'M03-AG-2784'}; % or manually select which subjects to process


%% Stage 1: preprocessing & downsampling

for i = 1:length(SubjectIDs)
    
    SubjectID = cell2mat(SubjectIDs(i));
    SubjectFolder = [DataFolder SubjectID '\\'];
    S1_output_file = [SubjectFolder S1_output_filename];

    
    % if haven't already processed this stage before, do it now & save a copy
    if (exist(S1_output_file, 'file') ~= 2)    
        [all_blocks, trialinfo_b] = preprocessing(SubjectFolder);
        %TEMP: reject trigger artefact on raw continuous data
        %[all_blocks, trialinfo_b] = preprocessing2(SubjectFolder);

        % downsample the data for saving
        %all_blocks.time(1:end) = all_blocks.time(1); % this avoids numeric round off issues in the time axes upon resampling
        cfg            = [];
        cfg.resamplefs = 200; % sampling freq was 1000Hz, best to use a divisor of it (200Hz is commonly used)
        cfg.detrend    = 'no';
        all_blocks     = ft_resampledata(cfg, all_blocks);

        % SAVE preprocessed data - takes a while!!
        save(S1_output_file, 'all_blocks', 'trialinfo_b');
    end
end


%%  Stage 2: trial exclusions

for k = 1:length(SubjectIDs)
    
    SubjectID = cell2mat(SubjectIDs(k));
    SubjectFolder = [DataFolder SubjectID '\\'];
    S2_output_file = [SubjectFolder S2_output_filename];

    % if haven't already processed this stage before, do it now & save a copy
    if (exist(S2_output_file, 'file') ~= 2)    

        load([SubjectFolder S1_output_filename]);

        events_allBlocks = identify_event_types(SubjectID, trialinfo_b);
        
        events_allBlocks = exclude_beh_errors(SubjectID, events_allBlocks);

        
        % === ICA artefact removal ===

        % remove mouth-movement artefact by extracting main components from the "response" epochs
        % and projecting these out of all trials
        load([ResultsFolder 'lay.mat']);
        [all_blocks_clean, response_comp] = remove_artefact_ICA(all_blocks, events_allBlocks, lay, 'response');
        
        % remove trigger-leak artefact (if needed)
        % Note: if we put this here, then have to redo visual rejection
        %       alternatively, we put this in Stage 3 (results looked similar)
        %if (REMOVE_TRIGGER_ARTEFACT_ON_INDI_EPOCHS)
        %    load([ResultsFolder 'lay.mat']);
        %    [all_blocks_clean, trigger_comp] = remove_artefact_ICA(all_blocks_clean, events_allBlocks, lay, 'trigger');
        %end

        
        % === Reject Outlier Trials ===

        % Print out SubjectID so we know which subject we are working on
        fprintf(['\nCURRENT SUBJECT: ' SubjectID '\n\n']); 

        % Display visual trial summary to reject deviant trials.
        cfg              = [];
        cfg.feedback     = 'no'; % suppress console output (so that it's easy to find the SubjectID we printed out above)
        cfg.method       = 'summary';
        cfg.keepchannel  = 'no';
        cfg.keeptrial    = 'nan'; % we keep the rejected trials as 'NaN' here,
            % because if we remove them, that will change the indices of all subsequent trials,
            % which will no longer match the indices we are using in events_allBlocks
        all_blocks_clean = ft_rejectvisual(cfg, all_blocks_clean);
        
    
        save([SubjectFolder S2_output_filename], 'all_blocks_clean', 'events_allBlocks', 'response_comp'); 
        % 'all_blocks' was not changed in Stage 2, so don't need to save again
    end
end


%% Stage 3: time-domain analysis (i.e. compute erf)

for i = 1:length(SubjectIDs)
    
    SubjectID = cell2mat(SubjectIDs(i));
    SubjectFolder = [DataFolder SubjectID '\\'];
    S3_output_file = [ResultsFolder SubjectID S3_output_filename];

    % if haven't already processed this stage before, do it now & save a copy
    if (exist(S3_output_file, 'file') ~= 2)    
        
        % make sure we have a clean start (i.e. no leftover var contents from last subject)
        clear erf; clear erf_clean;
        clear trials; clear trials_clean;

        load([SubjectFolder S2_output_filename]);
        
        % remove trigger-leak artefact (if needed)
        if (REMOVE_TRIGGER_ARTEFACT_ON_INDI_EPOCHS)
            load([ResultsFolder 'lay.mat']);
            [all_blocks_clean, trigger_comp] = remove_artefact_ICA(all_blocks_clean, events_allBlocks, lay, 'trigger');
        end

        if (REMOVE_TRIGGER_LEAK_CHANNELS)
            cfg         = [];
            cfg.channel = {'all', '-AG083', '-AG087', '-AG088', '-AG082', '-AG084', '-AG086', '-AG081', '-AG085', '-AG089'}; % {'MEG'};
            all_blocks_clean = ft_selectdata(cfg, all_blocks_clean);
        end

        % perform channel repair if needed
        if (CHANNEL_REPAIR)
            load([ResultsFolder 'neighbours.mat']);
            all_labels = all_blocks_clean.cfg.channel; % full list of 160 labels
            all_blocks_clean = repair_bad_channels(all_blocks_clean, neighbours, all_labels);
        end    


        % === ft_redefine all event types (i.e. 8 real conditions + 'response' event) ===

        % in uncleaned data
        if (CALC_UNCLEANED_ERF)
            % load 'all_blocks'
            load([SubjectFolder S1_output_filename]);

            for j = 1:length(eventnames)
                cfg = [];
                cfg.trials = events_allBlocks.(eventnames{j});
                trials.(eventnames{j}) = ft_redefinetrial(cfg, all_blocks);
            end
        end
        
        % in cleaned data
        for j = 1:length(eventnames)
            cfg = [];
            cfg.trials = events_allBlocks.(eventnames{j});
            trials_clean.(eventnames{j}) = ft_redefinetrial(cfg, all_blocks_clean);
        end

        % delete vars that are no longer needed (to clear up some memory)
        %clear all_blocks;
        %clear all_blocks_clean;
        %clear trialinfo_b;

        %save([SubjectFolder 'before_computing_erf.mat'], 'trials', 'trials_clean', 'response_comp');


        % === compute ERFs ===

        % in uncleaned data (just for quality check of PCA component rejection)
        if (CALC_UNCLEANED_ERF)
            for j = 1:length(eventnames)
                cfg         = [];
                %cfg.nanmean = 'yes';
                %trials.(eventnames{j}) = ft_selectdata(cfg, trials.(eventnames{j})); % Do this because we kept bad trials as NaN
                erf.(eventnames{j}) = ft_timelockanalysis(cfg, trials.(eventnames{j})); % Do this to create average field and timelock struct
            end
        end
        
        % in cleaned data (compute erfs & cov matrices)
        [erf_clean, erf_cue_combined, erf_target_combined] = compute_ERF(trials_clean);

        % Do not perform baseline correction here, because it will remove the cov matrix from the timelock output:
        % "baseline correction invalidates previous covariance estimate, removing cov"
        %{
        % baseline correction
        for j = 1:length(eventnames)
            cfg = [];
            cfg.baseline = [-0.2 0];
            erf.(eventnames{j}) = ft_timelockbaseline(cfg, erf.(eventnames{j})); % uncleaned data
            erf_clean.(eventnames{j}) = ft_timelockbaseline(cfg, erf_clean.(eventnames{j})); % cleaned data
        end
        %}
        
        % remove trigger artefact on averaged ERF
        if REMOVE_TRIGGER_ARTEFACT_ON_AVG_ERF
            % here we use remove_artefact_ICA just to compute the trigger_comp
            % (we ignore the cleaned data it returns),
            % then manually call ft_rejectcomponent on the avg ERF for each cond
            load([ResultsFolder 'lay.mat']);
            [~, trigger_comp] = remove_artefact_ICA(all_blocks_clean, events_allBlocks, lay, 'trigger');
            
            % ft_rejectcomponent automatically converts timelock data to raw (and then
            % converts back) if you pass in one condition at a time.
            % the "cov" and "var" fields get removed, but we only need "cov" from 
            % erf_cue_combined & erf_target_combined anyway, not from the indi conditions.

            for j = 1:length(eventnames)
                cfg              = [];
                cfg.demean       = 'no';
                cfg.component    = 1:1; % project out the 1st principal component
                erf_clean.(eventnames{j}) = ft_rejectcomponent(cfg, trigger_comp, erf_clean.(eventnames{j}));
            end

            % Note: the covmatrix computed above (in compute_ERF.m) 
            % should still be valid, because it is not based on any data 
            % containing the trigger artefact (we modified cfg.covariancewindow
            % to 0~650ms following cue onset)
        end
        
        %### TODO: insert TF analysis here ###
        
        
        % SAVE all relevant variables from the workspace
        save(S3_output_file, 'SubjectFolder', ...
            'erf_clean', 'erf_cue_combined', 'erf_target_combined'); %'erf',       
    end

    
    % === Plot ERF & GFP (can use this to regen all plots from saved erf results) ===

    %load([ResultsFolder SubjectID S3_output_filename]);
    %load([ResultsFolder 'lay.mat']); 
    
    if (CALC_UNCLEANED_ERF)
        %plot_ERF(erf, erf_clean, lay, true, true);
    else % clean erf only
        %plot_ERF([], erf_clean, lay, false, false);
    end
    
end
