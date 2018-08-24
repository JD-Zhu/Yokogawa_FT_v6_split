%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% stats_ROI_TFCE.m
%
% statistical analysis on reconstructed ROI activities using the TFCE method:
% https://github.com/Mensen/ept_TFCE-matlab
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stats_ROI_TFCE()
    
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
        ROI_activity = []; load(filename); % initialise the variable first to enable loading without error
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
    % each ROI|condition contains a 3d ("subject x channel x time") matrix
    for k = 1:length(ROIs_label)
        ROI_name = ROIs_label{k};

        for j = 1:length(eventnames_8) % loop thru each condition, convert it to eeglab format
            allSubjects_ROIs_eeglab.(ROI_name).(eventnames_8{j}) = convert_FT_to_eeglab(allSubjects_ROIs.(ROI_name).(eventnames_8{j})); 
            %allSubjects_ROIs_eeglab.(ROI_name).(eventnames_8{j}).label = ROI_name;
            %allSubjects_ROIs_eeglab.(ROI_name).(eventnames_8{j}).dimord = 'subj_chan_time';
        end
    end

    %% select the time interval for analysis (e.g. 0~750ms),
    % by cropping this portion from the original timecourse (-1~1s).
    time_field = allSubjects_ROIs.RIFG.cuechstay{1,1}.time; % get the original timecourse from any condition
    start_sample = find(time_field<0.0001, 1, 'last'); % find the index of the sample at 0ms (don't use the exact number '0', coz the actual time recorded is sth like -5.4845e-14 second
    end_sample = find(time_field>0.7499, 1, 'first'); % find the index of the sample at 750ms

    % crop the time_field
    time_field = time_field(start_sample:end_sample);
    % If you change the time interval above, should save a new copy:
    %save([ResultsFolder_ROI 'time_field.mat'], 'time_field');

    % crop the data accordingly 
    for k = 1:length(ROIs_label)
        ROI_name = ROIs_label{k};
        for j = 1:length(eventnames_8)
            allSubjects_ROIs_eeglab.(ROI_name).(eventnames{j}) = allSubjects_ROIs_eeglab.(ROI_name).(eventnames{j})(:,:,start_sample:end_sample);
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
                
        % Interaction (i.e. calc sw$ in each lang, then submit the 2 sw$ for comparison)
        fprintf('\nCUE window -> Testing lang x ttype interaction:\n');
        [timelock1, timelock2] = combine_conds_for_T_test('eeglab', 'interaction', data.cuechstay, data.cuechswitch, data.cueenstay, data.cueenswitch);
        [cue_interaction.(ROI_name)] = myWrapper_ept_TFCE(timelock1, timelock2);
        fprintf('\nTARGET window -> Testing lang x ttype interaction:\n');
        [timelock1, timelock2] = combine_conds_for_T_test('eeglab', 'interaction', data.targetchstay, data.targetchswitch, data.targetenstay, data.targetenswitch); %'2-1 vs 4-3');
        [target_interaction.(ROI_name)] = myWrapper_ept_TFCE(timelock1, timelock2);
    
        % Main effect of lang (collapse across stay-switch)
        fprintf('\nCUE window -> Main effect of lang:\n');
        [timelock1, timelock2] = combine_conds_for_T_test('eeglab', 'main_12vs34', data.cuechstay, data.cuechswitch, data.cueenstay, data.cueenswitch);
        [cue_lang.(ROI_name)] = myWrapper_ept_TFCE(timelock1, timelock2);
        fprintf('\nTARGET window -> Main effect of lang:\n');
        [timelock1, timelock2] = combine_conds_for_T_test('eeglab', 'main_12vs34', data.targetchstay, data.targetchswitch, data.targetenstay, data.targetenswitch); %'2-1 vs 4-3');
        [target_lang.(ROI_name)] = myWrapper_ept_TFCE(timelock1, timelock2);

        % Main effect of switch (collapse across langs)
        fprintf('\nCUE window -> Main effect of ttype:\n');
        [timelock1, timelock2] = combine_conds_for_T_test('eeglab', 'main_13vs24', data.cuechstay, data.cuechswitch, data.cueenstay, data.cueenswitch);
        [cue_ttype.(ROI_name)] = myWrapper_ept_TFCE(timelock1, timelock2);
        fprintf('\nTARGET window -> Main effect of ttype:\n');
        [timelock1, timelock2] = combine_conds_for_T_test('eeglab', 'main_13vs24', data.targetchstay, data.targetchswitch, data.targetenstay, data.targetenswitch); %'2-1 vs 4-3');
        [target_ttype.(ROI_name)] = myWrapper_ept_TFCE(timelock1, timelock2);

    end
    
    save([ResultsFolder_ROI 'stats_TFCE.mat'], 'cue_interaction', 'cue_lang', 'cue_ttype', 'target_interaction', 'target_lang', 'target_ttype');

    
    %% Find the effects & plot them

    % Automatically check all the stats output & read out the time interval
    % of each effect (from the stat.P_Values field)

    stats = load([ResultsFolder_ROI 'stats_TFCE.mat']);
    temp = load([ResultsFolder_ROI 'time_field.mat']); time_field = temp.time_field; % used for reading out the time interval of effect
    temp = load([ResultsFolder_ROI 'GA.mat']); GA = temp.GA; % if you don't have this file, run stats_ROI.m to obtain it
    fprintf('\nThe following effects were detected:\n');

    % loop thru all 6 stats output (cue/target lang/ttype/interxn) and loop thru all ROIs in each,
    % check if any p-values in the results are significant (these are the effects)
    stats_names = fieldnames(stats);
    for i = 1:length(stats_names) % each cycle handles one effect (e.g. cue_lang)
        stat_name = stats_names{i};
        ROIs_names = fieldnames(stats.(stat_name)); % get the list of ROI names
        for k = 1:length(ROIs_names) % each cycle handles one ROI
            ROI_name = ROIs_names{k};
            pvalues = stats.(stat_name).(ROI_name).P_Values; % get the p-values
            pvalues = pvalues(1,:); % look at the 1st channel only, since the 2nd channel is a "fake" one we created (identical to 1st channel)
            
            % if any p-value is sig, that's an effect
            effect = find(pvalues < 0.05);
            % if there is an effect, we print it out
            if ~isempty(effect) 
                time_points = sprintf(' %d', effect);
                start_time = time_field(effect(1));
                end_time = time_field(effect(end)); %NOTE: we are assuming the effect is continuous here (which is prob true in most cases). But really should check this!! (which is why we output the samples / time points below)
                fprintf('%s has an effect in %s, between %.f~%.f ms (significant at samples:%s). p = %f\n', ROI_name, stat_name, start_time*1000, end_time*1000, time_points, min(pvalues(effect))); % convert units to ms

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

end
