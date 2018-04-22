% === Identify the trial numbers belonging to each condition (ie. event type), as defined in the eventcodes{} array ===

% Chn stay: cue - channel 17, target - channel 18
% Chn sw:	cue - channel 19, target - channel 20
% Eng stay: cue - channel 21, target - channel 22
% Eng sw:	cue - channel 23, target - channel 24

function events_allBlocks = identify_event_types(SubjectID, trialinfo_b)
    
    % run the #define section
    global eventcodes; global eventnames;
    common();

    
    trialsgone = 0;

    % each cycle is one "block" (i.e. one '.con' file)
    for i = 1:length(trialinfo_b)
        for j = 1:length(eventnames)
            events.(eventnames{j}) = find(strcmp({trialinfo_b(i).event.value}, eventcodes{j,2})); 
            % 9 fields representing the 9 types of events
            % each field contains a list of all events belonging to this type (by matching event code)
            % NB: this is a tmp var, it gets overwritten in each cycle
        end

        if i == 1 % first block

            for j = 1:length(eventnames)
                % save the lists to a perm var, also transpose each list
                events_allBlocks.(eventnames{j}) = events.(eventnames{j})'; 
            end

        else % all other blocks
            trialsinblock = length(trialinfo_b(i-1).event); % how many "trials" (i.e. events) were identified in previous block
            trialsgone = trialsgone + trialsinblock; % add this number to the total number of "past" trials

            for j = 1:length(eventnames)
                events_allBlocks.(eventnames{j}) = [events_allBlocks.(eventnames{j}); events.(eventnames{j})' + trialsgone]; 
                % continue to append to the perm lists (stored in "events_allBlocks")                                                                                             
                % in the end, each perm list will contain all events of that type from all blocks
            end
        end
    end

    % special provision for SubjectID = M09-BL-2731
    % for whom I forgot to record MEG data at the very beginning of B1
    % (the first response may have been partially cut off, so just remove it)
    if strcmp(SubjectID, 'M09-BL-2731')
        events_allBlocks.response(1) = [];
    end

end