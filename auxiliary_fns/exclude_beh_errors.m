% === exclude trials involving beh errors (based on errorsheet from manual checking) ===

function events_allBlocks = exclude_beh_errors(SubjectID, events_allBlocks)
    
    % get trial info & error info for all crit trials
    allCritTrials_table = read_errorsheet_v2(SubjectID);

    % separate the trials into 4 lists, 1 list per condition
    chstay_rows = strcmp(allCritTrials_table.lang,'Chn') & strcmp(allCritTrials_table.t_type,'Stay');
    chswitch_rows = strcmp(allCritTrials_table.lang,'Chn') & strcmp(allCritTrials_table.t_type,'Switch');
    enstay_rows = strcmp(allCritTrials_table.lang,'Eng') & strcmp(allCritTrials_table.t_type,'Stay');
    enswitch_rows = strcmp(allCritTrials_table.lang,'Eng') & strcmp(allCritTrials_table.t_type,'Switch');
    chstay = allCritTrials_table(chstay_rows,:);
    chswitch = allCritTrials_table(chswitch_rows,:);
    enstay = allCritTrials_table(enstay_rows,:);
    enswitch = allCritTrials_table(enswitch_rows,:);


    % first check to ensure #trials in each cond are the same 
    % btwn errorsheet & MEG data
    assert(height(chstay) == length(events_allBlocks.cuechstay));
    assert(height(chswitch) == length(events_allBlocks.cuechswitch));
    assert(height(enstay) == length(events_allBlocks.cueenstay));
    assert(height(enswitch) == length(events_allBlocks.cueenswitch));
    assert(height(chstay) == length(events_allBlocks.targetchstay));
    assert(height(chswitch) == length(events_allBlocks.targetchswitch));
    assert(height(enstay) == length(events_allBlocks.targetenstay));
    assert(height(enswitch) == length(events_allBlocks.targetenswitch));

    % sort the trial list for each condition (just in case they are not already sorted)
    chstay = sortrows(chstay, [1 2]); % sort the trials based on round_n then trial_n
    chswitch = sortrows(chswitch, [1 2]);
    enstay = sortrows(enstay, [1 2]);
    enswitch = sortrows(enswitch, [1 2]);

    % for each cond, find the indices of trials needing to be excluded,
    % then exclude from both the cue-locked and target-locked events
    remove = find(chstay.error); % any non-zero entry indicates exclusion
    events_allBlocks.cuechstay(remove) = [];
    events_allBlocks.targetchstay(remove) = [];
    remove = find(chswitch.error);
    events_allBlocks.cuechswitch(remove) = [];
    events_allBlocks.targetchswitch(remove) = [];
    remove = find(enstay.error);
    events_allBlocks.cueenstay(remove) = [];
    events_allBlocks.targetenstay(remove) = [];
    remove = find(enswitch.error);
    events_allBlocks.cueenswitch(remove) = [];
    events_allBlocks.targetenswitch(remove) = [];
    
end
