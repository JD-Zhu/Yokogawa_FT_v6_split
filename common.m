% the common #define section for all scripts

% Warning: all of the global vars below are CONSTANTS.
% UNDER NO CIRCUMSTANCES should their values be assigned/modified in other scripts.

function [] = common()

    % specify all paths as absolute paths, to avoid any issues when we 'cd' into diff folders    
    global DataFolder; global ResultsFolder; global ResultsFolder_ROI; global ResultsFolder_Source;
    DataFolder = [pwd '\\..\\..\\RAW_DATA\\'];
    ResultsFolder = [pwd '\\..\\..\\results_ERF\\']; % all subjects' erf data will be stored here
    ResultsFolder_ROI = [pwd '\\..\\..\\results_ROI\\']; % all subjects' ROI source-reconstruction results will be stored here
    ResultsFolder_Source = [pwd '\\..\\..\\results_SOURCE\\']; % all subjects' source localisation results
    
    global filename_suffix; % select which pre-processing option: noPCA, reject components 1:3, or normal (reject components 1:5)
    % also need to change the last few lines in reject_response_component.m & load correct result files into the ResultsFolder
    filename_suffix = ''; % '_noPCA'; %'_rejectTop3'; %'';
    
    % =================================================================
    
    % trigger events (DO NOT change the order of this list)
    global eventcodes; global eventnames; global eventnames_8;
    eventcodes = {{'cuechstay'},{'17'};{'cuechswitch'},{'19'};{'cueenstay'},{'21'};{'cueenswitch'},{'23'}; ...
        {'targetchstay'},{'18'};{'targetchswitch'},{'20'};{'targetenstay'},{'22'};{'targetenswitch'},{'24'};{'response'},{'30'}};
    eventnames = eventcodes(:,1); % extract a list of all event names
    eventnames = [eventnames{:}]; % convert into strings
    
    % for ease of reference to the conditions in cue window & target window
    global conds_cue; global conds_target;
    conds_cue = 1:4;
    conds_target = 5:8;
    eventnames_8 = eventnames([conds_cue conds_target]); % 8 actual event types
    
    % =================================================================
    
    % addpath to access custom functions in all subfolders
    addpath(genpath(pwd));
    
    % =================================================================
    
    % = Settings = %

    % Plot shaded patch around ROI time course? 
    % Options: 'no', 'SEM', 'STDEV', 'CI_95'
    % (note: SEM < 95% CI < STDEV)
    global PLOT_SHADE;
    PLOT_SHADE = 'SEM';

    % colours for time course plots (one colour for each condition):
    % need to specify manually because we plot each cond separately, and simply 
    % using default colourmap makes all lines the same colour when calling boundedline()
    global colours;
    colours = ['b', 'r', 'y', 'm', 'b', 'r', 'y', 'm'];

    % toolbox to plot shaded boundary around each timecourse (diff paths for diff computers)
    addpath(genpath('C:\Users\Judy\Documents\MATLAB\kakearney-boundedline-pkg-50f7e4b'));
    addpath(genpath('C:\Users\43606024\Documents\MATLAB\kakearney-boundedline-pkg-50f7e4b'));
    
    % toolbox to save figure exactly as it appears on screen
    addpath(genpath('C:\Users\Judy\Documents\MATLAB\altmany-export_fig-9676767'));
    addpath(genpath('C:\Users\43606024\Documents\MATLAB\altmany-export_fig-9676767'));

end