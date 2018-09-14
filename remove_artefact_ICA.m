% === artefact removal using PCA (projecting the artefact components out of all trials) ===
%
% @param type_of_artefact:
%       'response' = mouth movement artefact during speech
%       'trigger'  = trigger leak artefact (spike around 55ms prior to cue onset & target onset)
%
function [data_clean, artefact_comp] = remove_artefact_ICA(data, events_allBlocks, lay, type_of_artefact)
    % run the #define section
    global eventnames_8; global conds_cue; 
    common();

    
    clear artefact; % ensures no contamination
    
    % grab the relevant "trials" which characterise the artefact
    if strcmp(type_of_artefact, 'response')
        % take all the "response" epochs
        cfg        = [];
        cfg.trials = events_allBlocks.response; % list of "artefact" events
        %cfg.toilim = [-0.2 0.5]; % -200 to 500ms around speech onset
        artefact   = ft_redefinetrial(cfg, data);
        
    elseif strcmp(type_of_artefact, 'trigger')
        % take all the cue epochs
        % and extract the section containing the trigger artefact (-65 to -45ms)
        eventnames_cue = eventnames_8(conds_cue);
        %eventnames_cue = eventnames_8; % if you want to use both pre-cue & pre-target as definition of trigger artefact
        for j = 1:length(eventnames_cue)
            cfg = [];
            cfg.trials = events_allBlocks.(eventnames_cue{j});
            cfg.toilim = [-0.065 -0.040]; % -65 to -45ms
            trials = ft_redefinetrial(cfg, data);

            % append to existing trials
            if (j == 1)
                artefact = trials;
            else
                artefact = ft_appenddata([], artefact, trials);
            end
        end
    else
        error('Error in remove_artefact_ICA: type of artefact was not specified correctly.')
    end

    
    % We kept bad trials as NaN, need to exclude them before SVD
    good_trials_idx = find(~isnan(cell2mat(cellfun(@(isgood)isgood(1),artefact.trial,'uni',0)))); % Just need to evaluate the first element as all samples in bad trial are NaN
    cfg             = [];
    cfg.trials      = good_trials_idx;
    artefact        = ft_redefinetrial(cfg, artefact);
        
    % calc erf
    cfg          = [];
    artefact_erf = ft_timelockanalysis(cfg, artefact);

    %Run ICA on the "artefact" erf
    disp('About to run ICA using the SVD method')
    cfg           = [];
    cfg.method    = 'svd';
    artefact_comp = ft_componentanalysis(cfg, artefact_erf);
    %save([SubjectFolder 'artefact_comp.mat'],'artefact_comp','-v7.3')

    %Change the colourmap
    ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
    colormap(flipud(brewermap(64, 'RdBu'))) % change the colormap

    %Display the components identified by PCA - change layout as needed
    cfg          = [];
    cfg.viewmode = 'component';
    cfg.layout   = lay;
    %ft_databrowser(cfg, artefact_comp)
    %drawnow; pause;

    % project certain components in the artefact erf out of all trials
    if strcmp(type_of_artefact, 'response')
        cfg              = [];
        cfg.component    = 1:5; % select top 5 components in the artefact erf
        data_clean = ft_rejectcomponent(cfg, artefact_comp, data); % reject these comps from all trials
        %data_clean = data; % noPCA version: cleaned data is same as uncleaned data
    elseif strcmp(type_of_artefact, 'trigger')
        cfg              = [];
        cfg.component    = 1:1; % project out the 1st principal component
        data_clean = ft_rejectcomponent(cfg, artefact_comp, data); % reject these comps from all trials
    end    
end