% === clean the data (by running PCA and projecting the response component out of all trials) ===

function [all_blocks_clean, response_comp] = reject_response_component(all_blocks, events_allBlocks, lay)
    
    % ft_redefine the "response" trials & calc erf
    cfg        = [];
    cfg.trials = events_allBlocks.response; % list of "response" events
    response   = ft_redefinetrial(cfg, all_blocks);

    cfg          = [];
    response_erf = ft_timelockanalysis(cfg, response);

    %Run ICA on the "response" trials
    disp('About to run ICA using the SVD method')
    cfg           = [];
    cfg.method    = 'svd';
    response_comp = ft_componentanalysis(cfg, response_erf);
    %save([SubjectFolder 'response_comp.mat'],'response_comp','-v7.3')

    %Change the colourmap
    ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
    colormap(flipud(brewermap(64, 'RdBu'))) % change the colormap

    %Display the components identified by PCA - change layout as needed
    cfg          = [];
    cfg.viewmode = 'component';
    cfg.layout   = lay;
    %ft_databrowser(cfg, response_comp)

    % project certain components in the response erf out of all trials
    cfg              = [];
    cfg.component    = 1:3; % select top 5 components in the response erf
    all_blocks_clean = ft_rejectcomponent(cfg, response_comp, all_blocks); % reject these comps from all trials
    %all_blocks_clean = all_blocks; % noPCA version: cleaned data is same as uncleaned data
    
end