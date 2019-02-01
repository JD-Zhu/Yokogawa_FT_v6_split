% NOTE: This fn is still WIP! may or may not work as expected...
%
% ICA to remove eye blinks & other large muscle artefacts

function [data_clean] = remove_artefact_ICA(data, lay)

    disp('About to run ICA using the Runica method')
    cfg            = [];
    cfg.method     = 'fastica';
    comp           = ft_componentanalysis(cfg, data);
    %save('comp.mat','comp','-v7.3')

    %Change the colourmap
    ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
    colormap(flipud(brewermap(64, 'RdBu'))) % change the colormap

    %Display the components identified by ICA
    %
    cfg          = [];
    cfg.viewmode = 'component';
    cfg.layout   = lay;
    ft_databrowser(cfg, comp)
    drawnow; pause;
    %

    
    diary on;

    % Ask user to specify the components to be removed
    disp('Enter components in the form [1 2 3]')
    comps2remove = input('Which components would you like to remove?\n');

    % Remove the components
    cfg           = [];
    cfg.component = [comps2remove]; %these are the components to be removed
    data_clean    = ft_rejectcomponent(cfg, comp, data);

    % Save the clean data
    %disp('Saving data_clean...');
    %save('data_clean','data_clean','-v7.3')
    
    diary off
end