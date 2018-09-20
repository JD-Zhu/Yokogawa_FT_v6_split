 function [all_blocks, trialinfo_b] = preprocessing2(SubjectFolder)
    global ResultsFolder;
    common();
 
    % find all .con files
    files = dir([SubjectFolder '*.con']);
    % skip "_test.con" file
    for i = 1:length(files)
        if (strcmp(files(i).name(end-8:end), '_test.con') == 1)
            files(i) = [];
        end
    end
    
    for i = 1:length(files)
        rawfile = [SubjectFolder files(i).name];

        % ft_definetrial: defines the segments of data that will be read in by FT_PREPROCESSING
        cfg                      = [];
        cfg.headerfile           = rawfile;
        cfg.datafile             = rawfile;
        cfg.trialdef.triallength = Inf;
        cfg.trialdef.ntrials     = 1; % read in all data as a single segment, coz filtering should be done on continuous data
        cfg = ft_definetrial(cfg);
        
        cfg.continuous = 'yes';
        all_blocks = ft_preprocessing(cfg);
        
        % Select channels 1-160 (i.e. MEG data)
        cfg         = [];
        cfg.channel = all_blocks.label(1:160);
        all_blocks = ft_selectdata(cfg, all_blocks);
        
        block(i) = all_blocks;
    end
    
    % combine data from all the blocks into one dataset
    all_blocks = block(1);
    for i = 2:length(block) % each cycle appends one block
        cfg = [];
        all_blocks = ft_appenddata(cfg, all_blocks, block(i));
    end

    
    % Run PCA on the raw continuous data
    disp('About to run PCA using the SVD method')
    cfg           = [];
    cfg.method    = 'svd';
    comp = ft_componentanalysis(cfg, all_blocks);
    %save([SubjectFolder 'artefact_comp.mat'],'artefact_comp','-v7.3')

    %Change the colourmap
    ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
    colormap(flipud(brewermap(64, 'RdBu'))) % change the colormap

    %Display the components identified by PCA - change layout as needed
    cfg          = [];
    cfg.viewmode = 'component';
    %cfg.channels = {'AG083', 'AG087', 'AG088', 'AG082', 'AG084', 'AG086'};
    %cfg.plotlabels              = 'yes';
    load([ResultsFolder 'lay.mat']);
    cfg.layout   = lay;
    ft_databrowser(cfg, comp)
    drawnow; pause;
    
    
    % highlight channens on the component topo
    figure;
    %Change the colourmap
    ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
    colormap(flipud(brewermap(64, 'RdBu'))) % change the colormap


    cfg           = [];
    cfg.component = 1:15;       % specify the component(s) that should be plotted
    cfg.layout    = lay; % specify the layout file that should be used for plotting
    cfg.comment   = 'no';
    cfg.highlight = 'on';
    cfg.highlightchannel   = {'AG083', 'AG087', 'AG088', 'AG082', 'AG084', 'AG086', 'AG081', 'AG085','AG089'};
    ft_topoplotIC(cfg, comp)
    
    % TEMP: visually identify the components to reject, enter those into
    % the section below:
    
    % project the identified components out of raw continuous data
    cfg              = [];
    cfg.component    = 4:6;
    all_blocks = ft_rejectcomponent(cfg, comp, all_blocks);
    
    
    
    
    % ft_preprocessing: reads in MEG data
    cfg = [];
    cfg.bpfilter   = 'yes';
    cfg.bpfreq     = [0.5 40]; % bandpass filter
    all_blocks = ft_preprocessing(cfg);

    % deal with 50Hz line noise (necessary even after bandpass filter, coz the 50Hz noise is huge)
    cfg          = [];
    cfg.bsfilter = 'yes';
    cfg.bsfreq   = [49.5 50.5];
    all_blocks = ft_preprocessing(cfg, all_blocks);

    %     Define trials using custom trialfun
    cfg                   = [];
    cfg.dataset           = rawfile;
    cfg.continuous        = 'yes';
    cfg.trialfun          = 'trig_fun_160_basic_v2';
    cfg.trialdef.prestim  = 1;         % pre-stimulus interval
    cfg.trialdef.poststim = 1;        % post-stimulus interval
    trialinfo_b(i) = ft_definetrial(cfg);
    all_blocks = ft_redefinetrial(trialinfo_b(i), all_blocks);

    cfg         = [];
    cfg.demean  = 'yes'; % subtracts the mean of the time window from all samples (i.e. centres the waveform on 0)
    cfg.detrend = 'yes'; % removes low-frequency drift
    all_blocks = ft_preprocessing(cfg, all_blocks);
 end