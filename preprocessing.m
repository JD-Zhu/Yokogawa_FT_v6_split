%{
Performs preprocessing for Yokogawa MEG data:
- filtering
- trigger-based trial definition (i.e. epoching)
%}

function [all_blocks, trialinfo_b] = preprocessing(SubjectFolder)
    
    % find all .con files
    files = dir([SubjectFolder '*.con']);
    % skip "_test.con" file
    for i = 1:length(files)
        if (strcmp(files(i).name(end-8:end), '_test.con') == 1)
            files(i) = [];
        end
    end

    % each cycle processes one '.con' file.
    % we do it in multiple steps (multiple calls to ft_prepreocessing)
    % to ensure things happen in the desired order
    for i = 1:length(files)
        rawfile = [SubjectFolder files(i).name];

        % ft_definetrial: defines the segments of data that will be read in by FT_PREPROCESSING
        cfg                      = [];
        cfg.headerfile           = rawfile;
        cfg.datafile             = rawfile;
        cfg.trialdef.triallength = Inf;
        cfg.trialdef.ntrials     = 1; % read in all data as a single segment, coz filtering should be done on continuous data
        cfg = ft_definetrial(cfg);

        % ft_preprocessing: reads in MEG data
        cfg.continuous = 'yes';
        cfg.bpfilter   = 'yes';
        cfg.bpfreq     = [0.5 40]; % bandpass filter
        alldata = ft_preprocessing(cfg);

        % Select channels 1-160 (i.e. MEG data)
        cfg         = [];
        cfg.channel = alldata.label(1:160);
        alldata = ft_selectdata(cfg, alldata);

        % deal with 50Hz line noise (necessary even after bandpass filter, coz the 50Hz noise is huge)
        cfg          = [];
        cfg.bsfilter = 'yes';
        cfg.bsfreq   = [49.5 50.5];
        alldata = ft_preprocessing(cfg, alldata);

        
        % browse raw data to remove noisy segments (e.g. clenched jaw)
        % http://www.fieldtriptoolbox.org/walkthrough/#visual-data-inspection
        %
        % Note: doing 'partial' rejection here because there is only 1 trial,
        % but this gives weird partial trials later on, and the final avg timecourse
        % (ERFs) are also partial - which is bad. 
        % Alternative is to do this step after epoching, however in that case
        % the data we inspect will be after filtering, so it's harder to
        % identify the artefacts (coz muscle artefacts are usually high freq)
        %
        % Maybe best to use ICA later instead of doing manual selection here
        % see remove_artefact_ICA.m (called by Yokogawa_FT_v6.m)
        %{
        cfg = [];
        cfg.viewmode = 'vertical';
        cfg.continous = 'yes';
        cfg.blocksize = 10; % display 10-sec segments
        cfg = ft_databrowser(cfg, alldata);
        cfg.artfctdef.reject = 'partial'; % this rejects complete trials, use 'partial' if you want to do partial artifact rejection
        alldata = ft_rejectartifact(cfg, alldata);
        %}
        
        
        %     Define trials using custom trialfun
        cfg                   = [];
        cfg.dataset           = rawfile;
        cfg.continuous        = 'yes';
        cfg.trialfun          = 'trig_fun_160_basic_v2';
        cfg.trialdef.prestim  = 1;         % pre-stimulus interval
        cfg.trialdef.poststim = 1;        % post-stimulus interval
        trialinfo_b(i) = ft_definetrial(cfg);
        alldata = ft_redefinetrial(trialinfo_b(i), alldata);

        cfg         = [];
        cfg.demean  = 'yes'; % subtracts the mean of the time window from all samples (i.e. centres the waveform on 0)
        cfg.detrend = 'yes'; % removes low-frequency drift
        block(i) = ft_preprocessing(cfg, alldata);


        % Create layout file for later (plotting)
        %{
        cfg      = [];
        cfg.grad = alldata.grad; % struct containing gradiometer definition
        lay      = ft_prepare_layout(cfg, alldata); % creates a 2-D layout of the channel locations
        %save([ResultsFolder 'lay.mat'], 'lay');
        %}
        
        % Prepare neighbours & save for use in stat analysis
        %{
        cfg_neighb        = [];
        cfg_neighb.method = 'triangulation'; % 'distance' may have some issues
                   % 'triangulation' seems to give each channel more neighbours 
                   % (about twice as many as 'distance' method gives)
        neighbours        = ft_prepare_neighbours(cfg_neighb, all_blocks);
        save([ResultsFolder 'neighbours.mat'], 'neighbours');
        %}

    end

    % combine data from all the blocks into one dataset
    all_blocks = block(1);
    for i = 2:length(block) % each cycle appends one block
        all_blocks = ft_appenddata([], all_blocks, block(i));
    end
    
end