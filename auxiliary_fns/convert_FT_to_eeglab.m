% Convert data from fieldtrip format to eeglab format, so that they can be
% accepted by various stats toolbox (e.g. TFCE)
%
%    'eeglab' format = 3d matrix of values (subj x chan x time)
% 'fieldtrip' format = cell array of subjects, each cell containing the timelock struct for one subject (.avg field contains chan x time matrix)
%
% @param data: can be ERF result or ROI result
%
function subj_chan_time = convert_FT_to_eeglab(data)
    subj_chan_time = []; % initialise the 3d matrix

    for subject = 1:length(data) % loop thru all subjects
        % add this subject's "chan x time" matrix to the 3d matrix
        chan_time = data{subject}.avg;
        subj_chan_time = cat(3, subj_chan_time, chan_time); % concatenate along the 3rd dimension
    end
    
    % subj_chan_time is now the 3d matrix containing all subjects ("subject" being the 3rd dimension)
    % change the order of matrix dimensions to: subj x chan x time
    subj_chan_time = permute(subj_chan_time, [3 1 2]); 
end