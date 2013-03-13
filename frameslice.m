function frameslice(dir_in, file_in, num_slice, dir_out)
% FRAMESLICE Slice KU radar data frames into smaller pieces.
%
%   FRAMESLICE(DIR_IN, FILE_IN, NUM_SLICE, DIR_OUT) slices the original KU
%   data frames that match FILE_IN, which can contain wildcards (*), and
%   stored in directory DIR_IN, into NUM_SLICE slices and saves the output
%   in DIR_OUT.
%
% Joe MacGregor (UTIG)
% Last updated: 02/01/13

% warning/error checks
if (nargin < 4)
    error('frameslice:inputs', 'Not enough input arguments.')
end
if ~ischar(dir_in)
    error('frameslice:dirinstr', 'Frame input directory (DIR_IN) is not a string.')
end
if (~isempty(dir_in) && ~exist(dir_in, 'dir'))
    error('frameslice:nodirin', 'Frame input directory (DIR_IN) does not exist.')
end
if ~ischar(file_in)
    error('frameslice:fileinstr', 'Frame filename (FILE_IN) is not a string.');
end
if (~isnumeric(num_slice) || (length(num_slice) > 1))
    error('frameslice:numslicenum', 'Number of slices (NUM_SLICE) is not a scalar.')
end
if mod(num_slice, 1)
    num_slice               = round(num_slice);
    warning('frameslice:roundnumslice', ['Number of slices (NUM_SLICE) rounded to ' num2str(num_slice) '.'])
end
if (num_slice <= 1)
    error('frameslice:numslicesize', 'Number of slices (NUM_SLICE) must be greater than one.')
end
if ~ischar(dir_out)
    error('frameslice:diroutstr', 'Output directory (DIR_OUT) is not a string.')
end
if (~isempty(dir_out) && ~exist(dir_out, 'dir'))
    error('frameslice:nodirout', 'Output directory (DIR_OUT) does not exist.')
end

file_in_all                 = dir([dir_in file_in '.mat']);
file_in_all                 = {file_in_all.name};
num_in                      = length(file_in_all); % number of files
if ~num_in
    error('frameslice:nofiles', ['No files matching ' file_in '.mat found in ' dir_in '.'])
end

disp(['Slicing ' num2str(num_in) ' files in ' dir_in ' into ' num2str(num_slice) ' slices each...'])

vars2slice                  = {'Bottom' 'Data' 'Elevation' 'GPS_time' 'Latitude' 'Longitude' 'Surface'}; % old format

for ii = 1:num_in
    
    tmp1                    = load([dir_in file_in_all{ii}]);
    
    % remove overlapping data
    if (ii > 1)
        if any(tmp1.GPS_time < tmp1_old.GPS_time(end));
            ind_tmp         = find(tmp1.GPS_time > tmp1_old.GPS_time(end)); %#ok<EFIND>
            if ~isempty(ind_tmp)
                for jj = 1:length(vars2slice)
                    eval(['tmp1.' vars2slice{jj} ' = tmp1.' vars2slice{jj} '(:, ind_tmp);'])
                end
            else
                disp(['Skipping ' file_in_all{ii} ' because it overlaps completely in time with the previous data file (' file_in_all{ii - 1} ').'])
                continue % no points in this data file that don't overlap with previous file (whoa)
            end
        end
    end
    
    % slicing prep
    num_trace               = length(tmp1.Latitude); % length of dataset
    ind_slice               = round(linspace(1, num_trace, (num_slice + 1)));
    ind_slice(end)          = num_trace + 1; % adjust the last index to simplify loop
    
    Time                    = tmp1.Time; % two-way traveltime doesn't need slicing
    Depth                   = tmp1.Depth; % same for depth
    
    % loop through slices
    for jj = 1:num_slice
        
        curr_ind            = ind_slice(jj):(ind_slice(jj + 1) - 1); % indices within original dataset to keep for this loop iteration
        
        for kk = 1:length(vars2slice) % go through each variable and get the values of the current indices
            eval([vars2slice{kk} ' = tmp1.' vars2slice{kk} '(:, curr_ind);'])
        end
        
        if isfield(tmp1, 'param_records') % new format
            
            if isfield(tmp1, 'array_param')
                array_param = tmp1.array_param; %#ok<*NASGU>
            else
                array_param = [];
            end
            if isfield(tmp1, 'param_combine_wf_chan')
                param_combine_wf_chan ...
                            = tmp1.param_combine_wf_chan;
            else
                param_combine_wf_chan ...
                            = [];
            end
            if isfield(tmp1, 'param_csarp')
                param_csarp = tmp1.param_csarp;
            else
                param_csarp = [];
            end
            if isfield(tmp1, 'param_radar')
                param_radar= tmp1.param_radar;
            else
                param_radar =[];
            end
            if isfield(tmp1, 'param_records')
                param_records ...
                            = tmp1.param_records;
            else
                param_records ...
                            = [];
            end
            if (jj < 10)
                save([dir_out file_in_all{ii}(1:(end - 4)) '_slice_0' num2str(jj)], '-v7.3', 'Bottom', 'Data', 'Depth', 'Elevation', 'GPS_time', 'Latitude', 'Longitude', 'Surface', 'Time', 'array_param', ...
                                                                                             'param_combine_wf_chan', 'param_csarp', 'param_radar', 'param_records')
            else
                save([dir_out file_in_all{ii}(1:(end - 4)) '_slice_' num2str(jj)], '-v7.3', 'Bottom', 'Data', 'Depth', 'Elevation', 'GPS_time', 'Latitude', 'Longitude', 'Surface', 'Time', 'array_param', ...
                                                                                            'param_combine_wf_chan', 'param_csarp', 'param_radar', 'param_records')
            end
            
        else % old format
            
            % save the sliced variable in the output directory with the slice extension
            if (jj < 10)
                save([dir_out file_in_all{ii}(1:(end - 4)) '_slice_0' num2str(jj)], '-v7.3', 'Bottom', 'Data', 'Elevation', 'GPS_time', 'Latitude', 'Longitude', 'Surface', 'Time')
            else
                save([dir_out file_in_all{ii}(1:(end - 4)) '_slice_' num2str(jj)], '-v7.3', 'Bottom', 'Data', 'Elevation', 'GPS_time', 'Latitude', 'Longitude', 'Surface', 'Time')
            end
            
        end        
        
    end
    
    tmp1_old                = tmp1; % preserve previous frame
    
end

disp(['DONE slicing frames in ' dir_in ' and saving to ' dir_out '.'])