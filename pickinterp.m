function pickinterp(dir_merge, file_merge, dir_data, file_data_all, frame_or_pk, dir_pk_orig)
% PICKINTERP Interpolates merged picked layers back into the same format as the original CReSIS data frames.
%
%   PICKINTERP(DIR_MERGE, FILE_MERGE, DIR_DATA, FILE_DATA_ALL, FRAME_OR_PK, DIR_PK_ORIG)
%   loads the merged pick file FILE_MERGE from the directory DIR_MERGE and
%   interpolates them back onto the positions of the original CReSIS data
%   files for that transect stored in DIR_DATA that follow the naming
%   sequence FILE_DATA_ALL. If FRAME_OR_PK is set to 'frame', then the
%   picks are saved to the frame. If set to 'pk', the interpolated picks
%   are saved in DIR_PK_ORIG using the naming sequence of the merged pick
%   files.
% 
% Joe MacGregor (UTIG)
% Last updated: 09/30/14

if ~any(nargin == [5 6])
    error('pickinterp:nargin', 'Incorrect number of input arguments.')
end
if ~ischar(dir_merge)
    error('pickinterp:dirmergestr', 'Picks input directory (DIR_MERGE) is not a string.')
end
if (~isempty(dir_merge) && ~exist(dir_merge, 'dir'))
    error('pickinterp:nodirmerge', 'Picks input directory (DIR_MERGE) does not exist.')
end
if ~ischar(file_merge)
    error('pickinterp:namefilestr', 'Input merged picks filename (FILE_MERGE) is not a string.')
end
if ~ischar(dir_data)
    error('pickinterp:dirdatastr', 'Data directory (DIR_DATA) is not a string.')
end
if (~isempty(dir_data) && ~exist(dir_data, 'dir'))
    error('pickinterp:nodirdata', 'Data directory (DIR_DATA) does not exist.')
end
if ~ischar(file_data_all)
    error('pickinterp:namefiledataallstr', 'Data filename (FILE_DATA_ALL) is not a string.')
end
if ~ischar(frame_or_pk)
    error('pickinterp:frameorpkstr', 'FRAME_OR_PK is not a string.')
end
if ~any(strcmp(frame_or_pk, {'frame' 'pk'}))
    error('pickinterp:frameorpk', 'FRAME_OR_PK is not a ''frame'' or ''pk''.')    
end
if (strcmp(frame_or_pk, 'frame') && (nargin ~= 6))
    error('pickinterp:frameorpkdir', 'FRAME_OR_PK=''frame'' but DIR_PK_ORIG is not provided.')
end
if strcmp(frame_or_pk, 'frame')
    if ~ischar(dir_pk_orig)
        error('pickinterp:dirpkorigstr', 'Picks output directory (DIR_PK_ORIG) is not a string.')
    end
    if (~isempty(dir_pk_orig) && ~exist(dir_pk_orig, 'dir'))
        error('pickinterp:nodirpkorig', 'Picks output directory (DIR_PK_ORIG) does not exist.')
    end
end

pk_merge                    = load([dir_merge file_merge]);
try
    pk_merge                = pk_merge.pk;
catch %#ok<*CTCH>
    disp([file_merge ' does not contain a pk structure. Try again.'])
    return
end

disp(['Loaded ' file_merge ' from ' dir_merge '...'])

% blocks to get matched up
file_data                   = dir([dir_data file_data_all '.mat']);
file_data                   = {file_data.name};
num_data                    = length(file_data);

if ~num_data
    error('pickinterp:nodata', ['No original data frames matching ' file_data_all ' found in ' dir_data ' to interpolate picks onto.'])
end

var_layer                   = {'depth' 'depth_smooth' 'elev' 'elev_gimp' 'elev_smooth' 'elev_smooth_gimp' 'ind_y' 'ind_y_smooth' 'int' 'int_smooth' 'twtt' 'twtt_ice' 'twtt_ice_smooth' 'twtt_smooth'}; % layer variables
num_var_layer               = length(var_layer);

var_pos                     = {'dist' 'elev_air' 'elev_air_gimp' 'elev_bed' 'elev_bed_gimp' 'elev_surf' 'elev_surf_gimp' 'int_bed' 'int_surf' 'lat' 'lon' 'time' 'twtt_bed' 'twtt_surf' 'x' 'y'}; % position variables

disp(['Evaluating ' num2str(num_data) ' frame(s) in ' dir_data ' for overlap...'])

for ii = 1:num_data
    
    % load current block
    frame                   = load([dir_data file_data{ii}]);
    if ~isfield(frame, 'GPS_time')
        disp([file_data{ii}(1:(end - 4)) ' does not contain Time variable. Try again.'])
        return
    end
    
    % save some memory briefly
    if strcmp(frame_or_pk, 'pk')
        frame                   = rmfield(frame, 'Data');
    end
    
    % number of traces in current frame
    num_trace_curr          = length(frame.GPS_time);
    
    % median time difference
    time_diff_med           = median(diff(frame.GPS_time));
    
    % indices in current merged file that match current data frame
    ind_curr                = NaN(1, num_trace_curr);
    for jj = 1:num_trace_curr
        if ~isempty(find((frame.GPS_time(jj) == pk_merge.time), 1))
            ind_curr(jj)    = find((frame.GPS_time(jj) == pk_merge.time), 1);
        else
            ind_time        = interp1(pk_merge.time, 1:pk_merge.num_trace_tot, frame.GPS_time(jj), 'nearest');
            if isnan(ind_time)
                continue
            end
            if (abs(pk_merge.time(ind_time) - frame.GPS_time(jj)) < time_diff_med)
                ind_curr(jj)= ind_time;
            end
        end
    end
    
    % evaluate matching between merged file and frame
    if all(isnan(ind_curr))
        disp([file_data{ii}(1:(end - 4)) ' (' num2str(ii) ' / ' num2str(num_data) ') has no picks.'])
        continue
    else
        disp([num2str(round((100 * length(find(~isnan(ind_curr)))) / num_trace_curr)) '% of ' file_data{ii}(1:(end - 4)) ' is within range of merged picks file...'])
    end
    
    % initialize pk structure
    pk                      = struct;
    pk.layer                = struct;
    
    % reassign layer values
    for jj = 1:num_var_layer
        for kk = 1:pk_merge.num_layer
            eval(['pk.layer(kk).' var_layer{jj} ' = NaN(1, num_trace_curr);'])
            if isfield(pk_merge, var_layer{jj})
                eval(['pk.layer(kk).' var_layer{jj} '(~isnan(ind_curr)) = pk_merge.' var_layer{jj} '(kk, ind_curr(~isnan(ind_curr)));'])
            end
            if ~any(strcmp(var_layer(jj), {'int' 'int_smooth'}))
                eval(['pk.layer(kk).' var_layer{jj} '(pk.layer(kk).' var_layer{jj} ' <= 0) = NaN;']) % fix occasional layer oddities
            end
        end
    end
    
    pk.ind_match            = find(sum(~isnan(pk_merge.twtt(:, ind_curr(~isnan(ind_curr)))), 2)); % layers with any non-NaN values for this chunk of transect
    pk.num_layer            = length(pk.ind_match); % number of layers for current block
    
    if ~pk.num_layer % no non-NaN layers for this frame
        disp([file_data{ii}(1:(end - 4)) ' (' num2str(ii) ' / ' num2str(num_data) ') has no non-NaN picks.'])
        continue
    end
    pk.layer                = pk.layer(pk.ind_match); % only keep layers relevant to this chunk
    
    % remove empty layers
    tmp1                    = [];
    for jj = 1:pk.num_layer
        if all(isnan(pk.layer(jj).ind_y))
            tmp1            = [tmp1 jj]; %#ok<AGROW>
        end
    end
    if ~isempty(tmp1)
        pk.layer            = pk.layer(setdiff(1:pk.num_layer, tmp1));
        pk.num_layer        = pk.num_layer - length(tmp1);
    end
    
    % reassign position variables
    for jj = 1:length(var_pos)
        eval(['pk.' var_pos{jj} ' = NaN(1, num_trace_curr);']);
        if isfield(pk_merge, var_pos{jj})
            eval(['pk.' var_pos{jj} '(~isnan(ind_curr)) = pk_merge.' var_pos{jj} '(ind_curr(~isnan(ind_curr)));'])
        end
    end
    
    pk                      = orderfields(pk);
    pk.layer                = orderfields(pk.layer);
    
    % save broken-out pk structure for this frame
    switch frame_or_pk
        case 'frame'
            save([dir_data file_data{ii}], '-append', 'pk')
            disp(['Interpolated picks to ' file_data{ii}(1:(end - 4)) ' and saved to ' dir_data'.'])
        case 'pk'
            save([dir_pk_orig file_data{ii}(1:(end - 4)) '_pk'], 'pk')
            disp(['Interpolated picks to ' file_data{ii}(1:(end - 4)) ' and saved to ' dir_pk_orig '.'])
    end
    
end