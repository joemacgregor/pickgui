function pickinterp(dir_merge, file_merge, dir_data, file_data_all, dir_pk_orig)
% PICKINTERP Interpolates merged picked layers back onto the original CReSIS data files.
%
%   PICKINTERP(DIR_MERGE, FILE_MERGE, DIR_DATA, FILE_DATA, DIR_PK_ORIG)
%   loads the merged pick file FILE_MERGE from the directory DIR_MERGE and
%   interpolates them back onto the positions of the original CReSIS data
%   files for that transect stored in DIR_DATA that follow the naming
%   sequence FILE_DATA. The interpolated picks are saved in DIR_PK_ORIG
%   using the naming sequence of the merged pick files.
% 
% Joe MacGregor (UTIG)
% Last updated: 11/26/12

% clear
% nargin                      = 5;
% dir_merge                   = '~/Desktop/test2/merge/';
% file_merge                  = '20110502_02_pk_merge';
% dir_data                    = '~/Desktop/test2/data/';
% file_data_all               = '*';
% dir_pk_orig                 = '~/Desktop/test2/pk_orig/';

if (nargin ~= 5)
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
    error('pickinterp:namefiledatastr', 'Data filename (FILE_DATA) is not a string.')
end
if ~ischar(dir_pk_orig)
    error('pickinterp:dirpkorigstr', 'Picks output directory (DIR_PK_ORIG) is not a string.')
end
if (~isempty(dir_pk_orig) && ~exist(dir_pk_orig, 'dir'))
    error('pickinterp:nodirpkorig', 'Picks output directory (DIR_PK_ORIG) does not exist.')
end

pk_merge                    = load([dir_merge file_merge]);
try
    pk_merge                = pk_merge.pk;
    disp(['Loaded ' file_merge ' from ' dir_merge '...'])
catch %#ok<*CTCH>
    disp([file_merge ' does not contain a pk structure. Try again.'])
    return
end

% blocks to get matched up
file_data                   = dir([dir_data file_data_all '.mat']);
file_data                   = {file_data.name};
num_data                    = length(file_data);

if ~num_data
    error('pickinterp:nodata', ['No original data blocks matching ' file_data_all ' found in ' dir_data ' to interpolate picks onto.'])
end

var_layer                   = {'ind_y' 'twtt' 'twtt_ice' 'int' 'ind_y_smooth' 'ind_y_flat_smooth' 'twtt_smooth' 'twtt_ice_smooth' 'int_smooth' 'depth' 'depth_smooth' 'elev' 'elev_smooth'}; % layer variables
var_pos                     = {'dist' 'lat' 'lon' 'x' 'y' 'elev_air' 'twtt_surf' 'elev_surf' 'twtt_bed' 'elev_bed'}; % position variables
var_pos_orig                = {'' 'Latitude' 'Longitude' '' '' 'Elevation' 'Surface' '' 'Bottom' ''};
num_var_layer               = length(var_layer);

% check to see if surface and bed picks are available
if (~isfield(pk_merge, 'twtt_surf') && ~isfield(pk_merge, 'twtt_bed'))
    ind_var_pos             = 1:6;
elseif (isfield(pk_merge, 'twtt_surf') && ~isfield(pk_merge, 'twtt_bed'))
    ind_var_pos             = 1:8;
elseif (~isfield(pk_merge, 'twtt_surf') && isfield(pk_merge, 'twtt_bed'))
    ind_var_pos             = [1:6 9 10];
else
    ind_var_pos             = 1:length(var_pos);
end

disp(['Evaluating ' num2str(num_data) ' block(s) in ' dir_data ' for overlap...'])
                        
for ii = 1:num_data
    
    % load current block
    data                    = load([dir_data file_data{ii}]);
    if ~isfield(data, 'Data')
        disp([file_data{ii}(1:(end - 4)) ' does not contain Data variable. Try again.'])
        return
    end
    data                    = rmfield(data, 'Data'); % save some memory
        
    num_trace_curr          = length(data.Latitude);
    
    ind_curr                = NaN(1, num_trace_curr);
    for jj = 1:num_trace_curr
        if ~isempty(find(((data.Latitude(jj) == pk_merge.lat) & (data.Longitude(jj) == pk_merge.lon)), 1))
            ind_curr(jj)    = find(((data.Latitude(jj) == pk_merge.lat) & (data.Longitude(jj) == pk_merge.lon)), 1);
        end
    end
    
    if all(isnan(ind_curr))
        disp([file_data{ii}(1:(end - 4)) ' (' num2str(ii) ' / ' num2str(num_data) ') has no picks.'])
        continue
    else
        disp([num2str(round((100 * length(find(~isnan(ind_curr)))) / num_trace_curr)) '% of ' file_data{ii}(1:(end - 4)) ' is within range of merged picks file...'])
    end
    
    pk                      = struct; % initialize pk structure
    pk.layer                = struct;
    
    % reassign layer, position and miscellaneous variables to current range
    curr_block              = find((find(~isnan(ind_curr), 1) >= pk_merge.ind_trace_start), 1, 'last'):find((find(~isnan(ind_curr), 1, 'last') > pk_merge.ind_trace_start), 1, 'last');
    
    [pk.ind_layer_keep, pk.ind_x_start, pk.ind_y_man, pk.ind_y_phase_max, pk.num_layer_keep, pk.num_man, pk.num_phase] ...
                            = deal(0);
    [pk.file_block, pk.file_in, pk.freq, pk.ind_trim_start, pk.length_smooth, pk.num_ind_mean, pk.num_sample, pk.num_trace, pk.num_win, pk.twtt_min_ref, pk.twtt_match, pk.twtt_max_ref] ...
                            = deal(file_data{ii}, file_data{ii}, pk_merge.freq(curr_block(1)), min(pk_merge.ind_trim_start(curr_block)), max(pk_merge.length_smooth(curr_block)), 5, length(data.Time), ...
                                   num_trace_curr, max(pk_merge.num_win(curr_block)), min(pk_merge.twtt_min_ref(curr_block)), max(pk_merge.twtt_match(curr_block)), max(pk_merge.twtt_max_ref(curr_block)));
    pk.ind_x_mean           = ceil((pk.num_ind_mean / 2) + 1):pk.num_ind_mean:(pk.num_trace - ceil(pk.num_ind_mean / 2));
    
    for jj = 1:num_var_layer
        for kk = 1:pk_merge.num_layer
            eval(['pk.layer(kk).' var_layer{jj} ' = NaN(1, num_trace_curr);'])
            eval(['pk.layer(kk).' var_layer{jj} '(~isnan(ind_curr)) = pk_merge.' var_layer{jj} '(kk, ind_curr(~isnan(ind_curr)));'])
            eval(['pk.layer(kk).' var_layer{jj} '(pk.layer(kk).' var_layer{jj} ' <= 0) = NaN;']) % fix layer oddities
        end
    end
    pk.ind_match            = find(sum(~isnan(pk_merge.twtt(:, ind_curr(~isnan(ind_curr)))), 2)); % layers with any non-NaN values for this chunk of transect
    pk.num_layer            = length(pk.ind_match); % number of layers for current block
    if ~pk.num_layer % no non-NaN layers for this
        disp([file_data{ii}(1:(end - 4)) ' (' num2str(ii) ' / ' num2str(num_data) ') has no non-NaN picks.'])
        continue
    end
    pk.layer                = pk.layer(pk.ind_match); % only keep layers relevant to this chunk
    
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
    
    for jj = ind_var_pos
        if ~isempty(var_pos_orig{jj})
            eval(['pk.' var_pos{jj} ' = data.' var_pos_orig{jj} ';'])
        else
            eval(['pk.' var_pos{jj} ' = NaN(1, num_trace_curr);']);
        end
    end
    
    pk                      = orderfields(pk);
    pk.layer                = orderfields(pk.layer);
    
    % save broken-out pk structure for this block
    save([dir_pk_orig file_data{ii}(1:(end - 4)) '_pk'], 'pk')
    disp(['Interpolated picks to ' file_data{ii}(1:(end - 4)) ' and saved to ' dir_pk_orig '.'])
    
end