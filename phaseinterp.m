function phaseinterp(phase_type, dir_data_old, lat_std, decim, dir_block_new, file_block_new, param_phase)
% PHASEINTERP Horizontal phase gradient from complex data and interpolation onto data blocks.
%
%   PHASEINTERP(PHASE_TYPE,DIR_DATA_OLD,LAT_STD,DECIM,DIR_BLOCK_NEW,FILE_BLOCK_NEW,PARAM_PHASE)
%   calculates the horizontal phase gradient of complex radar data using
%   the method specified by PHASE_TYPE.
%   
%   PHASEINTERP loads the single channel data from DIR_DATA_OLD, determines
%   their horizontal phase gradient, decimates these data to every DECIM
%   indices and interpolates and saves this phase gradient onto the
%   horizontal positions of the focused data blocks named FILE_BLOCK_NEW in
%   DIR_BLOCK_NEW, based on the standard latitude LAT_STD.
% 
%   For PHASE_TYPE='direct', the horizontal gradient of the phase is
%   directly calculated and filtered. PARAM_PHASE must be a one-element
%   cell that contains either a numeric scalar or two-element vector, which
%   determines the filter size. If it is a scalar then it must be the
%   integer number of indices by which the horizontal phase difference is
%   filtered along columns and then rows. If it is a two-element vector,
%   then the phase difference is median-filtered along columns using the
%   first element and along rows using its second element.
%   
%   For PHASE_TYPE='Doppler', the horizontal gradient of the phase is
%   calculated from the Dopper centroid wavenumber. PARAM_PHASE must be a
%   two-element cell whose two elements are size of the FFT filter (a
%   scalar) and the wavenumber threshold (also a scalar).
%   
% Joe MacGregor (UTIG)
% Last updated: 04/09/14

if (nargin ~= 7)
    error('phaseinterp:nargin', 'Incorrect number of input arguments (must be 7).')
end
if ~ischar(phase_type)
    error('phaseinterp:phasetypechar', 'PHASE_TYPE is not a string.')
end
if ~any(strcmp(phase_type, {'direct' 'Doppler'}))
    error('phaseinterp:phasetypetype', 'PHASE_TYPE is not either ''direct'' or ''Doppler''.')
end
if ~ischar(dir_data_old)
    error('phaseinterp:dirdataoldstr', 'Old data directory (DIR_DATA_OLD) is not a string.')
end
if (~isempty(dir_data_old) && ~exist(dir_data_old, 'dir'))
    error('phaseinterp:nodirdataold', 'Old data directory (DIR_DATA_OLD) does not exist.')
end
if (~isnumeric(lat_std) || ~isscalar(lat_std))
    error('phaseinterp:latstdscalar', 'Standard parallel (LAT_STD) is not a scalar.')
end
if (abs(lat_std) > 90)
    error('phaseinterp:latstdrange', ['Absolute value of standard parallel (LAT_STD = ' num2str(abs(lat_std)) ') is larger than 90 degrees.'])
end
if (~isnumeric(decim) || ~isscalar(decim))
    error('phaseinterp:decimscalar', 'Decimation number (DECIM) is not a scalar.')
end
if any(mod(decim, 1))
    decim                   = round(decim);
    warning('phaseinterp:rounddecim', ['Decimation number (DECIM) rounded to ' num2str(decim) '.'])
end
if ~ischar(dir_block_new)
    error('phaseinterp:dirblocknewstr', 'New blocks directory (DIR_BLOCK_NEW) is not a string.')
end
if (~isempty(dir_block_new) && ~exist(dir_block_new, 'dir'))
    error('phaseinterp:nodirblocknew', 'Block directory (DIR_BLOCK_NEW) does not exist.')
end
if ~ischar(file_block_new)
    error('phaseinterp:namefileblockstr', 'New blocks filename (FILE_BLOCK_NEW) is not a string.')
end
if ~iscell(param_phase)
    error('phaseinterp:paramphasecell', 'PARAM_PHASE is not a cell.')
end
switch phase_type
    case 'direct'
        if (length(param_phase) ~= 1)
            error('phaseinterp:paramphasedirect', 'PARAM_PHASE is not one-element cell for PHASE_TYPE=''direct''.')
        end
        sz_filt             = param_phase{1};
        if (~isnumeric(sz_filt) || ~isvector(sz_filt) || (length(sz_filt) > 2))
            error('phaseinterp:szfilttype', 'Filter size (SZ_FILT) is not a scalar or two-element vector.')
        end
        if any(mod(sz_filt, 1))
            sz_filt                 = round(sz_filt);
            warning('phaseinterp:roundfilt', ['Filter size (SZ_FILT) rounded to ' num2str(sz_filt) '.'])
        end
        if any(sz_filt < 1)
            error('phaseinterp:nofilt', 'Filter size (SZ_FILT) must be greater than zero.')
        end
    case 'Doppler'
        if (length(param_phase) ~= 2)
            error('phaseinterp:paramphasedoppler', 'PARAM_PHASE is not two-element cell for PHASE_TYPE=''Doppler''.')
        end
        if ~exist('smoothn', 'file')
            error('phaseinterp:smoothn', 'Function SMOOTHN is not available within this user''s path.')
        end
        if ~exist('dctn', 'file')
            error('phaseinterp:dctn', 'Function DCTN is not available within this user''s path.')
        end
        if ~exist('idctn', 'file')
            error('phaseinterp:idctn', 'Function IDCTN is not available within this user''s path.')
        end
        [sz_filt, wavenum_max] ...
                            = deal(param_phase{1}, param_phase{2});
        if (~isnumeric(sz_filt) || ~isscalar(sz_filt))
            error('phaseinterp:szfilttype', 'FFT filter size (SZ_FILT) is not a scalar.')
        end
        if any(mod(sz_filt, 1))
            sz_filt                 = round(sz_filt);
            warning('phaseinterp:roundfilt', ['Filter size (SZ_FILT) rounded to ' num2str(sz_filt) '.'])
        end
        if any(sz_filt < 2)
            error('phaseinterp:nofilt', 'Filter size (SZ_FILT) must be greater than one.')
        end
        if (~isnumeric(wavenum_max) || ~isscalar(wavenum_max))
            error('phaseinterp:wavenummaxtype', 'Maximum wavenumber size (WAVENUM_MAX) is not a scalar.')
        end
        if ((wavenum_max <= 0) || (wavenum_max >= 1))
            error('phaseinterp:wavenummaxval', 'Maximum wavenumber size (WAVENUM_MAX) must be betwen 0 and 1.')
        end
end

% use GIMP projection if available
if (license('test', 'map_toolbox') && exist('mat/gimp_90m.mat' , 'file'))
    load mat/gimp_90m gimp_info
    gimp_avail              = true;
else
    gimp_avail              = false;
end

% preparing to load single channel data files
file_old                    = dir([dir_data_old '*.mat']);
file_old                    = {file_old.name};
num_data                    = length(file_old);
if ~num_data
    error('phaseinterp:nodatafiles', ['No data files in ' dir_data_old '.'])
end

disp(['Loading ' num2str(num_data) ' original data files in ' dir_data_old ', determining their horizontal phase gradient...'])

% preallocate key variables that will grow
[x_all, y_all, phase_diff_filt] ...
                            = deal(single([]));

if strcmp(phase_type, 'Doppler')
    wavenum_rad             = fftshift((2 * pi) .* (mod(((1 / 2) + ((0:(sz_filt - 1)) / sz_filt)), 1) - (1 / 2))); % dimensionalized wavenumbers in radians for filtering
    ind_filt                = find(abs(wavenum_rad) < (wavenum_max * 2 * pi)); % indices of wavenumbers below threshold
end

% loop through data_old/ files
for ii = 1:num_data
    
    disp([num2str(ii) ' / ' num2str(num_data) '...'])
    
    % load data extract key parts
    try
        data                = load([dir_data_old file_old{ii}]);
    catch
        disp('Could not open file for unknown reason...')
        continue
    end
    [lat, lon, time]        = deal(data.hdr.lat, data.hdr.lon, data.hdr.gps_time);
    if ((ii == 1) || ~exist('twtt_old', 'var'))
       twtt_old             = data.hdr.wfs(2).time;
    elseif (data.hdr.wfs(2).time(1) < twtt_old(1))
       twtt_old             = data.hdr.wfs(2).time;
    end
    data                    = data.data;
    
    % ignore overlapping data
    if ((ii > 1) && exist('time_old', 'var'))
        if any(time < time_old(end))
            ind_tmp         = find(time > time_old(end));
            if ~isempty(ind_tmp)
                [lat, lon, data] ...
                            = deal(lat(ind_tmp), lon(ind_tmp), data(:, ind_tmp));
            else
                disp(['Skipping ' file_old{ii} ' because it overlaps completely in time with the previous data file (' file_old{ii - 1} ').'])
                continue % no points in this data file that don't overlap with previous file (whoa)
            end
        end
    end
    
    [num_sample, num_trace] = size(data); % current size of data matrix
    if (num_trace < 2)
        disp(['Not enough traces ' num2str(num_trace) ' to be worthwhile...'])
        continue
    end
    
    % calculate x/y positions
    if (license('checkout', 'map_toolbox') && (lat(1) > 0)) % in Greenland and can do better x/y coordinates
        if gimp_avail
            [x_curr, y_curr]= projfwd(gimp_info, lat, lon);
        else            
            wgs84           = almanac('earth', 'wgs84', 'meters');
            ps_struct       = defaultm('ups');
            [ps_struct.geoid, ps_struct.mapparallels, ps_struct.falsenorthing, ps_struct.falseeasting, ps_struct.origin] ...
                            = deal(wgs84, lat_std, 0, 0, [90 -45 0]);
            [x_curr, y_curr]= mfwdtran(ps_struct, lat, lon);
        end
    else
        [x_curr, y_curr]    = ll2ps(lat, lon, lat_std); % convert lat/lon to polar stereographic x/y (better for Antarctica then Greenland)
    end
    [x_curr, y_curr]        = deal((1e-3 .* (x_curr + ([diff(x_curr) diff(x_curr((end - 1):end))] ./ 2))), (1e-3 .* (y_curr + ([diff(y_curr) diff(y_curr((end - 1):end))] ./ 2)))); % account for x positions being at difference (not on the dot)
    
    % prepare filter and decimation vectors
    if (decim > 1)
        decim_vec           = (1 + ceil(decim / 2)):decim:(num_trace - ceil(decim / 2));
    else
        decim_vec           = 1:num_trace;
    end
    
    switch phase_type
        
        case 'direct'
            
            phase_diff      = diff(unwrap(angle(data), [], 2), 1, 2); % horizontal difference of horizontally unwrapped phase
            phase_diff(phase_diff > pi) ...
                            = pi;
            phase_diff(phase_diff < -pi) ...
                            = -pi;
            switch length(sz_filt)
                case 1
                    phase_diff_filt_curr ...
                            = medfilt2(medfilt2(phase_diff, [(sz_filt + 1) 1], 'symmetric'), [1 (sz_filt + 1)], 'symmetric'); % median filtering
                case 2
                    phase_diff_filt_curr ...
                            = medfilt2(medfilt2(phase_diff, [(sz_filt(1) + 1) 1], 'symmetric'), [1 (sz_filt(2) + 1)], 'symmetric');
            end
            phase_diff_filt_curr ...
                            = phase_diff_filt_curr(:, decim_vec(decim_vec <= size(phase_diff_filt_curr, 2)));
            
        case 'Doppler'
            
            filt_vec        = (1 + ceil(sz_filt / 2)):ceil(sz_filt / 2):(num_trace - ceil(sz_filt / 2));
            num_filt        = length(filt_vec);
            
            if (num_filt <= 1) % otherwise a bit pointless
                continue
            end
            phase_diff_filt_curr ...
                            = zeros(num_sample, num_filt, 'single');
            
            % calculate fft, filter/narrow and find peak
            for jj = 1:num_filt
                curr_vec    = (filt_vec(jj) - floor(sz_filt / 2)):(filt_vec(jj) + floor(sz_filt / 2) - 1);
                data_fft    = fftshift(fft(data(:, curr_vec), [], 2), 2);
                [int_fft_max, ind_fft_max] ...
                            = max(abs(data_fft(:, ind_filt)), [], 2);
                tmp1        = ind_fft_max;
                tmp1((int_fft_max ./ median(abs(data_fft(:, ind_filt)), 2)) < 5) ...
                            = NaN; % ignore data with poor S/N
                phase_diff_filt_curr(:, jj) ...
                            = interp1(1:sz_filt, wavenum_rad, smoothn((tmp1 + ind_filt(1) - 1), 1), 'linear'); % adjust peak indices, smooth them, and get their wavenumber values
            end
            
            % interpolate phase difference onto current x/y
            [tmp1, tmp2]    = meshgrid(filt_vec, 1:num_sample);
            [tmp3, tmp4]    = meshgrid(decim_vec, 1:num_sample);
            
            warning('off', 'MATLAB:interp2:NaNstrip')
            try
                phase_diff_filt_curr ...
                            = interp2(tmp1, tmp2, phase_diff_filt_curr, tmp3, tmp4, 'spline'); % interpolate horizontal phase difference to decimation vector
            catch
                disp(['Skipping ' num2str(ii) ' due to a 2-D interpolation error...'])
                continue
            end
            warning('on', 'MATLAB:interp2:NaNstrip')
            phase_diff_filt_curr(phase_diff_filt_curr > pi) ...
                            = pi;
            phase_diff_filt_curr(phase_diff_filt_curr < -pi) ...
                            = -pi;
    end
        
    % decimate x/y and preserve x/y/phase_diff_filt
    [x_all, y_all]          = deal([x_all; single(x_curr(decim_vec))'], [y_all; single(y_curr(decim_vec))']);
    phase_diff_filt         = [phase_diff_filt phase_diff_filt_curr]; %#ok<AGROW>
    
    clear data phase_diff_filt_curr x_curr y_curr dist data_fft tmp*
    
    time_old                = time;
end

if isempty(phase_diff_filt)
    disp('No phase to interpolate...')
    return
end

dist_all                    = cumsum([0; sqrt((diff(x_all) .^ 2) + (diff(y_all) .^ 2))]); % distance vector for concatened/decimated single channel data

% prepare to load blocks to get matched up with their phase
file_new                    = dir([dir_block_new file_block_new '.mat']);
file_new                    = {file_new.name};
num_block                   = length(file_new);
if ~num_block
    error('phaseinterp:nofiles', ['No files matching ' file_block_new '.mat in ' dir_block_new '.'])
end

disp(['Loading ' num2str(num_block) ' focused blocks in ' dir_block_new ' and interpolating horizontal phase gradient onto them...'])

for ii = 1:num_block
    
    disp([num2str(ii) ' / ' num2str(num_block) '...'])
    
    % load current block
    block                   = load([dir_block_new file_new{ii}]);
    try
        block               = block.block;
    catch %#ok<*CTCH>
        disp([file_new{ii}(1:(end - 4)) ' does not contain a block structure. Try again.'])
        return
    end
    
    if isfield(block, 'phase_diff_filt')
        block               = rmfield(block, 'phase_diff_filt');
    end
    
    if (all(isfield(block, {'x_gimp', 'y_gimp'})) && gimp_avail)
        [x_block, y_block]  = deal(single(block.x_gimp), single(block.y_gimp)); % x/y positions of blocks        
    elseif (~all(isfield(block, 'x_gimp', 'y_gimp')) && gimp_avail)
        disp('Original data in GIMP projection but block does not have not have GIMP projection so skipping...')
        continue
    else
        [x_block, y_block]  = deal(single(block.x), single(block.y)); % x/y positions of blocks
    end
    
    % range within old single channel data to evaluate differenced phase
    dist_range              = [(block.dist(1) - (0.5 * diff(block.dist([1 end])))) (block.dist(end) + (0.5 * diff(block.dist([1 end]))))];
    
    try
        ind_range           = [find((dist_all > dist_range(1)), 1) find((dist_all < dist_range(2)), 1, 'last')];
    catch
        disp([file_new{ii}(1:(end - 4)) ' (' num2str(ii) ' / ' num2str(num_block) ') is being skipped because its positions would not be matched across the entire block.'])
        continue
    end
    if (any(isnan(ind_range)) || (length(ind_range) ~= 2))
        disp([file_new{ii}(1:(end - 4)) ' (' num2str(ii) ' / ' num2str(num_block) ') is being skipped because its positions would not be matched across the entire block.'])
        continue
    end
    
    % find current indices in the merged file that connect with current block
    [dist_curr, ind_curr]   = min(sqrt(((x_all(ind_range(1):ind_range(2), ones(1, block.num_trace)) - x_block(ones(1, (diff(ind_range) + 1)), :)) .^ 2) + ((y_all(ind_range(1):ind_range(2), ones(1, block.num_trace)) - y_block(ones(1, (diff(ind_range) + 1)), :)) .^ 2)));
    ind_curr(dist_curr > (10 * mean(diff(block.dist)))) ...
                            = NaN;
    if any(isnan(ind_curr)) % need phase across entire transect otherwise using it to predict layers in pickgui becomes too complex
        disp([file_new{ii}(1:(end - 4)) ' (' num2str(ii) ' / ' num2str(num_block) ') is being skipped because its positions would not be matched across the entire block.'])
        continue
    end
    ind_curr                = ind_curr + ind_range(1) - 1;
    
    % map phase gradient onto current block
    try
        block.phase_diff_filt ...
                            = phase_diff_filt(interp1(twtt_old, 1:length(twtt_old), block.twtt, 'nearest', 'extrap'), ind_curr);
    catch
        disp('Block overlap failed.')
        continue
    end
    block                   = orderfields(block); %#ok<NASGU>
    
    % save block with phase added in
    save([dir_block_new file_new{ii}], 'block')
    disp(['Filtered and interpolated horizontal phase gradient for ' file_new{ii}(1:(end - 4)) ' and saved in ' dir_block_new '.'])
end