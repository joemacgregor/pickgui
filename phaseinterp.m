function phaseinterp(dir_data_old, sz_filt, wavenum_max, freq, lat_std, decim, dir_block_new, file_block_new)
% PHASEINTERP Determines horizontal phase gradient from complex data and
% interpolates onto data blocks.
%
%   PHASEINTERP(DIR_DATA_OLD,SZ_FILT,WAVENUM_MAX,LAT_STD,DECIM,DIR_BLOCK_NEW,FILE_BLOCK_NEW)
%   loads the single channel data in DIR_DATA_OLD, determine the Doppler
%   centroid frequency, aka their horizontal phase gradient, using a FFT
%   filter of size SZ_FILT and wavenumber threshold of WAVENUM_MAX,
%   interpolates that phase gradient onto the horizontal positions of the
%   focused blocks in DIR_BLOCK_NEW based on the standard latitude LAT_STD
%   after decimation by DECIM indices, and saves the focused blocks back
%   into DIR_BLOCK_NEW.
% 
% Joe MacGregor (UTIG)
% Last updated: 12/14/12

% clear
% dir_data_old                = '~/Desktop/test7/data_old/';
% % dir_data_old                = '2011_p3/data_old/20110502_01/';
% sz_filt                     = 512;
% wavenum_max                 = 0.1;
% freq                        = 195e6;
% lat_std                     = 70;
% % decim                       = 1;
% decim                       = 10;
% dir_block_new               = '~/Desktop/test7/block/';
% % dir_block_new               = '2011_p3/block/20110502_01/';
% file_block_new              = '*';

if ~exist('smoothn', 'file')
    error('phaseinterp:smoothn', 'Function SMOOTHN is not available within this user''s path.')
end
if ~exist('dctn', 'file')
    error('phaseinterp:dctn', 'Function DCTN is not available within this user''s path.')
end
if ~exist('idctn', 'file')
    error('phaseinterp:idctn', 'Function IDCTN is not available within this user''s path.')
end
if (nargin < 8)
    error('phaseinterp:nargin', 'Not enough input arguments.')
end
if ~ischar(dir_data_old)
    error('phaseinterp:dirdataoldstr', 'Old data directory (DIR_DATA_OLD) is not a string.')
end
if (~isempty(dir_data_old) && ~exist(dir_data_old, 'dir'))
    error('phaseinterp:nodirdataold', 'Old data directory (DIR_DATA_OLD) does not exist.')
end
if (~isnumeric(sz_filt) || (length(sz_filt) > 2))
    error('phaseinterp:szfilttype', 'Filter size (SZ_FILT) is not a scalar or two-element vector.')
end
if any(mod(sz_filt, 1))
    sz_filt                 = round(sz_filt);
    warning('phaseinterp:roundfilt', ['Filter size (SZ_FILT) rounded to ' num2str(sz_filt) '.'])
end
if any(sz_filt < 1)
    error('phaseinterp:nofilt', 'Filter size (SZ_FILT) must be greater than zero.')
end
if (~isnumeric(wavenum_max) || (length(wavenum_max) > 1))
    error('phaseinterp:wavenummaxtype', 'Maximum wavenumber size (WAVENUM_MAX) is not a scalar.')
end
if ((wavenum_max <= 0) || (wavenum_max >= 1))
    error('phaseinterp:wavenummaxval', 'Maximum wavenumber size (WAVENUM_MAX) must be betwen 0 and 1.')
end
if (~isnumeric(freq) || (length(freq) > 1))
    error('phaseinterp:freqtype', 'Center frequency (FREQ) is not a scalar.')
end
if (~isnumeric(lat_std) || (length(lat_std) > 1))
    error('phaseinterp:latstdscalar', 'Standard parallel (LAT_STD) is not a scalar.')
end
if (abs(lat_std) > 90)
    error('phaseinterp:latstdrange', ['Absolute value of standard parallel (LAT_STD = ' num2str(abs(lat_std)) ') is larger than 90 degrees.'])
end
if (~isnumeric(decim) || (length(decim) > 1))
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

% preparing to load single channel data files
file_old                    = dir([dir_data_old '*.mat']);
file_old                    = {file_old.name};
num_data                    = length(file_old);
if ~num_data
    error('phaseinterp:nodatafiles', ['No data files in ' dir_data_old '.'])
end

disp(['Loading ' num2str(num_data) ' original data files in ' dir_data_old ', determining their horizontal phase gradient...'])

% preallocate key variables that will grow
[x_old, y_old, phase_diff_filt_old] ...%, specularity_old, phase_old, phase_diff_alt] ...
                            = deal(single([]));

wavenum_rad                 = fftshift((2 * pi) .* (mod(((1 / 2) + ((0:(sz_filt - 1)) / sz_filt)), 1) - (1 / 2))); % dimensionalized wavenumbers in radians for filtering
ind_filt                    = find(abs(wavenum_rad) < (wavenum_max * 2 * pi)); % indices of wavenumbers below threshold
% speed_vacuum                = 299792458; % m/s
% permitt_ice                 = 3.15; % CReSIS standard value
% wavelength_air              = speed_vacuum / freq;
% wavelength_ice              = speed_vacuum / (freq * sqrt(permitt_ice));

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
       twtt_old             = data.hdr.wfs(1).time;
    elseif (data.hdr.wfs(1).time(1) < twtt_old(1))
       twtt_old             = data.hdr.wfs(1).time;
    end
%     if ((ii == 1) || ~exist('twtt_old', 'var'))
%         twtt_old             = data.hdr.wfs(end).time;
%     elseif (data.hdr.wfs(end).time(1) < twtt_old(1))
%         twtt_old             = data.hdr.wfs(end).time;
%     end
%     twtt_surf               = data.hdr.surface;
%     ind_surf                = interp1(twtt_old, 1:length(twtt_old), twtt_surf, 'nearest', 'extrap');
%     ind_twtt_0              = interp1(twtt_old, 1:length(twtt_old), 0, 'nearest', 'extrap');
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
    if (license('test', 'map_toolbox') && (lat(1) > 0)) % in Greenland and can do better x/y coordinates
        wgs84               = almanac('earth', 'wgs84', 'meters');
        ps_struct           = defaultm('ups');
        [ps_struct.geoid, ps_struct.mapparallels, ps_struct.falsenorthing, ps_struct.falseeasting, ps_struct.origin] ...
                            = deal(wgs84, lat_std, 0, 0, [90 -45 0]);
        [x, y]              = mfwdtran(ps_struct, lat, lon);
    else
        [x, y]              = ll2ps(lat, lon, lat_std); % convert lat/lon to polar stereographic x/y (better for Antarctica then Greenland)
    end
    [x, y]                  = deal((1e-3 .* (x + ([diff(x) diff(x((end - 1):end))] ./ 2))), (1e-3 .* (y + ([diff(y) diff(y((end - 1):end))] ./ 2)))); % account for x positions being at difference (not on the dot)
%     dist_diff_mean          = mean(diff(1e3 .* cumsum([0 sqrt((diff(x) .^ 2) + (diff(y) .^ 2))])));
%     wavenum                 = wavenum_rad ./ (2 * pi * dist_diff_mean);
    
    % prepare filter and decimation vectors
    filt_vec                = (1 + ceil(sz_filt / 2)):ceil(sz_filt / 2):(num_trace - ceil(sz_filt / 2));
    if (decim > 1)
        decim_vec           = (1 + ceil(decim / 2)):decim:(num_trace - ceil(decim / 2));
    else
        decim_vec           = 1:num_trace;
    end
    num_filt                = length(filt_vec);
    
    if (num_filt > 1) % otherwise a bit pointless
        
        phase_diff_filt     = zeros(num_sample, num_filt, 'single');
%         [phase_diff_filt, specularity] ...
%                             = deal(zeros(num_sample, num_filt, 'single'));
        
        % calculate fft, filter/narrow and find peak
        for jj = 1:num_filt
            
            curr_vec        = (filt_vec(jj) - floor(sz_filt / 2)):(filt_vec(jj) + floor(sz_filt / 2) - 1);
%             ind_surf_curr   = round(mean(ind_surf(curr_vec)));
            data_fft        = fftshift(fft(data(:, curr_vec), [], 2), 2);
            [int_fft_max, ind_fft_max] ...
                            = max(abs(data_fft(:, ind_filt)), [], 2);
            tmp1            = ind_fft_max;
            tmp1((int_fft_max ./ median(abs(data_fft(:, ind_filt)), 2)) < 5) ...
                            = NaN; % ignore data with poor S/N
            
            phase_diff_filt(:, jj) ...
                            = interp1(1:sz_filt, wavenum_rad, smoothn((tmp1 + ind_filt(1) - 1), 1), 'linear'); % adjust peak indices, smooth them, and get their wavenumber values
            
%             frac_air        = (round(mean(ind_surf(curr_vec))) - ind_twtt_0) ./ ((1:num_sample)' - ind_twtt_0);
%             frac_air(1:(ind_surf_curr - 1)) ...
%                             = 1;
%             
%             twtt_moveout    = zeros(num_sample, sz_filt); % moveout shape
%             for kk = (ind_twtt_0 + 1):num_sample
%                 twtt_moveout(kk, :) ...
%                             = twtt_old(kk) ./ sqrt(1 - (((0.5 * ((wavelength_air * frac_air(kk)) + ((1 - frac_air(kk)) * wavelength_ice))) .* (wavenum - wavenum(ind_fft_max(kk) + ind_filt(1) - 1))) .^ 2));
%             end
%             twtt_moveout(logical(imag(twtt_moveout))) ...
%                             = 0;
%             ind_diff        = interp1(twtt_old, 1:num_sample, twtt_moveout, 'nearest', 'extrap');
%             
%             for kk = ind_surf_curr:num_sample
%                 specularity(kk, jj) ...
%                             = int_fft_max(kk) ./ sum(abs(data_fft(sub2ind([num_sample sz_filt], ind_diff(kk, :), 1:sz_filt))));
%             end
            
        end
        
%         phase_diff          = diff(unwrap(angle(data), [], 2), 1, 2);
%         phase_diff(phase_diff > pi) ...
%                             = pi;
%         phase_diff(phase_diff < -pi) ...
%                             = -pi;
%         phase_diff_filt_alt = medfilt2(medfilt2(phase_diff, [(25 + 1) 1], 'symmetric'), [1 (25 + 1)], 'symmetric'); % old way using median filtering
        
        % interpolate phase difference and specularity onto current x/y
        [tmp1, tmp2]        = meshgrid(filt_vec, 1:num_sample);
        [tmp3, tmp4]        = meshgrid(decim_vec, 1:num_sample);
        
        warning('off', 'MATLAB:interp2:NaNstrip')
        try
            phase_diff_filt = interp2(tmp1, tmp2, phase_diff_filt, tmp3, tmp4, 'spline'); % interpolate horizontal phase difference to decimation vector
        catch
            disp(['Skipping ' num2str(ii) ' due to a 2-D interpolation error...'])
            continue
        end
        warning('on', 'MATLAB:interp2:NaNstrip')
        phase_diff_filt(phase_diff_filt > pi) ...
                            = pi;
        phase_diff_filt(phase_diff_filt < -pi) ...
                            = -pi;
%         specularity         = interp2(tmp1, tmp2, specularity, tmp3, tmp4, 'spline');
%         specularity(specularity < 0) ...
%                             = 0;
        
        % decimate x/y and preserve x/y/phase_diff_filt
        [x_old, y_old]      = deal([x_old; single(x(decim_vec))'], [y_old; single(y(decim_vec))']);
        phase_diff_filt_old = [phase_diff_filt_old phase_diff_filt]; %#ok<AGROW>
%         specularity_old     = [specularity_old specularity]; %#ok<AGROW>
%         phase_old           = [phase_old angle(data)];
%         phase_diff_alt      = [phase_diff_alt phase_diff_filt_alt];
        
    end
    clear data phase_diff_filt specularity x y dist data_fft tmp*
    
    time_old                = time;
    
end

clear time_old

if isempty(phase_diff_filt_old)
    disp('No phase to interpolate...')
    return
end

dist_old                    = cumsum([0; sqrt((diff(x_old) .^ 2) + (diff(y_old) .^ 2))]); % distance vector for concatened/decimated single channel data

% prepare to load blocks to get matched up with their phase
file_new                    = dir([dir_block_new file_block_new '.mat']);
file_new                    = {file_new.name};
num_block                   = length(file_new);
if ~num_block
    error('phaseinterp:nofiles', ['No files matching ' file_block_new '.mat in ' dir_block_new '.'])
end

disp(['Loading ' num2str(num_block) ' focused blocks in ' dir_block_new ' and interpolating horizonal phase gradient and specularity onto them...'])

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
        block                = rmfield(block, 'phase_diff_filt');
    end
    
    [x_b, y_b]              = deal(single(block.x), single(block.y)); % x/y positions of blocks
    
    % range within old single channel data to evaluate differenced phase
    
    dist_range              = [(block.dist(1) - (0.5 * diff(block.dist([1 end])))) (block.dist(end) + (0.5 * diff(block.dist([1 end]))))];
    
%     ind_range               = NaN(1, 2);
    try
        ind_range           = [find((dist_old > dist_range(1)), 1) find((dist_old < dist_range(2)), 1, 'last')];
    catch
        disp([file_new{ii}(1:(end - 4)) ' (' num2str(ii) ' / ' num2str(num_block) ') is being skipped because its positions would not be matched across the entire block.'])
        continue
    end
    if (any(isnan(ind_range)) || (length(ind_range) ~= 2))
        disp([file_new{ii}(1:(end - 4)) ' (' num2str(ii) ' / ' num2str(num_block) ') is being skipped because its positions would not be matched across the entire block.'])
        continue
    end
    
%     ind_range               = [1 112834];

    % find current indices in the merged file that connect with current block
    [dist_curr, ind_curr]   = min(sqrt(((x_old(ind_range(1):ind_range(2), ones(1, block.num_trace)) - x_b(ones(1, (diff(ind_range) + 1)), :)) .^ 2) + ...
                                       ((y_old(ind_range(1):ind_range(2), ones(1, block.num_trace)) - y_b(ones(1, (diff(ind_range) + 1)), :)) .^ 2)));
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
                            = phase_diff_filt_old(interp1(twtt_old, 1:length(twtt_old), block.twtt, 'nearest', 'extrap'), ind_curr);
%         block.specularity   = specularity_old(interp1(twtt_old, 1:length(twtt_old), block.twtt, 'nearest', 'extrap'), ind_curr;
%         tmp1                = phase_old(interp1(twtt_old, 1:length(twtt_old), block.twtt, 'nearest', 'extrap'), ind_curr);
%         tmp2                = phase_diff_alt(interp1(twtt_old, 1:length(twtt_old), block.twtt, 'nearest', 'extrap'), ind_curr);
    catch
        disp('Block overlap failed.')
        continue
    end
    block                   = orderfields(block); %#ok<NASGU>
    
    % save block with phase added in
%     save([dir_block_new file_new{ii}], 'block')
    disp(['Filtered and interpolated horizontal phase gradient and calculated specularity to ' file_new{ii}(1:(end - 4)) ' and saved in ' dir_block_new '.'])
    
end

% %%
% 
% figure('position',[100 100 1200 800])
% imagesc(dist_old(1:10313), (1e6 .* twtt_old), phase_old(:, 1:10313), [-pi pi])
% colormap(jet)
% set(gca, 'fontsize', 20)
% xlabel('Distance (km)')
% ylabel('Traveltime ({\mu}s)')
% text(5.02, 1.8, '(rad)', 'fontsize', 20)
% colorbar('fontsize', 20)
% box on
% axis([0 5 2.7524 31.579])
% 
% %%
% figure('position',[100 100 1200 800])
% imagesc(dist_old(1:5:end), (1e6 .* twtt_old), phase_diff_alt(:, 1:5:end), [-0.25 0.25])
% colormap(jet)
% set(gca, 'fontsize', 20)
% xlabel('Distance (km)')
% ylabel('Traveltime ({\mu}s)')
% text(54, 1.6, '(rad trace^{-1})', 'fontsize', 20)
% colorbar('fontsize', 20)
% box on
% ylim([2.7524 31.579])
% 
% 
% %%
% 
% figure('position', [100 100 1200 800])
% imagesc(wavenum_rad, (1e6 .* twtt_old), (10 .* log10(abs(data_fft))), [-70 -30])
% colormap(jet)
% hold on
% vl = vline(0, 'w--');
% set(vl, 'linewidth', 3)
% set(gca, 'fontsize', 20)
% xlabel('Wavenumber (rad)')
% ylabel('Traveltime ({\mu}s)')
% text(2.05, 1.6, '(dB)', 'fontsize', 20)
% colorbar('fontsize', 20)
% box on
% axis([-2 2 2.7524 31.579])