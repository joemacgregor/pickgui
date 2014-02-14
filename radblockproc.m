function radblockproc(dir_in, file_in, file_block, num_file_block, num_overlap, lat_std, dir_out)
% RADBLOCKPROC Consolidate and pre-process CReSIS ice-penetrating radar data files into manageable blocks.
% 
%   RADBLOCKPROC(DIR_IN,FILE_IN,FILE_BLOCK,NUM_FILE_BLOCK,NUM_OVERLAP,LAT_STD,DIR_OUT)
%   pre-processes raw radar data stored in DIR_IN whose names are of the
%   form FILE_IN, which can contain wildcards (*), and saves them in
%   DIR_OUT with the form of FILE_BLOCK. The ".mat" extension is implicit.
%   NUM_FILE_BLOCK is the maximum number of files per block, which is used
%   for all blocks except the last one, which may have contain fewer files.
%   NUM_OVERLAP is the number of overlapping files per block. LAT_STD is
%   the standard parallel in a polar stereographic projection in decimal
%   degrees (N: positive; S: negative).
% 
%   If the Mapping Toolbox is licensed and available and the radar transect
%   is north of the equator, polar stereographic coordinates will be
%   calculated assuming a standard meridian of 45W. Otherwise, the function
%   LL2PS must be available within the user's path.
% 
% Joe MacGregor (UTIG), Mark Fahnestock (UAF-GI)
% Last updated: 01/19/13

if (nargin ~= 7)
    error('radblockproc:inputs', ['Number input arguments (' numstr(nargin) ') should be 7.'])
end
if ~ischar(dir_in)
    error('radblockproc:dirinstr', 'Data input directory (DIR_IN) is not a string.')
end
if (~isempty(dir_in) && ~exist(dir_in, 'dir'))
    error('radblockproc:nodirin', 'Data input directory (DIR_IN) does not exist.')
end
if ~ischar(file_in)
    error('radblockproc:namefilestr', 'Data filename (FILE_IN) is not a string.')
end
if ~ischar(file_block)
    error('radblockproc:nameblockstr', 'Ouput block filename (FILE_BLOCK) is not a string.')
end
if (~isnumeric(num_file_block) || ~isscalar(num_file_block))
    error('radblockproc:numfileblockposscalar', 'Number of files per block (NUM_FILE_BLOCK) is not a scalar.')
end
if mod(num_file_block, 1)
    num_file_block          = round(num_file_block);
    warning('radblockproc:roundblock', ['Number of files per block (NUM_FILE_BLOCK) rounded to ' num2str(num_file_block) '.'])
end
if (num_file_block < 1)
    error('radblockproc:noblock', 'Number of files per block (NUM_FILE_BLOCK) must be greater than zero.');
end
if (~isnumeric(num_overlap) || ~isscalar(num_overlap))
    error('radblockproc:numfileblockposscalar', 'Number of overlapping files (NUM_OVERLAP) is not a positive scalar.');
end
if mod(num_overlap, 1)
    num_overlap             = round(num_overlap);
    warning('radblockproc:roundoverlap', ['Number of overlapping files (NUM_OVERLAP) rounded to ' num2str(num_overlap) '.'])
end
if (num_overlap < 0)
    error('radblockproc:nooverlap', 'Number of overlapping files (NUM_OVERLAP) must be positive.')
end
if (num_overlap > num_file_block)
    error('radblockproc:overlapsize', ['Number of overlapping files (NUM_OVERLAP = ' num2str(num_overlap) ') is larger than the number of files per block (NUM_FILE_BLOCK = ' num2str(num_file_block) ').'])
end
if (~isnumeric(lat_std) || ~isscalar(lat_std))
    error('radblockproc:latstdscalar', 'Standard parallel (LAT_STD) is not a scalar.')
end
if (abs(lat_std) > 90)
    error('radblockproc:latstdrange', ['Absolute value of standard parallel (LAT_STD = ' num2str(abs(lat_std)) ') is larger than 90 degrees.'])
end
if ~ischar(dir_out)
    error('radblockproc:diroutstr', 'Block output directory (DIR_OUT) is not a string.')
end
if (~isempty(dir_out) && ~exist(dir_out, 'dir'))
    error('radblockproc:nodirout', 'Block output directory (DIR_OUT) does not exist.')
end
if (~license('test', 'map_toolbox') && ~exist('ll2ps', 'file'))
    error('radblockproc:ll2ps', 'Function LL2PS is not available within this user''s path and Mapping Toolbox is not available.')
end

% add GIMP projection if available
if (license('test', 'map_toolbox') && exist('mat/gimp_90m.mat' , 'file'))
    load mat/gimp_90m gimp_info
    gimp_avail              = true;
else
    gimp_avail              = false;
end
    
% begin file setup
file_in_all                 = dir([dir_in file_in '.mat']); % structure containing all file names matching pattern
num_file                    = length(file_in_all); % number of files
if ~num_file
    error('radblockproc:nofiles', ['No files matching "' file_in '" found in ' dir_in '.'])
end
tmp1                        = 1:num_file_block;
num_block                   = 1;
while true % sort out number of blocks based on number of files and desired overlap
    if (tmp1(end) >= num_file) % stop once number of files in tmp exceeds correct number
        break
    end
    tmp1                    = (tmp1(end) + 1 - num_overlap):(tmp1(end) - num_overlap + num_file_block);
    num_block               = num_block + 1;
end

for ii = 1:num_block
    
    tic
    
    block                   = struct;
    block.call              = struct;
    block.param             = struct;
    
    % preserve function call
    [block.call.dir_in, block.call.file_basename, block.call.block_basename, block.call.num_file_block, block.call.num_overlap, block.call.lat_std, block.call.dir_out] ...
                            = deal(dir_in, file_in, file_block, num_file_block, num_overlap, lat_std, dir_out); % save function call
    
    % get ready to load by determining the breakout between files and blocks
    jj                      = 0; % curr_files local counter
    if (ii == 1)
        num_start           = 1;
    else
        num_start           = 1 + ((ii - 1) * (num_file_block - num_overlap)); % file number to start at, accounting for overlap
    end
    curr_files              = num_start:(num_start + num_file_block - 1); % current file numbers to be blocked
    if any(curr_files > num_file) % last block may need shortening
        curr_files          = num_start:num_file;
    end
        
    [block.param.array_param, block.param.param_combine_wf_chan, block.param.param_csarp, block.param.param_radar, block.param.param_records] ...
                            = deal(cell(1, length(curr_files)));
    
    % read in each data file
    length_in               = zeros(length(curr_files), 1);
    for kk = curr_files
        jj                  = jj + 1;
        tmp1                = load([dir_in file_in_all(kk).name]); % data structure
        [tmp1.data, tmp1.elev_air, tmp1.time, tmp1.lat, tmp1.lon, tmp1.twtt_surf, tmp1.twtt] ...
                            = deal(tmp1.Data, tmp1.Elevation, tmp1.GPS_time, tmp1.Latitude, tmp1.Longitude, tmp1.Surface, tmp1.Time);
        if isfield(tmp1, 'Bottom') % necessary if Bottom field does not exist in, e.g., accumulation radar
            tmp1.twtt_bed   = tmp1.Bottom;
        else
            tmp1.twtt_bed   = NaN(size(tmp1.Surface));
        end
        tmp1                = rmfield(tmp1, {'Data' 'Elevation' 'GPS_time' 'Latitude' 'Longitude' 'Surface' 'Time'});
        if isfield(tmp1, 'Bottom')
            tmp1            = rmfield(tmp1, 'Bottom');
        end
        if isfield(tmp1, 'array_param')
            block.param.array_param{jj} ...
                            = tmp1.array_param;
        end
        if isfield(tmp1, 'param_combine_wf_chan')
            block.param.param_combine_wf_chan{jj} ...
                            = tmp1.param_combine_wf_chan;
        end
        if isfield(tmp1, 'param_csarp')
            block.param.param_csarp{jj} ...
                            = tmp1.param_csarp;
        end
        if isfield(tmp1, 'param_radar')
            block.param.param_radar{jj} ...
                            = tmp1.param_radar;
        end
        if isfield(tmp1, 'param_records')
            block.param.param_records{jj} ...
                            = tmp1.param_records;
        end
        if isfield(tmp1, 'param_get_heights')
            block.param.param_get_heights{jj} ...
                            = tmp1.param_get_heights;
        end
        if isfield(tmp1, 'param_qlook')
            block.param.param_qlook{jj} ...
                            = tmp1.param_qlook;
        end
        
        tmp1.data(~tmp1.data) ...
                            = NaN; % deal with zeros
        tmp1.data           = tmp1.data ./ max(tmp1.data(:)); % normalize data
        load_struct(jj)     = tmp1; %#ok<AGROW>
        clear tmp1
        length_in(jj)       = length(load_struct(jj).lat(~isnan(load_struct(jj).lat))); % length of structure
    end
       
    % begin distributing block pieces, starting with latitude and longitude
    block.lat               = [load_struct(:).lat];
    block.lon               = [load_struct(:).lon];
    
    if (sign(block.lat(find(~isnan(block.lat), 1))) ~= sign(lat_std))
        error('radblockproc:latmatch', 'Standard parallel is in the wrong hemisphere.')
    end
        
    try
        block.amp           = single([load_struct(:).data]); % amplitude
        block.twtt          = [load_struct(1).twtt]; % traveltime
    catch %#ok<CTCH> % different block sizes
        tmp2                = zeros(length(curr_files), 2);
        for jj = 1:length(curr_files)
            tmp2(jj, :)     = size(load_struct(jj).data);
        end
        [tmp4, tmp5]        = max(tmp2(:, 1));
        for jj = 1:length(curr_files)
            if (tmp2(jj, 1) < tmp4)
                load_struct(jj).data ...
                            = [load_struct(jj).data; NaN((tmp4 - tmp2(jj, 1)), tmp2(jj, 2))]; %#ok<AGROW>
            end
        end
        block.amp           = single([load_struct(:).data]);
        block.twtt          = [load_struct(tmp5).twtt];
    end
    [block.num_sample, block.num_trace] ...
                            = size(block.amp);
    
    % individual file names
    block.file_in           = {file_in_all(curr_files).name}';
    
    % address missing gps data by dropping those points altogether
    if any(isnan(block.lat))
        ind_nonan           = find(~isnan(block.lat));
        block.num_trace     = length(ind_nonan);
        [block.amp, block.lat, block.lon] ...
                            = deal(block.amp(:, ind_nonan), block.lat(ind_nonan), block.lon(ind_nonan));
        gps_nan             = true;
    else
        gps_nan             = false;
    end
    
    % record indices at which overlapping with previous/next blocks begins/ends
    tmp2                    = sum(length_in(1:(length(curr_files) - num_overlap)));
    block.ind_overlap       = NaN(1, 2);
    if num_overlap
        if (ii == 1)
            block.ind_overlap(2)    = tmp2 + 1;
        elseif ((ii > 1) && (ii < num_block))
            block.ind_overlap       = [(sum(length_in(1:num_overlap))) (tmp2 + 1)];
        elseif (ii == num_block)
            block.ind_overlap(1)    = sum(length_in(1:num_overlap));
        end
    end
    
    if (license('checkout', 'map_toolbox') && (block.lat(1) > 0)) % in Greenland and can do better x/y coordinates
        if (ii == 1)
            wgs84           = almanac('earth', 'wgs84', 'meters');
            ps_struct       = defaultm('ups');
            [ps_struct.geoid, ps_struct.mapparallels, ps_struct.falsenorthing, ps_struct.falseeasting, ps_struct.origin] ...
                            = deal(wgs84, lat_std, 0, 0, [90 -45 0]);
        end
        [block.x, block.y]  = mfwdtran(ps_struct, block.lat, block.lon);
        if gimp_avail
            [block.x_gimp, block.y_gimp] ...
                            = projfwd(gimp_info, block.lat, block.lon); % GIMP projection
            [block.x_gimp, block.y_gimp] ...
                            = deal((1e-3 .* block.x_gimp), (1e-3 .* block.y_gimp)); % m to km
        end
    else
        % convert lat/lon to polar stereographic x/y (better for Antarctica then Greenland)
        [block.x, block.y]  = ll2ps(block.lat, block.lon, lat_std);
    end
    [block.x, block.y]      = deal((1e-3 .* block.x), (1e-3 .* block.y)); % m to km
    
    block.dist              = cumsum([0 sqrt((diff(block.x) .^ 2) + (diff(block.y) .^ 2))]); % distance vector
    block.dist_lin          = interp1([1 block.num_trace], block.dist([1 end]), 1:block.num_trace);
    if (ii > 1) % increment distance vectors along transect, i.e., don't restart them for each block
        block.dist          = block.dist + tmp3(3) + sqrt(((block.x(1) - tmp3(1)) ^ 2) + ((block.y(1) - tmp3(2)) ^ 2));
        block.dist_lin      = block.dist_lin + tmp3(4) + sqrt(((block.x(1) - tmp3(1)) ^ 2) + ((block.y(1) - tmp3(2)) ^ 2)); % monotonically increasing distance vector that is easier for imagesc plots
    end
    
    % do the same for the GIMP projection if available
    if (license('checkout', 'map_toolbox') && (block.lat(1) > 0) && gimp_avail)
        block.dist_gimp     = cumsum([0 sqrt((diff(block.x_gimp) .^ 2) + (diff(block.y_gimp) .^ 2))]);
        block.dist_lin_gimp = interp1([1 block.num_trace], block.dist_gimp([1 end]), 1:block.num_trace);
        if (ii > 1)
            block.dist_gimp = block.dist_gimp + tmp3_gimp(3) + sqrt(((block.x_gimp(1) - tmp3_gimp(1)) ^ 2) + ((block.y_gimp(1) - tmp3_gimp(2)) ^ 2));
            block.dist_lin_gimp ...
                            = block.dist_lin_gimp + tmp3_gimp(4) + sqrt(((block.x_gimp(1) - tmp3_gimp(1)) ^ 2) + ((block.y_gimp(1) - tmp3_gimp(2)) ^ 2));
        end
    end
    
    % aircraft elevation
    block.elev_air          = [load_struct(:).elev_air];    
    if gps_nan
        block.elev_air      = block.elev_air(ind_nonan);
    end
    
    % make sure traveltime is a column vector
    if isrow(block.twtt)
        block.twtt          = block.twtt'; 
    end
    block.dt                = block.twtt(2) - block.twtt(1);
    
    % measurement time
    block.time              = [load_struct(:).time];
    if gps_nan
        block.time          = block.time(ind_nonan);
    end
    
    % surface arrival auto-pick
    block.twtt_surf         = [load_struct(:).twtt_surf];
    if gps_nan
        block.twtt_surf     = block.twtt_surf(ind_nonan);
    end
    
    % bed pick
    block.twtt_bed          = [load_struct(:).twtt_bed];
    if gps_nan
        block.twtt_bed      = block.twtt_bed(ind_nonan);
    end
    
    % order fields alphabetically
    block                   = orderfields(block);
    block.call              = orderfields(block.call);
    
    % save block
    if (ii < 10)
        save([dir_out file_block '_0' num2str(ii)], '-v7.3', 'block')
        disp(['Processed ' num2str(length(curr_files)) ' files for block ' num2str(ii) ' / ' num2str(num_block) ' and saved as ' file_block '_0' num2str(ii) ' in ' dir_out ' in ' num2str(toc, '%4.2f') ' s.'])
    else
        save([dir_out file_block '_' num2str(ii)], '-v7.3', 'block')
        disp(['Processed ' num2str(length(curr_files)) ' files for block ' num2str(ii) ' / ' num2str(num_block) ' and saved as ' file_block '_' num2str(ii) ' in ' dir_out ' in ' num2str(toc, '%4.2f') ' s.'])
    end
    
    % get x/y/distance values at end of data before overlap, so that distance vector can stay sane
    tmp3                    = [block.x(tmp2) block.y(tmp2) block.dist(tmp2) block.dist_lin(tmp2)];
    tmp3_gimp               = [block.x_gimp(tmp2) block.y_gimp(tmp2) block.dist_gimp(tmp2) block.dist_lin_gimp(tmp2)];
    
    % clear the big stuff
    clear load_struct block
end