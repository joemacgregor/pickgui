function data_cat			= radframeproc(dir_in, file_in, dir_cat, file_cat, do_norm, do_tol, x_dem, y_dem, z_dem, name_dem, dir_layer)
% RADFRAMEPROC Concatenate and pre-process CReSIS ice-penetrating radar data frames into larger continuous sets for PICKGUI v2.
%
%   RADFRAMEPROC(DIR_IN,FILE_IN,DIR_CAT,FILE_CAT,DO_NORM,XY_TOL,X_DEM,Y_DEM,Z_DEM,NAME_DEM)
%	pre-processes CReSIS radar data frames stored in DIR_IN whose
%   names are contained in the cell string FILE_IN and then saves them in
%   DIR_CAT with the name FILE_CAT. DO_NORM is a logical scalar that, if
%   true, normalizes the radar data by their maximum value. DO_TOL is a
%   logical scalar that, if true, removes non-unique positions using
%   UNIQUETOL; if false, then only exactly repeated locations are removed.
%   The local surface elevation Z_DEM in meters is extracted from digital
%   elevation model NAME_DEM, from DEM-corrected elevation variables are
%   generated for the radar data. RADFRAMEPROC assumes that vectors X_DEM
%   and Y_DEM are coordinates are in meters but otherwise the same
%   projection as that selected based on the hemisphere (EPSG:3413 for
%   positive latitudes, i.e., Greenland, and EPSG:3031 for negative
%   latitudes, i.e., Antarctica).
%
%	RADFRAMEPROC(...,DIR_LAYER) loads a filename matching FILE_CAT from
%	DIR_LAYER to overwrite existing surface and bed picks.
% 
%	DATA_CAT = RADFRAMPROC(...) outputs the concatenated data file
%	*instead* of saving it.
% 
% Joe MacGregor (NASA)
% Last updated: 12 August 2024

if ~license('test', 'map_toolbox')
	error('radframeproc:mapping', 'Mapping Toolbox license is required.')
end
if ~any(nargin == [10 11])
	error('radframeproc:inputs', 'Number of inputs should be 10 or 11.')
end
if ~ischar(dir_in)
	error('radframeproc:dirinstr', 'Data input directory (DIR_IN) is not a string.')
end
if (~isempty(dir_in) && ~exist(dir_in, 'dir'))
	error('radframeproc:nodirin', 'Data input directory (DIR_IN) does not exist.')
end
if ~iscell(file_in)
	error('radframeproc:fileincell', 'List of file names (FILE_IN) is not a cell.')
end
num_file					= length(file_in);
for ii = 1:num_file
	if ~ischar(file_in{ii})
		error('radframeproc:fileinstr', ['FILE_IN(' num2str(ii) ') is not a string.'])
	end
	if ~exist([dir_in file_in{ii}], 'file')
		error('radframeproc:fileinexist', ['FILE_IN(' num2str(ii) ') does not exist.'])
	end
end
if ~ischar(dir_cat)
	error('radframeproc:dircatstr', 'Processed set output directory (DIR_CAT) is not a string.')
end
if (~isempty(dir_cat) && ~exist(dir_cat, 'dir'))
	error('radframeproc:nodircat', 'Processed set output directory (DIR_CAT) does not exist.')
end
if ~ischar(file_cat)
	error('radframeproc:namecatstr', 'Ouput concatened filename (FILE_CAT) is not a string.')
end
if (~islogical(do_norm) || ~isscalar(do_norm))
	error('radframeproc:donorm', 'DO_NORM must be a logical scalar.')
end
if (~islogical(do_tol) || ~isscalar(do_tol))
	error('radframeproc:dotol', 'DO_TOL must be a logicalscalar.')
end
if (~isnumeric(x_dem) || ~isrow(x_dem))
	error('radframeproc:xdem', 'X_DEM is not a numeric row vector.')
end
if (~isnumeric(y_dem) || isrow(y_dem))
	error('radframeproc:ydem', 'Y_DEM is not a numeric column vector.')
end
if (~isnumeric(z_dem) || ~ismatrix(z_dem))
	error('radframeproc:dem1', 'Z_DEM is not a numeric matrix.')
end
if ((size(z_dem, 2) ~= length(x_dem)) || (size(z_dem, 1) ~= length(y_dem)))
	error('radframeproc:dem2', 'Z_DEM size is not compatible with X_DEM and Y_DEM.')
end
if ~ischar(name_dem)
	error('radframeproc:namedem', 'NAME_DEM is not a string.')
end
if (nargin == 11)
	if ~ischar(dir_layer)
		error('radframeproc:dirlayer', 'DIR_LAYER is not a string.')
	end
	if (~isempty(dir_layer) && ~exist(dir_layer, 'dir'))
		error('radframeproc:nodirlayer', 'DIR_LAYER does not exist.')
	end
	do_layer				= true;
else
	do_layer				= false;
end
if (nargout > 1)
	error('radframeproc:nargout', 'RADFRAMEPROC can only have one output (DATA_CAT).')
end

% DEM correction prep
dist_pad					= 2e3; % pad distance for spline, m
speed_vacuum				= 299792458; % speed of light in the vacuum, m/s
permitt_ice					= 3.15; % permittivity of ice, dimensionless
speed_ice					= speed_vacuum / sqrt(permitt_ice); % speed in ice, m/s

% preserve function call and radar parameters
radframeproc_call			= struct('dir_in', dir_in, 'file_in', {file_in}, 'file_cat', file_cat, 'dir_cat', dir_cat, 'do_norm', do_norm, 'do_tol', do_tol, 'name_dem', name_dem);
if (nargin == 11)
	radframeproc_call(1).dir_layer ...
							= dir_layer;
end
radar_param					= struct;

% read in each data file
for ii = 1:num_file
	
	radar_data_tmp			= load([dir_in file_in{ii}]); % load frame into temporary structure
	
	if do_layer
		try
			if (length(file_in{ii}) == 24) % normal file
				layer_data_tmp ...
							= load([dir_layer file_in{ii}], 'layerData'); % load layerData associated with frame
				if isempty(fieldnames(layer_data_tmp))
					layer_data_tmp ...
							= load([dir_layer file_in{ii}], 'twtt', 'gps_time'); % load layerData associated with frame
					layer_data_tmp.layerData ...
							= cell(1, 2);
					layer_data_tmp.layerData{1}.value ...
							= cell(1, 2);
					layer_data_tmp.layerData{1}.value{2}.data ...
							= layer_data_tmp.twtt(1, :);
					layer_data_tmp.layerData{2}.value ...
							= cell(1, 2);
					layer_data_tmp.layerData{2}.value{2}.data ...
							= layer_data_tmp.twtt(2, :);
					disp('stupid 2018/9 layerData')
				end
			else
				layer_data_tmp ...
							= load([dir_layer file_in{ii}([1:5 13:end])], 'layerData'); % load layerData associated with img frame
			end
		catch
			disp(['!!! layerData file matching ' file_in{ii} ' does not exist !!!'])
			do_layer		= false;
		end
	end
	
	% new field names
	[radar_data_tmp.amp, radar_data_tmp.elev_air, radar_data_tmp.time, radar_data_tmp.lat, radar_data_tmp.lon, radar_data_tmp.twtt_surf, radar_data_tmp.twtt] ...
							= deal(radar_data_tmp.Data, radar_data_tmp.Elevation, radar_data_tmp.GPS_time, radar_data_tmp.Latitude, radar_data_tmp.Longitude, radar_data_tmp.Surface, radar_data_tmp.Time);
	
	if isfield(radar_data_tmp, 'Bottom') % Bottom field does not exist for some radar data types
		radar_data_tmp.twtt_bed	...
							= radar_data_tmp.Bottom;
	else
		radar_data_tmp.twtt_bed	...
							= NaN(1, length(radar_data_tmp.Surface));
	end
	
	% use layerData picks if available AND they make sense
	if do_layer
		if (length(layer_data_tmp.layerData{1}.value{2}.data) == length(radar_data_tmp.lat))
			radar_data_tmp.twtt_surf ...
							= layer_data_tmp.layerData{1}.value{2}.data;
		elseif isfield(layer_data_tmp, 'gps_time') % 2018 fix
			ind_surf_tmp	= interp1(radar_data_tmp.twtt, 1:length(radar_data_tmp.twtt), interp1(layer_data_tmp.gps_time, layer_data_tmp.layerData{1}.value{2}.data, radar_data_tmp.time, 'linear', 'extrap'), 'nearest');
			radar_data_tmp.twtt_surf(~isnan(ind_surf_tmp)) ...
							= radar_data_tmp.twtt(ind_surf_tmp(~isnan(ind_surf_tmp)))';
		else
			disp(['!!! layerData surface traveltime length does not match number of traces in ' file_in{ii} ' !!!'])
		end
		if (length(layer_data_tmp.layerData{2}.value{2}.data) == length(radar_data_tmp.lat))
			radar_data_tmp.twtt_bed ...
							= layer_data_tmp.layerData{2}.value{2}.data;
		elseif isfield(layer_data_tmp, 'gps_time') % 2018 fix
			ind_bed_tmp		= interp1(radar_data_tmp.twtt, 1:length(radar_data_tmp.twtt), interp1(layer_data_tmp.gps_time, layer_data_tmp.layerData{2}.value{2}.data, radar_data_tmp.time, 'linear', 'extrap'), 'nearest');
			radar_data_tmp.twtt_bed(~isnan(ind_bed_tmp)) ...
							= radar_data_tmp.twtt(ind_bed_tmp(~isnan(ind_bed_tmp)))';
		else
			disp(['!!! layerData bed traveltime length does not match number of traces in ' file_in{ii} ' !!!'])
		end
	end
	
	% remove old fields
	radar_data_tmp			= rmfield(radar_data_tmp, {'Data' 'Elevation' 'GPS_time' 'Latitude' 'Longitude' 'Surface' 'Time'});	
	if isfield(radar_data_tmp, 'Bottom')
		radar_data_tmp		= rmfield(radar_data_tmp, 'Bottom');
	end
	
	% check for other fields, record as possible and then delete
	if isfield(radar_data_tmp, 'Heading')
		radar_param.azimuth{ii} ...
							= radar_data_tmp.Heading;
		radar_data_tmp		= rmfield(radar_data_tmp, 'Heading');
	end
	if isfield(radar_data_tmp, 'Pitch')
		radar_param.pitch{ii} ...
							= radar_data_tmp.Pitch;
		radar_data_tmp		= rmfield(radar_data_tmp, 'Pitch');
	end
	if isfield(radar_data_tmp, 'Roll')
		radar_param.roll{ii}= radar_data_tmp.Roll;
		radar_data_tmp		= rmfield(radar_data_tmp, 'Roll');
	end
	if isfield(radar_data_tmp, 'array_param')
		radar_param.array_param{ii} ...
							= radar_data_tmp.array_param;
		radar_data_tmp		= rmfield(radar_data_tmp, 'array_param');
	end
	if isfield(radar_data_tmp, 'param_combine_wf_chan')
		radar_param.param_combine_wf_chan{ii} ...
							= radar_data_tmp.param_combine_wf_chan;
		radar_data_tmp		= rmfield(radar_data_tmp, 'param_combine_wf_chan');
	end
	if isfield(radar_data_tmp, 'param_csarp')
		radar_param.param_csarp{ii} ...
							= radar_data_tmp.param_csarp;
		radar_data_tmp		= rmfield(radar_data_tmp, 'param_csarp');
	end
	if isfield(radar_data_tmp, 'param_radar')
		radar_param.param_radar{ii} ...
							= radar_data_tmp.param_radar;
		radar_data_tmp		= rmfield(radar_data_tmp, 'param_radar');
	end
	if isfield(radar_data_tmp, 'param_records')
		radar_param.param_records{ii} ...
							= radar_data_tmp.param_records;
		radar_data_tmp		= rmfield(radar_data_tmp, 'param_records');
	end
	if isfield(radar_data_tmp, 'param_get_heights')
		radar_param.param_get_heights{ii} ...
							= radar_data_tmp.param_get_heights;
		radar_data_tmp		= rmfield(radar_data_tmp, 'param_get_heights');
	end
	if isfield(radar_data_tmp, 'param_qlook')
		radar_param.param_qlook{ii} ...
							= radar_data_tmp.param_qlook;
		radar_data_tmp		= rmfield(radar_data_tmp, 'param_qlook');
	end
	if isfield(radar_data_tmp, 'param_vectors')
		radar_param.param_vectors{ii} ...
							= radar_data_tmp.param_vectors;
		radar_data_tmp		= rmfield(radar_data_tmp, 'param_vectors');
	end
	
	radar_data_tmp.amp((radar_data_tmp.amp == 0) | isinf(radar_data_tmp.amp)) ...
							= NaN; % deal with zeros
	
	if do_norm
		radar_data_tmp.amp	= radar_data_tmp.amp ./ max(radar_data_tmp.amp, [], 'all', 'omitnan'); % normalize data
	end
	
	if (ii == 1)
		radar_data			= radar_data_tmp;
	else
		radar_data(ii)		= radar_data_tmp;
	end
end

clear radar_data_tmp

% begin distributing block pieces, starting with latitude and longitude, deg
lat							= [radar_data(:).lat];
lon							= [radar_data(:).lon];

% concatenate amplitude matrices
try
	amp						= single([radar_data(:).amp]); % amplitude
	twtt					= [radar_data(1).twtt]; % traveltime
catch %#ok<CTCH> % different data frame sizes uh oh
	warning('RADARGRAMS TO BE CONCATENATED HAVE DIFFERENT NUMBERS OF SAMPLES...COULD BE BAD...CERTAINLY ANNOYING...')
	dt_tmp					= zeros(1, num_file);
	for ii = 1:num_file
		dt_tmp(ii)			= radar_data(ii).twtt(2) - radar_data(ii).twtt(1);
	end
	if ~isempty(find((diff(dt_tmp) ~= 0), 1))
		warning('RADARGRAMS HAVE DIFFERENT FAST TIME SAMPLING INTERVALS...YIKES...CANNOT BE CONCATENATED SO ABORTING. REPORTING FILE NUMBERS AND INTERVALS (ns).')
		disp([(1:num_file)' (1e9 .* dt_tmp)'])
		return
	end
	sz_radar_data			= zeros(num_file, 2); % amplitude matrix sizes
	for ii = 1:num_file
		sz_radar_data(ii, :) ...
							= size(radar_data(ii).amp);
	end
	[tmp1, tmp2]			= max(sz_radar_data(:, 1));
	for ii = find(sz_radar_data(:, 1) < tmp1)' % fill in missing traveltimes in smaller matrices with NaNs, instead of just throwing NaNs on the end like I used to
		radar_data(ii).amp	= radar_data(ii).amp(interp1(radar_data(ii).twtt, 1:sz_radar_data(ii, 1), radar_data(tmp2).twtt, 'nearest', 'extrap'), :); % find nearby traveltime and link up...not ideal but not horrible
	end
	amp						= single([radar_data(:).amp]);
	twtt					= [radar_data(tmp2).twtt];
end
[num_sample, num_trace]		= size(amp);

% convert to dB
amp							= 10 .* log10(abs(amp));

% address missing GNSS data by dropping those points altogether
if (~isempty(find(isnan(lat), 1)) || ~isempty(find(isnan(lon), 1)))
	ind_good_gnss			= find(~isnan(lat) & ~isnan(lon));
	warning('NaN latitudes or longitudes detected!')
else
	ind_good_gnss			= 1:num_trace;
end

% measurement time, slow time GPS time (not UTC), s
time						= [radar_data(:).time];

% search for unique times/positions
do_unique					= false;
[~, ind_unique_time]		= unique(time, 'stable');
if (length(ind_unique_time) < num_trace) % parts of 19990513_01 particularly challenging on time
	do_unique				= true;
	warning('Non-unique measurement times detected!')
end
if do_tol % 2010-2013 looking at you
	[~, ind_unique_pos]		= uniquetol([lat' lon'], (0.5 * median(diff(cumsum([0 distance([lat(1:(end - 1))' lon(1:(end - 1))'], [lat(2:end)' lon(2:end)'], wgs84Ellipsoid)'])))), 'ByRows', true, 'DataScale', 1);
	ind_unique_pos			= sort(ind_unique_pos);
else
	[~, ind_unique_pos]		= unique([lat' lon'], 'rows', 'stable');
end
if (length(ind_unique_pos) < num_trace)
	if do_unique
		ind_unique			= intersect(ind_unique_time, ind_unique_pos);
	else
		ind_unique			= ind_unique_pos;
	end
	do_unique				= true;
	disp(length(ind_unique))
	disp(num_trace)	
	warning('Non-unique positions detected!')
end
if ~do_unique
	ind_unique				= 1:num_trace;
end
ind_good					= intersect(ind_good_gnss, ind_unique);

if (length(ind_good) < num_trace)
	do_fix					= true;
	num_trace				= length(ind_good);	
else
	do_fix					= false;
end

if do_fix
	[amp, lat, lon, time]	= deal(amp(:, ind_good), lat(ind_good), lon(ind_good), time(ind_good));	
end

% distance
dist						= cumsum([0 distance([lat(1:(end - 1))' lon(1:(end - 1))'], [lat(2:end)' lon(2:end)'], wgs84Ellipsoid)']); % great-circle distance, m
dist_lin					= interp1([1 num_trace], dist([1 end]), 1:num_trace); % simplified monotonic distance, m

% project lat/lon to polar stereographic
if (lat(1) > 0)
	[x, y]					= projfwd(projcrs(3413), lat, lon); % EPSG:3413 for Greenland
else	
	[x, y]					= projfwd(projcrs(3031), lat, lon); % EPSG:3031 for Antarctica
end

% aircraft elevation
elev_air					= [radar_data(:).elev_air];
if do_fix
	elev_air				= elev_air(ind_good);
end

% make sure traveltime is a column vector
if isrow(twtt)
	twtt					= twtt';
end

% surface arrival auto-pick
twtt_surf					= [radar_data(:).twtt_surf];
if do_fix
	twtt_surf				= twtt_surf(ind_good); % s
end

% bed pick
twtt_bed					= [radar_data(:).twtt_bed];
if do_fix
	twtt_bed				= twtt_bed(ind_good); % s
end

elev_surf					= elev_air - (twtt_surf .* (speed_vacuum / 2));
thick						= (twtt_bed - twtt_surf) .* (speed_ice / 2);
elev_bed					= elev_surf - thick;

% determine DEM elevations
ind_x						= find((x_dem >= (min(x) - dist_pad)) & (x_dem <= (max(x) + dist_pad)));
ind_y						= find((y_dem >= (min(y) - dist_pad)) & (y_dem <= (max(y) + dist_pad)));

if ((length(ind_x) > 1) && (length(ind_y) > 1))
	elev_surf_corr			= interp2(x_dem(ind_x), y_dem(ind_y), z_dem(ind_y, ind_x), x, y, 'spline'); % DEM elevation
	elev_bed_corr			= elev_bed + (elev_surf_corr - elev_surf); % DEM-corrected bed elevation
elseif ~isempty(find(~isnan(z_dem), 1))
	error('radframeproc:demarea', 'DEM area does not include radar data coverage!')
else
	[elev_bed_corr, elev_surf_corr]	...
							= deal(NaN(1, num_trace));
end

% trim dataset above and below surface and bed picks if those are continuous
if isempty(find(isnan(twtt_surf), 1))
	ind_trim_surf			= min(interp1(twtt', 1:num_sample, twtt_surf, 'nearest', 'extrap'));
	[amp, twtt]				= deal(amp(ind_trim_surf:end, :), twtt(ind_trim_surf:end));
	num_sample				= size(amp, 1);
else
	ind_trim_surf			= NaN;
end
if isempty(find(isnan(twtt_bed), 1))
	ind_trim_bed			= max(interp1(twtt', 1:num_sample, twtt_bed, 'nearest', 'extrap'));
	[amp, twtt]				= deal(amp(1:ind_trim_bed, :), twtt(1:ind_trim_bed));
	num_sample				= size(amp, 1);
else
	ind_trim_bed			= NaN;
end
ind_trim_bed				= ind_trim_bed + ind_trim_surf; % account for surface offset also

% order fields alphabetically
radframeproc_call           = orderfields(radframeproc_call);
radar_param					= orderfields(radar_param);

% save concatenated radar data
if nargout
	data_cat				= struct('amp', amp, 'dist', dist, 'dist_lin', dist_lin, 'elev_air', elev_air, 'elev_bed', elev_bed, 'elev_bed_corr', elev_bed_corr, 'elev_surf', elev_surf, 'elev_surf_corr', elev_surf_corr, ...
									 'ind_trim_bed', ind_trim_bed, 'ind_trim_surf', ind_trim_surf, 'lat', lat, 'lon', lon, 'num_sample', num_sample, 'num_trace', num_trace, 'radframeproc_call', radframeproc_call, ...
									 'radar_param', radar_param, 'time', time, 'thick', thick, 'twtt', twtt, 'twtt_bed', twtt_bed, 'twtt_surf', twtt_surf, 'x', x, 'y', y);
else
	save([dir_cat file_cat], '-v7.3', 'amp', 'dist', 'dist_lin', 'elev_air', 'elev_bed', 'elev_bed_corr', 'elev_surf', 'elev_surf_corr', 'ind_trim_bed', 'ind_trim_surf', 'lat', 'lon', 'num_sample', 'num_trace', ...
									  'radframeproc_call', 'radar_param', 'time', 'thick', 'twtt', 'twtt_bed', 'twtt_surf', 'x', 'y')
	disp(['Processed ' num2str(num_file) ' files for concatenated file saved as ' file_cat ' in ' dir_cat '.'])
end