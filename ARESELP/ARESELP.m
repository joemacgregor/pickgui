function ARESELP(frame_or_cat, dir_in, file_in, parallel_check, wavelet_scale, wavelet_shp, pk_prom_min, twtt_sep_min, num_pk_max, len_blk, twtt_jump_max, angle_diff, twtt_proxim_max, len_layer_min)
% ARESELP Modified version of Automated RES Englacial Layer-tracing Package (ARESELP) developed by Xiong et al. (2018).
%   
%   ARESELP uses wavelet and Hough transforms to automatically identify
%   candidate layers in a radargram. In this implementation, concatenated
%   radargrams outputted by RADFRAMEPROC are loaded, analyzed and then
%   saved with an additional M-object structure, layer_ARESELP, which are
%   the indices associated with the M candidate layers. These candidate
%   layers can be used by PICKGUI to predict the radiostratigraphy and then
%   flatten the radargram for initial processing.
%   
%   ARESELP was originally developed and described by:
%   Xiong, S., Muller, J.P. and Carretero, R.C. (2018), A New Method for
%   Automatically Tracing Englacial Layers from MCoRDS Data in NW
%   Greenland, Remote Sensing, 10(43), doi:10.3390/rs10010043.
%	
%   ARESELP is based on code archived at:
%   https://github.com/xiongsiting/ARESELP/
%   
% Joe MacGregor (NASA)
% Last updated: 12 August 2024

% EXAMPLE SETTING BELOW
% clear
% nargin						= 12;
% frame_or_cat				= 'cat'; % use either individual CReSIS data frame ('frame') or radframeproc-concatenated set of frames ('cat')
% dir_in						= '/Users/jamacgre/Desktop/work/'; % directory where data are located
% file_in						= '*'; % files to process
% parallel_check				= true; % check for Parallel Computing Toolbox intialization
% ind_thick_min				= 50; % minimum number of vertical samples for wavelet transform peak finding
% wavelet_scale				= 3:15; % wavelet scales
% wavelet_shp					= 'mexh'; % wavelet shape
% pk_prom_min					= 5; % minimum peak prominence, dimensionless
% seed_pt_obj_cutoff			= 10; % cutoff number of pixels for removing small objects
% twtt_sep_min				= 0.5e-6; % minimum vertical sample separation of layers, s
% num_pk_max					= 25; % maximum number of wavelet peaks per trace, #
% len_blk						= 1e3; % length of horizontal block for Hough slope calculation, m
% twtt_jump_max				= 0.3e-6; % maximum vertical jump between blocks, s
% angle_diff					= 20; % slope angle difference, deg
% twtt_proxim_max				= 0.05e-6; % maximum proximity of a new layer before merging with an old one, s
% len_layer_min				= 10e3; % minimum length of layer, m
% % num_test					= 3e3; % number of horizontal samples to use when testing/debugging

if ~license('test', 'wavelet_toolbox')
    error('areselp:wavelet', 'Wavelet Toolbox license is not available but needed for continuous wavelet transforms.')
end
if ~license('test', 'statistics_toolbox')
    error('areselp:stats', 'Statistics and Machine Learning Toolbox license is not available but needed for wavelet distribution fitting.')
end
if ~license('test', 'image_toolbox')
    error('areselp:image', 'Image Processing Toolbox license is not available but needed for Hough transforms.')
end
if (nargin ~= 12)
    error('areselp:nargin', 'Not enough input arguments (need 12).')
end
if ~ischar(dir_in)
    error('areselp:dirinstr', 'File directory (DIR_IN) is not a string.')
end
if (~isempty(dir_in) && ~exist(dir_in, 'dir'))
    error('areselp:nodirin', 'File directory (DIR_IN) does not exist.')
end
if ~ischar(file_in)
    error('areselp:fileinstr', 'Filename (FILE_IN) is not a string.')
end
if (~islogical(parallel_check) || ~isscalar(parallel_check))
    error('aresp_obj:parallelcheck', 'True/false check for Parallel Computing Toolbox license (PARALLEL_CHECK) is not a logical scalar.')
end
if (~isnumeric(wavelet_scale) || ~isvector(wavelet_scale))
    error('areselp:waveletscale', 'Wavelet scale (WAVELET_SCALE) is not a numeric vector.')
end
if (~ischar(wavelet_shp) || ~any(strcmp(wavelet_shp, {'mexh' 'morl'})))
    error('areselp:waveletshp', 'Wavelet shape (WAVELET_SHP) is not a character string that is either ''mexh'' or ''morl''.')
end
if (~isnumeric(pk_prom_min) || ~isscalar(pk_prom_min))
    error('areselp:pkprommin', 'Minimum peak prominence of a CWT peak (PK_PROM_MIN) is not a numeric scalar.')
end
if (~isnumeric(twtt_sep_min) || ~isscalar(twtt_sep_min))
    error('areselp:twttsepmin', 'Minimum traveltime separation between layers (TWTT_SEP_MIN) is not a numeric scalar.')
end
if (~isnumeric(num_pk_max) || ~isscalar(num_pk_max))
    error('areselp:numpkmax', 'Maximum number of peaks in CWT to preserve per trace (NUM_PK_MAX) is not a numeric scalar.')
end
if (~isnumeric(len_blk) || ~isscalar(len_blk))
    error('areselp:lenblk', 'Length of slope blocks (LEN_BLK) is not a numeric scalar.')
end
if (~isnumeric(twtt_jump_max) || ~isscalar(twtt_jump_max))
    error('areselp:twttjumpmax', 'Maximum traveltime jump between slopes (TWTT_JUMP_MAX) is not a numeric scalar.')
end
if (~isnumeric(angle_diff) || ~isscalar(angle_diff))
    error('areselp:anglediff', 'Maximum slope angle difference (ANGLE_DIFF) is not a numeric scalar.')
end
if (~isnumeric(twtt_proxim_max) || ~isscalar(twtt_proxim_max))
    error('areselp:twttproximmax', 'Maximum traveltime proximity of a new layer before merging with an old one (TWTT_PROXIM_MAX) is not a numeric scalar.')
end
if (~isnumeric(len_layer_min) || ~isscalar(len_layer_min))
    error('areselp:lenlayermin', 'Minimum layer length (LEN_LAYER_MIN) is not a numeric scalar.')
end

% filenames in dir_in
name_file                   = dir([dir_in file_in '.mat']);
name_file                   = {name_file.name}';
num_file                    = length(name_file);
sum_sq_layer_len			= NaN(1, num_file);

disp(['Starting ARESELP for ' dir_in '...'])

% loop through each file in dir_in
for ii = 1:num_file
	
    disp([name_file{ii}(1:(end - 4)) ' (' num2str(ii) ' / ' num2str(num_file) ')...'])
	
	% load only necessary variables from current file
	switch frame_or_cat
		case 'frame'
			data_cat		= radframeproc(dir_in, name_file(ii), dir_in, [name_file{ii}(1:(end - 4)) '_for_areselp'], false, [-1e9 1e9], [-1e9 1e9]', NaN(2), 'null');
			fieldnames_curr = fieldnames(data_cat);
			data_cat        = rmfield(data_cat, fieldnames_curr(~contains(fieldnames_curr, {'amp' 'dist_lin' 'ind_trim_surf' 'num_sample' 'num_trace' 'twtt' 'twtt_bed' 'twtt_surf'})));
		case 'cat'
			data_cat        = load([dir_in name_file{ii}], 'amp', 'dist_lin', 'ind_trim_surf', 'num_sample', 'num_trace', 'twtt', 'twtt_bed', 'twtt_surf');
	end
	
	% % trim horizontally for testing only
	% [data_cat.amp, data_cat.dist_lin, data_cat.num_trace, data_cat.twtt_bed, data_cat.twtt_surf] ...
	% 						= deal(data_cat.amp(:, 1:num_test), data_cat.dist_lin(1:num_test), num_test, data_cat.twtt_bed(1:num_test), data_cat.twtt_surf(1:num_test));
	
	% indices of surface and bed picks
	ind_surf				= interp1(data_cat.twtt', 1:data_cat.num_sample, data_cat.twtt_surf, 'nearest', 'extrap');
	ind_bed					= interp1(data_cat.twtt', 1:data_cat.num_sample, data_cat.twtt_bed, 'nearest', 'extrap');
	ind_thick				= ind_bed - ind_surf;
	ind_twtt_sep_min		= ceil(twtt_sep_min / diff(data_cat.twtt(1:2)));
	ind_thick_enough		= find(ind_thick > (ind_thick_min + (2 * ind_twtt_sep_min))); % traces where ice is thick enough for further evaluation
	
	pk_cwt					= zeros(data_cat.num_sample, data_cat.num_trace);
	
	% wavelet transform on each trace (column)
	if parallel_check
		[tmp1, tmp2, tmp3, tmp4] ...
							= deal(data_cat.amp(:, ind_thick_enough), data_cat.amp(:, ind_thick_enough), ind_surf(ind_thick_enough), ind_bed(ind_thick_enough));
		parfor jj = 1:length(ind_thick_enough)
			tmp1(:, jj)		= ARESELP_find_pk_cwt(tmp2(:, jj), wavelet_scale, wavelet_shp, tmp3(jj), tmp4(jj), pk_prom_min, ind_twtt_sep_min, num_pk_max);
		end
		pk_cwt(:, ind_thick_enough) ...
							= tmp1;
	else
		for jj = ind_thick_enough
			pk_cwt(:, jj)	= ARESELP_find_pk_cwt(data_cat.amp(:, jj), wavelet_scale, wavelet_shp, ind_surf(jj), ind_bed(jj), pk_prom_min, ind_twtt_sep_min, num_pk_max);
		end
	end
	
	% only keep peaks greater than or equal to mean of lognormal distribution fit
	ind_pk_cwt				= find(pk_cwt > 0); % only care about positive peaks
	seed_pt					= [pk_cwt(ind_pk_cwt) ind_pk_cwt];
	pk_cwt_dist				= fitdist(seed_pt(:, 1), 'lognormal'); % fit peaks to a lognormal distribution for further culling
	seed_pt					= seed_pt((seed_pt(:, 1) >= exp(pk_cwt_dist.mu + ((pk_cwt_dist.sigma .^ 2) / 2))), 2);
	
	% remove isolated seed pts then bridge small gaps then remove small blobs
	seed_pt_mat				= false(data_cat.num_sample, data_cat.num_trace);
	seed_pt_mat(seed_pt)	= true; % logical matrix of seed pts	
	seed_pt_mat				= bwmorph(bwmorph(seed_pt_mat, 'clean'), 'bridge'); % remove isolated pts and then connect nearby pts
	seed_pt_obj				= bwconncomp(seed_pt_mat); % connected components of binary image
	num_seed_obj			= cellfun(@numel, seed_pt_obj.PixelIdxList); % number of pixels in each seed pt component
	num_seed_obj(num_seed_obj <= seed_pt_obj_cutoff) ...
							= 0;
	
	[num_seed_obj, ind_sort_obj] ...
							= sort(num_seed_obj, 'descend'); % reorder to descending object size
	seed_pt_obj.PixelIdxList= seed_pt_obj.PixelIdxList(ind_sort_obj); % reordered object list
	num_seed_obj(cumsum(num_seed_obj) > (2 * data_cat.num_trace)) ... 
							= 0; % limit cumulative number of seed pts to double number of traces
	
	for jj = find(num_seed_obj <= seed_pt_obj_cutoff) % remove small blobs or excessive seed pts
		seed_pt_mat(seed_pt_obj.PixelIdxList{jj}) ...
							= false;
	end
	
	[seed_pt_clean_i, seed_pt_clean_j] ...
							= find(seed_pt_mat);
	
	% sort seed pts by peak in wavelet transform
	seed_pt					= sortrows([pk_cwt(sub2ind([data_cat.num_sample data_cat.num_trace], seed_pt_clean_i, seed_pt_clean_j)) seed_pt_clean_i seed_pt_clean_j], 'descend');
	seed_pt					= seed_pt(:, 2:3);
	
	% propagate layers based on seed pts
	dt						= diff(data_cat.twtt(1:2)); % time interval, s
	[layer_ARESELP_mat, num_layer_ARESELP] ...
							= ARESELP_seed_layer(logical(pk_cwt), data_cat.num_sample, data_cat.num_trace, seed_pt, ...
												 ceil(twtt_sep_min / dt), ceil(len_blk / diff(data_cat.dist_lin(1:2))), ceil(twtt_jump_max / dt), angle_diff, ceil(twtt_proxim_max / dt), ceil(len_layer_min / dt)); 
	
	% generate layer structure to be saved
	layer_ARESELP			= struct;
	for jj = 1:num_layer_ARESELP
		layer_ARESELP(jj).ind_z	...
							= NaN(1, data_cat.num_trace);
		[ind_i_curr, ind_j_curr] ...
							= find(layer_ARESELP_mat == jj);
		layer_ARESELP(jj).ind_z(ind_j_curr) ...
							= ind_i_curr';
	end
	
	sum_sq_layer_len(ii)	= sum(sum(~isnan(reshape([layer_ARESELP(:).ind_z], data_cat.num_trace, num_layer_ARESELP)), 2) .^ 2); % sum of square of layer lengths (candidate ARESELP quality metric)
	disp(['Sum of square of layer lengths: ' num2str(sum_sq_layer_len(ii))])
	
	zscope((1e-3 .* data_cat.dist_lin), (1e6 .* data_cat.twtt), data_cat.amp)
	hold on
	for jj = 1:num_layer_ARESELP
		line((1e-3 .* data_cat.dist_lin(~isnan(layer_ARESELP(jj).ind_z))), (1e6 .* data_cat.twtt(layer_ARESELP(jj).ind_z(~isnan(layer_ARESELP(jj).ind_z)))), 'LineWidth', 2)
	end
	xlabel('Distance (km)')
	ylabel(['Traveltime (' char(181) 's)'])
	
	pause(1)
	
	% append dataset with ARESELP layers
    save([dir_in name_file{ii}], 'layer_ARESELP', 'num_layer_ARESELP', '-append')
end

disp('DONE processing these files.')