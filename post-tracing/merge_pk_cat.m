% MERGE_CAT_PK Merge picks into a single .mat file.
% 
% Joe MacGregor (NASA/GSFC)
% Last updated: 17 July 2024

clear

do_save						= true;
plotting					= false;

dir_in                      = '/Users/jamacgre/OneDrive - NASA/data/gris_strat2/arctic_cat_pk_rev_r1/';
dir_cat						= '/Users/jamacgre/OneDrive - NASA/data/gris_strat2/arctic_cat_r1/';
dir_mat						= '/Users/jamacgre/OneDrive - NASA/research/matlab/pickgui_v2/mat/';
dir_shp						= '/Users/jamacgre/OneDrive - NASA/research/matlab/pickgui_v2/shp/';
dist_int_decim              = 1e3; % decimation index distance, m
num_ind_layer_min			= 5; % minimum number of samples to define a layer
num_win_pk					= 1; % number of samples to search +/- for 
frac_std_adj				= 0.2; % maximum value of index adjustment beyond which assume it's a round vertical move for layer

speed_vacuum                = 299792458; % m/s
permitt_ice                 = 3.15; % dimensionless
speed_ice                   = speed_vacuum / sqrt(permitt_ice); % m/s

load([dir_mat 'xyz_all.mat'], 'campaign', 'num_campaign', 'num_segment', 'segment')
load('/Users/jamacgre/OneDrive - NASA/research/matlab/greenland/mat/grl_coast.mat', 'num_coast', 'x_coast', 'y_coast')

% BedMachine v5
BM5                         = struct;
BM5.x                       = double(ncread('/Users/jamacgre/OneDrive - NASA/research/data/greenland/BedMachine/BedMachineGreenland-v5.nc', 'x'))'; % projected x, m
BM5.y                       = double(flipud(ncread('/Users/jamacgre/OneDrive - NASA/research/data/greenland/BedMachine/BedMachineGreenland-v5.nc', 'y'))); % projected y, m
BM5.mask_gris               = double(rot90(ncread('/Users/jamacgre/OneDrive - NASA/research/data/greenland/BedMachine/BedMachineGreenland-v5.nc', 'mask'))); % ice mask

% simplified masks for whole Greenland maps
BM5.mask_combo              = BM5.mask_gris;
BM5.mask_gris(BM5.mask_gris ~= 2) ...
                            = 0;
BM5.mask_gris(BM5.mask_gris == 2) ...
                            = 1;
BM5.mask_gris               = logical(BM5.mask_gris); % now mask_gris is for ice sheet only
BM5.mask_combo(BM5.mask_combo == 3) ...
                            = 0;
BM5.mask_combo(BM5.mask_combo == 4) ...
                            = 1; % mask_combo includes 0/ocean 1/land 2/ice
BM5.mask_combo              = uint8(BM5.mask_combo); % reduce size

% Mouginot 2019 ice-sheet drainage polygon basins, with peripheral ice masses manually removed in QGIS, then rasterized to BM5 grid, took me forever to figure out but got it with SAGA in QGIS
M19_ice_mask				= readgeoraster('/Users/jamacgre/OneDrive - NASA/data/GreenValley/Greenland_Basins_PS_v1.4.2_simple_raster.tif');
M19_ice_mask				= flipud(logical(M19_ice_mask));

% now BM5 mask 2=peripheral ice masses and 3=ice sheet (previously 2)
BM5.mask_combo_plot			= BM5.mask_combo;
BM5.mask_combo((BM5.mask_combo == 2) & M19_ice_mask) ...
							= 3;
BM5.mask_gris(BM5.mask_combo == 2) ...
							= false;

% GrIMP DEM (Howat et al., 2022)
DEM							= struct;
[DEM.elev_surf, DEM.R]		= readgeoraster('/Users/jamacgre/OneDrive - NASA/data/GrIMP/GrIMP_100m.tif'); % surface elevation from a QGIS-merged/resampled GeoTIFF of the original 30-m tiles, m
DEM.elev_surf               = flipud(DEM.elev_surf);
DEM.elev_surf(DEM.elev_surf == -9999) ...
                            = NaN; % ugh GeoTIFF don't make me figure out your silly bad value, just start out as NaNs
[DEM.x, DEM.y]				= deal((DEM.R.XWorldLimits(1):DEM.R.CellExtentInWorldX:(DEM.R.XWorldLimits(2) - DEM.R.CellExtentInWorldX)), ((DEM.R.YWorldLimits(1) + DEM.R.CellExtentInWorldY):DEM.R.CellExtentInWorldY:DEM.R.YWorldLimits(2))');
DEM.elev_surf(interp2(BM5.x, BM5.y, BM5.mask_combo, DEM.x, DEM.y, 'nearest') ~= 3) ...
                            = NaN; % restrict DEM to ice sheet only

%%

% cycle through picks
disp(['Scanning ' dir_in '...'])

file_pk						= dir([dir_in '*_pk.mat']);
file_pk						= {file_pk.name};
num_file_pk					= length(file_pk);

[depth, dist, elev_wgs84, frame_cat, gps_time, ind_decim, ind_decim_mid, ind_trace_layer, ind_z, int, int_bed, int_surf, lat, lon, num_layer_max_decim, thick, thick_decim, twtt_ice, x, y] ...
							= deal(cell(1, num_file_pk));
[ind_campaign, ind_segment, num_decim, num_frame, num_layer, num_sample, num_trace, num_trace_layer] ...
							= deal(zeros(1, num_file_pk));
pk_empty					= false(1, num_file_pk);

%%
for ii = 1:num_file_pk
	
    disp([file_pk{ii}(1:(end - 4)) ' (' num2str(ii) ' / ' num2str(num_file_pk) ')...'])
    
	% current picks
	pk_curr					= load([dir_in file_pk{ii}]);
	pk_curr					= pk_curr.pk;
	
	% campaign and segment index assignment (usually obvious from filename but simplifies sorting later on)
	ind_campaign(ii)		= 1;
	while (ind_campaign(ii) <= num_campaign)
		if ~isempty(find(contains(segment{ind_campaign(ii)}, file_pk{ii}(6:16)), 1))
			ind_segment(ii) = find(contains(segment{ind_campaign(ii)}, file_pk{ii}(6:16)));
			disp(['Identified as segment ' segment{ind_campaign(ii)}{ind_segment(ii)} ' from campaign ' campaign{ind_campaign(ii)}])
			break
		else
			ind_campaign(ii)= ind_campaign(ii) + 1;
		end
	end
	
	% assign non-layer variables from pk file
	[dist{ii}, frame_cat{ii}, gps_time{ii}, int_bed{ii}, int_surf{ii}, lat{ii}, lon{ii}, num_frame(ii), num_layer(ii), num_sample(ii), num_trace(ii), thick{ii}, x{ii}, y{ii}] ...
							= deal(pk_curr.dist, pk_curr.file_in, pk_curr.time, pk_curr.int_bed, pk_curr.int_surf, pk_curr.lat, pk_curr.lon, length(pk_curr.file_in), pk_curr.num_layer, pk_curr.num_sample, pk_curr.num_trace, ...
								   ((pk_curr.twtt_bed - pk_curr.twtt_surf) .* (speed_ice / 2)), pk_curr.x, pk_curr.y);
	
	% load data file just need amplitude and traveltime
	data_cat_curr			= load([dir_cat file_pk{ii}(1:(end - 7))], 'amp', 'twtt');
	
	% update layer positions based on amplitude
	for jj = 1:num_layer(ii)
		pk_curr.layer(jj).ind_z_adj ...
							= pk_curr.layer(jj).ind_z;
		iter_count			= 0;
		while (iter_count <= 1)
			for kk = find(~isnan(pk_curr.layer(jj).ind_z_adj)) % loop through each sample of each layer and look around num_win_pk as to whether there's a better maximum
				[~, ind_z_adj] ...
							= max(data_cat_curr.amp((pk_curr.layer(jj).ind_z_adj(kk) - num_win_pk):(pk_curr.layer(jj).ind_z_adj(kk) + num_win_pk), kk)); % i/z index of nearest max
				pk_curr.layer(jj).ind_z_adj(kk) ...
							= pk_curr.layer(jj).ind_z_adj(kk) + (ind_z_adj - (num_win_pk + 1)); % correct y index of max because search was done in a narrow window
			end
			ind_diff_mean	= mean((pk_curr.layer(jj).ind_z_adj - pk_curr.layer(jj).ind_z), 'omitnan');
			if (((mod(ind_diff_mean, 1) <= frac_std_adj) || (mod(ind_diff_mean, 1) >= (1 - frac_std_adj))) && (abs(round(ind_diff_mean)) > 0)) % close to unity
				pk_curr.layer(jj).ind_z_adj ...
							= pk_curr.layer(jj).ind_z + round(ind_diff_mean);
				disp(['Layer #' num2str(jj) ' moved by ' num2str(round(ind_diff_mean))])
				iter_count	= iter_count + 1;
			else
				break
			end
		end
		[pk_curr.layer(jj).int_adj, pk_curr.layer(jj).twtt_adj] ...
							= deal(NaN(1, num_trace(ii)));
		pk_curr.layer(jj).int_adj(~isnan(pk_curr.layer(jj).ind_z_adj)) ...
							= data_cat_curr.amp(sub2ind([num_sample(ii) num_trace(ii)], pk_curr.layer(jj).ind_z_adj(~isnan(pk_curr.layer(jj).ind_z_adj)), find(~isnan(pk_curr.layer(jj).ind_z_adj)))); % adjusted echo intensity 
        pk_curr.layer(jj).twtt_adj(~isnan(pk_curr.layer(jj).ind_z_adj)) ...
                            = data_cat_curr.twtt(pk_curr.layer(jj).ind_z_adj(~isnan(pk_curr.layer(jj).ind_z_adj)))';
        pk_curr.layer(jj).twtt_ice_adj ...
                            = pk_curr.layer(jj).twtt_adj - pk_curr.twtt_surf;
        pk_curr.layer(jj).depth_adj ...
                            = pk_curr.layer(jj).twtt_ice_adj .* (speed_ice / 2);
	end

	% ind_diff_all			= reshape([pk_curr.layer(:).ind_z_adj], num_trace(ii), num_layer(ii))' - reshape([pk_curr.layer(:).ind_z], num_trace(ii), num_layer(ii))';
	% disp(['Layer adjustment: ' sprintf('%1.2f', mean(ind_diff_all, 'all', 'omitnan')) ' ± ' sprintf('%1.2f', std(ind_diff_all, 0, 'all', 'omitnan'))])
	
	% matrix-ify layer variables
	ind_z{ii}				= reshape([pk_curr.layer(:).ind_z_adj], num_trace(ii), num_layer(ii))';
	
	% check for overlapping layers post-adjustment (bad), following pickgui's pk_cross function
	overlap_check			= false;
	overlap_set				= [];
	if (num_layer(ii) > 1)
		tmp2				= find(sum(~isnan(ind_z{ii})) > 1); % traces with >1 layer
		tmp1				= ind_z{ii}(:, tmp2);
		for jj = 1:(num_layer(ii) - 1)
        	tmp3			= find(sum(tmp1((jj + 1):end, ~isnan(tmp1(jj, :))), 2, 'omitnan')); % skip this layer if subsequent layers have no overlap where this layer exists
        	if isempty(tmp3)
            	continue
        	end
        	for kk = (jj + tmp3')
            	tmp4		= intersect(find(~isnan(tmp1(jj, :))), find(~isnan(tmp1(kk, :))));
            	if (~isempty(find(diff(sign(tmp1(jj, tmp4) - tmp1(kk, tmp4))), 1)) || all(diff(tmp1([jj kk], tmp4), 1, 1) == 0))
					overlap_check ...
							= true;
					overlap_set ...
							= [overlap_set; jj kk max([0 round(1e-3 .* pk_curr.dist_lin(tmp2(tmp4(find(diff(sign(tmp1(jj, tmp4) - tmp1(kk, tmp4))), 1)))))])]; %#ok<AGROW>
            	end
        	end
		end
	end
	if overlap_check
		disp(overlap_set)		
		error('OVERLAP!')
	end

	depth{ii}				= reshape([pk_curr.layer(:).depth_adj], num_trace(ii), num_layer(ii))';
	elev_wgs84{ii}			= repmat(double(interp2(DEM.x, DEM.y, DEM.elev_surf, x{ii}, y{ii})), num_layer(ii), 1) - depth{ii}; % referenced to GrIMP
	twtt_ice{ii}			= reshape([pk_curr.layer(:).twtt_ice_adj], num_trace(ii), num_layer(ii))';
	int{ii}					= reshape([pk_curr.layer(:).int_adj], num_trace(ii), num_layer(ii))';
	
	% identify small layers (<num_ind_layer_min), of which there should be none now
	len_layer				= sum(~isnan(depth{ii}), 2);
	if ~isempty(find((len_layer < num_ind_layer_min), 1))
		ind_layer_long		= find(len_layer >= num_ind_layer_min)';
		error('SHORT LAYER(S)!')
	end
	
	% decimation setup for easier plotting/analysis later on
    dist_decim_tmp			= dist{ii}(1):dist_int_decim:dist{ii}(end);
    dist_decim_tmp(end)		= dist{ii}(end);
	ind_decim{ii}			= unique(interp1(unique(dist{ii}), 1:length(unique(dist{ii})), dist_decim_tmp, 'nearest', 'extrap'));
	ind_decim_mid{ii}		= round(ind_decim{ii}(1:(end - 1)) + (diff(ind_decim{ii}) ./ 2));
	num_decim(ii)			= length(ind_decim{ii}) - 1;
	
	% maximum number of layers and thickness in each decimated segment
	[num_layer_max_decim{ii}, thick_decim{ii}] ...
							= deal(NaN(1, num_decim(ii)));
	for jj = 1:num_decim(ii)
		num_layer_max_decim{ii}(jj) ...
							= max(sum(~isnan(depth{ii}(:, ind_decim{ii}(jj):(ind_decim{ii}(jj + 1) - 1)))));
		thick_decim{ii}(jj) = mean(thick{ii}(ind_decim{ii}(jj):(ind_decim{ii}(jj + 1) - 1)), 'omitnan');
	end
	
	% trim layer matrices to only traces with layers (sometimes a big volume reduction)
	ind_trace_layer{ii}		= find(sum(~isnan(depth{ii}), 1));
	if (length(ind_trace_layer{ii}) < num_trace(ii))
		[ind_z{ii}, depth{ii}, elev_wgs84{ii}, int{ii}, twtt_ice{ii}] ...
							= deal(ind_z{ii}(:, ind_trace_layer{ii}), depth{ii}(:, ind_trace_layer{ii}), elev_wgs84{ii}(:, ind_trace_layer{ii}), int{ii}(:, ind_trace_layer{ii}), twtt_ice{ii}(:, ind_trace_layer{ii}));
	end
	
	% sparsify layer variables
	ind_nan					= find(isnan(depth{ii}));
	[ind_z{ii}(ind_nan), depth{ii}(ind_nan), elev_wgs84{ii}(ind_nan), int{ii}(ind_nan), twtt_ice{ii}(ind_nan)] ...
							= deal(0); % sparse needs zeros
	[ind_z{ii}, depth{ii}, elev_wgs84{ii}, int{ii}, twtt_ice{ii}] ...
							= deal(sparse(ind_z{ii}), sparse(depth{ii}), sparse(elev_wgs84{ii}), sparse(int{ii}), sparse(twtt_ice{ii}));
	
	num_trace_layer(ii)		= length(ind_trace_layer{ii});
	
	if isempty(depth{ii}) % uh oh no layers left (now all fixed)
		pk_empty(ii)		= true;
	end
end

%%
% fix 2011 ambiguities (ugh never fly two different aircraft/systems on the same day)
ind_campaign(187)			= 19;
for ii = 187
	ind_segment(ii)			= find(contains(segment{ind_campaign(ii)}, file_pk{ii}(6:16)));
end

ind_pk_merge				= [ind_campaign' ind_segment' ones(length(ind_campaign), 1)]; % merge+traced set
for ii = 2:size(ind_pk_merge, 1)
	if isequal(ind_pk_merge((ii - 1), 1:2), ind_pk_merge(ii, 1:2))
		ind_pk_merge(ii, 3)	= ind_pk_merge((ii - 1), 3) + 1;
	end
end
ind_subsegment				= ind_pk_merge(:, 3);

if do_save
	save([dir_mat 'pk_cat.mat'], 'depth', 'dist', 'elev_wgs84', 'frame_cat', 'gps_time', 'file_pk', 'ind_campaign', 'ind_decim', 'ind_decim_mid', 'ind_pk_merge', 'ind_segment', 'ind_subsegment', 'ind_trace_layer', ...
								 'ind_z', 'int', 'int_bed', 'int_surf', 'lat', 'lon', 'num_decim', 'num_file_pk', 'num_frame', 'num_layer', 'num_layer_max_decim', 'num_sample', 'num_trace', 'num_trace_layer', 'pk_empty', 'thick', ...
								 'thick_decim', 'twtt_ice', 'x', 'y')
    disp(['Saved merged picks in ' dir_mat ' as pk_cat.mat.'])
end

%%
if plotting
%%
    set(0, 'DefaultFigureWindowStyle', 'default') %#ok<UNRCH>

%%
    set(0, 'DefaultFigureWindowStyle', 'docked')

%%
end