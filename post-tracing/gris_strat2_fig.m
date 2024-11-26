% GRIS_STRAT2_FIG Figures for Greenland radiostratigraphy v2 manuscript.
% 
% Joe MacGregor (NASA)
% Last updated: 26 November 2024

clear

plotting                    = false;

dir_mat						= '/Users/jamacgre/OneDrive - NASA/research/matlab/pickgui_v2/mat/';
dir_gl						= '/Users/jamacgre/OneDrive - NASA/research/matlab/greenland/mat/';

load([dir_mat 'xyz_all.mat'], 'campaign', 'num_campaign', 'num_segment', 'segment', 'x', 'y')
pk_cat						= load([dir_mat 'pk_cat.mat'], 'depth', 'file_pk', 'ind_campaign', 'ind_decim', 'ind_decim_mid', 'ind_trace_layer', 'num_decim', 'num_file_pk', 'num_layer', 'num_trace', 'thick', 'x', 'y');
load([dir_gl 'grl_coast.mat'], 'num_coast', 'x_coast', 'y_coast')
load([dir_mat 'core_int.mat'], 'name_core', 'num_core', 'x_core', 'y_core')
load([dir_mat 'date_all.mat'], 'age', 'age_type', 'age_uncert')
load([dir_mat 'age_grd2_clean.mat'], 'age_iso', 'age_norm_smooth', 'age_norm_uncert_tot_smooth', 'depth_iso_smooth', 'depth_iso_uncert_tot_smooth', 'depth_norm', 'num_age_iso', 'num_depth_norm', 'x_grd', 'y_grd')
[xx_grd, yy_grd]			= meshgrid(x_grd, y_grd);

% v1 stratigraphy
age_grd2_v1					= load([dir_gl 'age_grd2_krige.mat'], 'age_iso', 'depth_iso2', 'depth_norm', 'num_age_iso', 'num_depth_norm', 'x_grd', 'y_grd');
age_grd2_v1_alt				= load([dir_gl 'age_grd2_krige_alt.mat'], 'age_norm2_alt');
age_grd2_v1.age_norm2		= age_grd2_v1_alt.age_norm2_alt;
clear age_grd2_v1_alt

[x_min, x_max, y_min, y_max]= deal(-632e3, 846e3, -3344e3, -670e3);
letters						= 'a':'z';

[x_core_label, y_core_label]= deal((1e3 .* [-560 -40 155 25 279 -95 -230]), (1e3 .* [-1260 -2790 -1475 -1870 -1856 -1330 -1615]));

FDM							= load([dir_mat 'gsfc_fdm_121.mat'], 'FAC', 'time_dt', 'x', 'y');
FDM.num_time_dt				= length(FDM.time_dt);
[FDM.x, FDM.y]				= deal(double(FDM.x), double(FDM.y));

% BedMachine v5
BM5                         = struct;
BM5.x                       = double(ncread('/Users/jamacgre/OneDrive - NASA/research/data/greenland/BedMachine/BedMachineGreenland-v5.nc', 'x'))'; % projected x, m
BM5.y                       = double(flipud(ncread('/Users/jamacgre/OneDrive - NASA/research/data/greenland/BedMachine/BedMachineGreenland-v5.nc', 'y'))); % projected y, m
BM5.mask_gris               = double(rot90(ncread('/Users/jamacgre/OneDrive - NASA/research/data/greenland/BedMachine/BedMachineGreenland-v5.nc', 'mask'))); % ice mask
BM5.elev_surf				= double(rot90(ncread('/Users/jamacgre/OneDrive - NASA/research/data/greenland/BedMachine/BedMachineGreenland-v5.nc', 'surface'))); % surface elevation, m
BM5.thick					= double(rot90(ncread('/Users/jamacgre/OneDrive - NASA/research/data/greenland/BedMachine/BedMachineGreenland-v5.nc', 'thickness'))); % ice thickness, m

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

BM5.mask_gris_grd			= interp2(BM5.x, BM5.y, BM5.mask_gris, xx_grd, yy_grd, 'nearest');
BM5.mask_combo_plot_grd		= interp2(BM5.x, BM5.y, BM5.mask_combo_plot, xx_grd, yy_grd, 'nearest');

% Mouginot et al. (2019) large-scale drainage basins
M19                         = shaperead('/Users/jamacgre/OneDrive - NASA/data/GBaTSv2/Greenland_Basins_PS_v1.4.2_simple.shp');
num_M19                     = length(M19);

[paral, merid, x_paral, y_paral, x_merid, y_merid] ...
                            = graticule_greenland([x_min x_max], [y_min y_max], 5, 10);

age_type_lim				= [0 1 10 Inf];

% simplify filenames and unsparsify depth (much simpler), and correct for firn
[depth_pk, depth_norm_pk, num_layer_pk, num_layer_type, FDM_corr] ...
							= deal(cell(1, pk_cat.num_file_pk));
for ii = 1:pk_cat.num_file_pk
	depth_tmp				= full(pk_cat.depth{ii});
	depth_tmp(~depth_tmp)	= NaN;
	depth_pk{ii}			= NaN(pk_cat.num_layer(ii), pk_cat.num_trace(ii));
	depth_pk{ii}(:, pk_cat.ind_trace_layer{ii}) ...
							= depth_tmp;
	num_layer_pk{ii}		= NaN(1, pk_cat.num_decim(ii));
	num_layer_type{ii}		= zeros(length(age_type_lim), pk_cat.num_decim(ii));
	depth_pk_tmp			= ~isnan(depth_pk{ii});
	ind_age_type			= cell(length(age_type_lim), 1);
	for jj = 1:length(age_type_lim)
		if (jj == 1)
			ind_age_type{jj}= find(age_type{ii} == age_type_lim(jj));
		elseif (jj == length(age_type_lim))
			ind_age_type{jj}= find(isnan(age_type{ii}));
		else
			ind_age_type{jj}= find((age_type{ii} > age_type_lim(jj - 1)) & (age_type{ii} <= age_type_lim(jj)));
		end
	end
	for jj = 1:pk_cat.num_decim(ii)
		num_layer_pk{ii}(jj)= length(find(sum(depth_pk_tmp(:, pk_cat.ind_decim{ii}(jj):pk_cat.ind_decim{ii}(jj + 1)), 2)));
		for kk = find(cellfun(@numel, ind_age_type))'
			num_layer_type{ii}(kk, jj) ...
							= length(find(sum(depth_pk_tmp(ind_age_type{kk}, pk_cat.ind_decim{ii}(jj):pk_cat.ind_decim{ii}(jj + 1)), 2)));
		end
	end
	ind_FDM_curr			= interp1(FDM.time_dt, 1:FDM.num_time_dt, datetime(str2double(pk_cat.file_pk{ii}(6:13)), 'ConvertFrom', 'yyyymmdd'), 'nearest'); % time index for flight in FDM time vector
	FDM_corr{ii}			= interp2(FDM.x, FDM.y, squeeze(FDM.FAC(:, :, ind_FDM_curr)), pk_cat.x{ii}(ones(pk_cat.num_layer(ii), 1), :), pk_cat.y{ii}(ones(pk_cat.num_layer(ii), 1), :), 'linear', 0);
	depth_pk{ii}			= depth_pk{ii} + FDM_corr{ii}; % add correction for firn air content (FAC)
	FDM_corr{ii}			= FDM_corr{ii}(1, :);
	depth_norm_pk{ii}		= depth_pk{ii} ./ pk_cat.thick{ii};
end
clear depth

campaign_priority			= [1 1 1 2 2 2 2 2 2 1 1 2 NaN 1 1 2 2 3 1 3 2 3 2 NaN 2 NaN NaN 3 2 2];
ind_campaign_ignore			= [13 24 26 27]; % campaigns to ignore (not NASA or NSF)

% clean-up v1
age_grd2_v1.depth_iso2		= age_grd2_v1.depth_iso2 .* repmat(interp2(BM5.x, BM5.y, BM5.thick, (1e3 .* age_grd2_v1.x_grd), (1e3 .* age_grd2_v1.y_grd)), 1, 1, age_grd2_v1.num_age_iso); % dimensionalize
[age_grd2_v1.age_norm2(age_grd2_v1.age_norm2 == 0), age_grd2_v1.depth_iso2(age_grd2_v1.depth_iso2 == 0)] ...
                            = deal(NaN);

%%

DTU							= shaperead('/Users/jamacgre/OneDrive - NASA/research/funding/greenland_layers_v2/misc/all_radar_flightlines_map_ai_20200708/shapefile/all_radar_flightlines_map_geo_vec_ai_01.shp');

path_AWI					= '/Users/jamacgre/OneDrive - NASA/research/funding/greenland_layers_v2/misc/awi_lines/shp/';
dir_AWI						= dir([path_AWI '*.shp']);
dir_AWI						= {dir_AWI.name};
num_AWI						= length(dir_AWI);
[AWI, x_AWI, y_AWI]			= deal(cell(1, num_AWI));
for ii = 1:num_AWI
	AWI{ii}					= shaperead([path_AWI dir_AWI{ii}]);
		[x_AWI{ii}, y_AWI{ii}] ...
							= projfwd(projcrs(3413), [AWI{ii}(:).Y], [AWI{ii}.X]);
end

path_hiawatha				= '/Users/jamacgre/OneDrive - NASA/research/matlab/hiawatha/shapefile/';
dir_hiawatha				= dir([path_hiawatha 'awi_day*.shp']);
dir_hiawatha				= {dir_hiawatha.name};
num_hiawatha				= length(dir_hiawatha);
[hiawatha, x_hiawatha, y_hiawatha] ...
							= deal(cell(1, num_hiawatha));
for ii = 1:num_hiawatha
	hiawatha{ii}			= shaperead([path_hiawatha dir_hiawatha{ii}]);
		[x_hiawatha{ii}, y_hiawatha{ii}] ...
							= projfwd(projcrs(3413), [hiawatha{ii}(:).Y], [hiawatha{ii}.X]);
end

%%

data_cat					= load('/Users/jamacgre/OneDrive - NASA/data/gris_strat2/arctic_cat_r1/Data_20170413_01_035-058.mat');
ind_pk_cat					= find(contains(pk_cat.file_pk, 'Data_20170413_01_035-058'));
data_cat_extra				= load('/Users/jamacgre/OneDrive - NASA/research/matlab/pickgui_v2/mat/Data_20170413_01_035-058_pk_plots.mat');

depth_norm_plot				= linspace(0, 1, data_cat.num_sample);

speed_vacuum                = 299792458; % m/s
permitt_ice                 = 3.15; % dimensionless
speed_ice                   = speed_vacuum / sqrt(permitt_ice); % m/s

dt							= median(diff(data_cat.twtt)); % traveltime sampling interval, s
depth_vec					= (speed_ice / 2) .* (0:dt:((data_cat.num_sample - 1) * dt))'; % simple monotonically increasing depth vector, m

num_color_db				= 2 ^ 8; % number of colors for display

% frames 49 to 56
ind_range					= [46691 73360];
[db_mean, db_std]			= deal(mean(data_cat.amp(:, ind_range(1):ind_range(2)), 'all'), std(data_cat.amp(:, ind_range(1):ind_range(2)), [], 'all'));

% topographically correct data and flip
elev_surf_grd_curr			= interp2(BM5.x, BM5.y, BM5.elev_surf, data_cat.x, data_cat.y);
elev_vec					= flipud(max(elev_surf_grd_curr) - depth_vec); % elevation vector, m
amp_elev					= single(flipud(topocorr(data_cat_extra.amp_depth, depth_vec, elev_surf_grd_curr, true)));

BM5_mask_color				= cell(1, 3);
color_ocean					= [135 206 235] ./ 255;
color_land					= [159 89 39] ./ 255;
for ii = 1:3
	BM5_mask_color{ii}		= 0.9 .* ones(length(BM5.y(1:10:end)), length(BM5.x(1:10:end)));
	BM5_mask_color{ii}((BM5.mask_combo_plot(1:10:end, 1:10:end) == 0)) ...
							= color_ocean(ii);
	BM5_mask_color{ii}((BM5.mask_combo_plot(1:10:end, 1:10:end) == 1)) ...
							= color_land(ii);
end
BM5_mask_color				= cat(3, BM5_mask_color{1}, BM5_mask_color{2}, BM5_mask_color{3});

%%

ind_pk_core					= find(contains(pk_cat.file_pk, '20140410_01_006-042'));
data_cat_core				= load('/Users/jamacgre/OneDrive - NASA/data/gris_strat2/arctic_cat_r1/Data_20140410_01_006-042.mat');
data_cat_core_extra			= load('/Users/jamacgre/OneDrive - NASA/research/matlab/pickgui_v2/mat/Data_20140410_01_006-042_pk_plots.mat');
dist_range_core				= [1100e3 data_cat_core.dist_lin(end)];
dt_core						= median(diff(data_cat_core.twtt)); % traveltime sampling interval, s
depth_vec_core				= (speed_ice / 2) .* (0:dt_core:((data_cat_core.num_sample - 1) * dt_core))'; % simple monotonically increasing depth vector, m

% do some basic amplitude corrections to amp_depth to get scaling ok, best to do before topographic correction
data_cat_core_extra.amp_depth ...
							= data_cat_core_extra.amp_depth + (3e-2 .* depth_vec_core(:, ones(1, data_cat_core.num_trace))); % correct for two-way attenuation (assume 15 dB/km one-way)
data_cat_core_extra.amp_depth ...
							= data_cat_core_extra.amp_depth + repmat((2e1 .* log10(data_cat_core.elev_air - data_cat_core.elev_surf)), data_cat_core.num_sample, 1); % correct for geometric spreading to ice surface
data_cat_core_extra.amp_depth(2:end, :) ...
							= data_cat_core_extra.amp_depth(2:end, :) + (2e1 .* log10(depth_vec_core(2:end, ones(1, data_cat_core.num_trace)) ./ sqrt(permitt_ice))); % correct for geometric spreading beneath ice surface, skip depth=0
if ~isreal(data_cat_core_extra.amp_depth)
	data_cat_core_extra.amp_depth ...
							= real(data_cat_core_extra.amp_depth);
end
data_cat_core_extra.amp_depth(isinf(data_cat_core_extra.amp_depth))	...
							= NaN; % should no longer be an issue but just in case
% topographically correct data and flip
elev_surf_grd_core			= interp2(BM5.x, BM5.y, BM5.elev_surf, data_cat_core.x, data_cat_core.y);
elev_vec_core				= flipud(max(elev_surf_grd_core) - depth_vec_core); % elevation vector, m
amp_elev_core				= single(flipud(topocorr(data_cat_core_extra.amp_depth, depth_vec_core, elev_surf_grd_core, true)));

%%

ind_pk_norm					= find(contains(pk_cat.file_pk, '20190418_01_007-007'));
data_cat_norm				= load('/Users/jamacgre/OneDrive - NASA/data/gris_strat2/arctic_cat_r1/Data_20190418_01_007-007.mat');
data_cat_norm_extra			= load('/Users/jamacgre/OneDrive - NASA/research/matlab/pickgui_v2/mat/Data_20190418_01_007-007_pk_plots.mat');
% dist_range_norm				= [1100e3 data_cat_norm.dist_lin(end)];
dt_norm						= median(diff(data_cat_norm.twtt)); % traveltime sampling interval, s
depth_vec_norm				= (speed_ice / 2) .* (0:dt_norm:((data_cat_norm.num_sample - 1) * dt_norm))'; % simple monotonically increasing depth vector, m
% % do some basic amplitude corrections to amp_depth to get scaling ok, best to do before topographic correction
% data_cat_norm_extra.amp_depth ...
% 							= data_cat_norm_extra.amp_depth + (3e-2 .* depth_vec_norm(:, ones(1, data_cat_norm.num_trace))); % correct for two-way attenuation (assume 15 dB/km one-way)
% data_cat_norm_extra.amp_depth ...
% 							= data_cat_norm_extra.amp_depth + repmat((2e1 .* log10(data_cat_norm.elev_air - data_cat_norm.elev_surf)), data_cat_norm.num_sample, 1); % correct for geometric spreading to ice surface
% data_cat_norm_extra.amp_depth(2:end, :) ...
% 							= data_cat_norm_extra.amp_depth(2:end, :) + (2e1 .* log10(depth_vec_norm(2:end, ones(1, data_cat_norm.num_trace)) ./ sqrt(permitt_ice))); % correct for geometric spreading beneath ice surface, skip depth=0
% if ~isreal(data_cat_norm_extra.amp_depth)
% 	data_cat_norm_extra.amp_depth ...
% 							= real(data_cat_norm_extra.amp_depth);
% end
% data_cat_norm_extra.amp_depth(isinf(data_cat_norm_extra.amp_depth))	...
% 							= NaN; % should no longer be an issue but just in case
% topographically correct data and flip
elev_surf_grd_norm			= interp2(BM5.x, BM5.y, BM5.elev_surf, data_cat_norm.x, data_cat_norm.y);
elev_vec_norm				= flipud(max(elev_surf_grd_norm) - depth_vec_norm); % elevation vector, m
amp_elev_norm				= single(flipud(topocorr(data_cat_norm_extra.amp_depth, depth_vec_norm, elev_surf_grd_norm, true)));

%%
if plotting
%%
    set(0, 'DefaultFigureWindowStyle', 'default') %#ok<UNRCH>

%%
    set(0, 'DefaultFigureWindowStyle', 'docked')

%% INTRO MAPS

    color_fig				= [([135 206 235] ./ 255); ([159 89 39] ./ 255); 0.9 0.9 0.9; parula(6)];
    figure('Position', [100 100 1200 1000], 'Color', 'w')
    colormap(color_fig)	
	ax						= gobjects(1, 2);
	color_priority			= [247 252 185; 173 221 142; 49 163 84] ./ 255;%{'r' 'b' 'g'};
	titles					= {'All campaigns by priority rating' 'Number of reflections traced'};
	for ii = 1:2
    	ax(ii)              = subplot('Position', [((ii - 1) * 0.50) 0.03 0.48 0.91]);
    	hold on
    	axis equal
    	axis([x_min x_max y_min y_max])
    	imagesc(BM5.x, BM5.y, BM5.mask_combo_plot)
    	for jj = 1:length(x_paral)
        	line(x_paral{jj}, y_paral{jj}, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
    	end
    	for jj = 1:length(x_merid)
        	line(x_merid{jj}, y_merid{jj}, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
    	end
		for jj = 1:num_coast
        	line((1e3 .* x_coast{jj}), (1e3 .* y_coast{jj}), 'Color', 'k', 'linewidth', 1)
		end
		switch ii
			case 1
				[~, ind_priority_sort] ...
							= sort(campaign_priority);
				for jj = setdiff(ind_priority_sort, ind_campaign_ignore, 'stable')
        			for kk = 1:num_segment(jj)
						line(x{jj}{kk}(1:20:end), y{jj}{kk}(1:20:end), 'LineWidth', 2, 'Color', color_priority(campaign_priority(jj), :))
        			end
				end
				pl			= gobjects(1, 3);
				for jj = 1:3
					pl(jj)	= line(NaN, NaN, 'LineWidth', 3, 'Color', color_priority(jj, :));
				end
			case 2
				for jj = setdiff(1:num_campaign, ind_campaign_ignore)
        			for kk = 1:num_segment(jj)
						line(x{jj}{kk}(1:20:end), y{jj}{kk}(1:20:end), 'LineWidth', 1, 'Color', 'k')
        			end
				end
				[x_cat, y_cat, num_layer_cat] ...
							= deal([]);
				for jj = 1:pk_cat.num_file_pk
					ind_good = find(num_layer_pk{jj} > 0);
					[x_cat, y_cat, num_layer_cat] ...
							= deal([x_cat pk_cat.x{jj}(pk_cat.ind_decim_mid{jj}(ind_good))], [y_cat pk_cat.y{jj}(pk_cat.ind_decim_mid{jj}(ind_good))], [num_layer_cat num_layer_pk{jj}(ind_good)]);
				end
				[~, ind_sort] ...
							= sort(num_layer_cat);
				scatter(x_cat(ind_sort), y_cat(ind_sort), 10, color_fig(interp1(((0:5:25) + (diff(0:5:30) ./ 2)), 4:9, num_layer_cat(ind_sort), 'nearest', 'extrap'), :), 'filled')
				pb			= line(NaN, NaN, 'LineWidth', 1, 'Color', 'k');
		end
		pic					= line(x_core, y_core, 'Color', 'k', 'Marker', '^', 'MarkerSize', 10, 'MarkerFaceColor', 'm', 'LineStyle', 'none', 'LineWidth', 1);
    	for jj = 1:num_core
			if (ii == 1)
				text(x_core_label(jj), y_core_label(jj), name_core{jj}, 'Color', 'm', 'FontSize', 18, 'FontWeight', 'bold', 'EdgeColor', 'k', 'BackgroundColor', 'w')
			end
    	end
    	ax(ii).FontSize		= 20;
		ax(ii).XTick		= [];
		ax(ii).YTick		= [];
		clim([0 9])
    	fill((1e3 .* [325 425 425 325]), (1e3 .* [-3230 -3230 -3260 -3260]), 'k')
    	fill((1e3 .* [425 525 525 425]), (1e3 .* [-3230 -3230 -3260 -3260]), 'w', 'EdgeColor', 'k')
    	text(300e3, -3180e3, '0', 'Color', 'k', 'FontSize', 18)
    	text(520e3, -3180e3, '200 km', 'Color', 'k', 'FontSize', 18)
    	title(['(' letters(ii) ') ' titles{ii}], 'Color', 'k', 'FontSize', 20, 'FontWeight', 'bold')
    	text(-542e3, -3238e3, '60\circN', 'Color', 'k', 'FontSize', 18, 'Rotation', -10)
    	text(-330e3, -3250e3, '50\circW', 'Color', 'k', 'FontSize', 18, 'Rotation', 82)
    	text(552e3, -2750e3, '65\circN', 'Color', 'k', 'FontSize', 18, 'Rotation', 15)
    	text(750e3, -2975e3, '30\circW', 'Color', 'k', 'FontSize', 18, 'Rotation', -75)
		switch ii
			case 1
				legend([fliplr(pl) pic], [fliplr({'low' 'medium' 'high'}) 'deep ice core'], 'FontSize', 20, 'Location', 'northeast')
			case 2
				legend(pb, 'none traced', 'FontSize', 20, 'Location', 'northeast')
				colorbar('FontSize', 20, 'Limits', [3 9], 'YTick', 3:9, 'YTickLabel', {'1' '5' '10' '15' '20' '25' '>30'}, 'TickLength', 0.0175)
		end
    	grid off
    	box on
	end

%% RADARGRAM CORNUCOPIA

	num_color_db			= 2 ^ 10;
	[db_min, db_max]		= deal(-180, -60);	
    figure('Position', [100 100 1300 1330], 'Color', 'w')
    colormap(bone(num_color_db))
	y_range					= [0 3000; 0 3000; 0 1; 0 1; -400 2600];
	ylabels					= {'Depth (m)' 'Depth (m)' {'Thickness-normalized'; 'depth'} {'Reflection-normalized'; 'depth'} 'Elevation (m)'};
	ax						= gobjects(1, 5);
    ind_surf				= interp1(data_cat.twtt, 1:data_cat.num_sample, data_cat.twtt_surf, 'nearest', 'extrap');
    ind_bed					= interp1(data_cat.twtt, 1:data_cat.num_sample, data_cat.twtt_bed, 'nearest', 'extrap');
	for ii = 1:5
    	ax(ii)              = subplot('Position', [0.07 (0.815 - (0.19 * (ii - 1))) 0.92 0.175]);
    	hold on
		switch ii
			case 1
				imagesc((1e-3 .* (data_cat.dist_lin - data_cat.dist_lin(ind_range(1)))), depth_vec, data_cat_extra.amp_depth, [db_min db_max])
				for jj = 1:data_cat.num_layer_ARESELP
					ind_good = find(~isnan(data_cat.layer_ARESELP(jj).ind_z));
					line((1e-3 .* (data_cat.dist_lin(ind_good) - data_cat.dist_lin(ind_range(1)))), depth_vec(data_cat.layer_ARESELP(jj).ind_z(ind_good) - ind_surf(ind_good) + 1), 'Color', 'c', 'LineWidth', 1, 'LineStyle', '--')
				end
				line((1e-3 .* (data_cat.dist_lin - data_cat.dist_lin(ind_range(1)))), depth_vec(ind_bed - ind_surf + 1), 'Color', 'r', 'LineWidth', 2)
			case 2
				imagesc((1e-3 .* (data_cat.dist_lin - data_cat.dist_lin(ind_range(1)))), depth_vec, data_cat_extra.amp_depth, [db_min db_max])
				for jj = 1:pk_cat.num_layer(ind_pk_cat)
					ind_good = find(~isnan(depth_pk{ind_pk_cat}(jj, :)));
					line((1e-3 .* (data_cat.dist_lin(ind_good) - data_cat.dist_lin(ind_range(1)))), (depth_pk{ind_pk_cat}(jj, ind_good) + FDM_corr{ind_pk_cat}(ind_good)), 'Marker', '.', 'Color', 'b', 'MarkerSize', 6, 'LineStyle', 'none')
				end
			case 3
				imagesc((1e-3 .* (data_cat.dist_lin - data_cat.dist_lin(ind_range(1)))), depth_norm_plot, data_cat_extra.amp_norm, [db_min db_max])
				% for jj = 1:pk_cat.num_layer(ind_pk_cat)
				% 	ind_good = find(~isnan(data_cat_extra.pk_ind_z_norm(jj, :)));
				% 	line((1e-3 .* (data_cat.dist_lin(ind_good) - data_cat.dist_lin(ind_range(1)))), depth_norm(data_cat_extra.pk_ind_z_norm(jj, ind_good)), 'Marker', '.', 'Color', 'b', 'MarkerSize', 6, 'LineStyle', 'none')
				% end
				line(0, 1, 'Marker', 'o', 'MarkerFaceColor', 'g', 'MarkerSize', 10, 'Color', 'k', 'LineWidth', 1)
				line((1e-3 * (data_cat.dist_lin(ind_range(2)) - data_cat.dist_lin(ind_range(1)))), 1, 'Marker', '^', 'MarkerFaceColor', 'r', 'MarkerSize', 10, 'Color', 'k', 'LineWidth', 1)
			case 4
				imagesc((1e-3 .* (data_cat.dist_lin - data_cat.dist_lin(ind_range(1)))), depth_norm_plot, data_cat_extra.amp_flat, [db_min db_max])
				for jj = 1:pk_cat.num_layer(ind_pk_cat)
					ind_good = find(~isnan(data_cat_extra.pk_ind_z_flat(jj, :)));
					line((1e-3 .* (data_cat.dist_lin(ind_good) - data_cat.dist_lin(ind_range(1)))), depth_norm_plot(data_cat_extra.pk_ind_z_flat(jj, ind_good)), 'Marker', '.', 'Color', 'b', 'MarkerSize', 6, 'LineStyle', 'none')
				end
			case 5
				imagesc((1e-3 .* (data_cat.dist_lin - data_cat.dist_lin(ind_range(1)))), elev_vec, amp_elev, [db_min db_max])
				% for jj = 1:pk_cat.num_layer(ind_pk_cat)
				% 	ind_good = find(~isnan(depth_pk{ind_pk_cat}(jj, :)));
				% 	line((1e-3 .* (data_cat.dist_lin(ind_good) - data_cat.dist_lin(ind_range(1)))), (elev_surf_grd_curr(ind_good) - depth_pk{ind_pk_cat}(jj, ind_good)), 'Marker', '.', 'Color', 'b', 'MarkerSize', 6, 'LineStyle', 'none')
				% end
		end
		switch ii
			case 1
				pl			= gobjects(1, 2);
				pl(1)		= line(NaN, NaN, 'Color', 'c', 'LineWidth', 1, 'LineStyle', '--');
				pl(2)		= line(NaN, NaN, 'Color', 'r', 'LineWidth', 2);
				legend(pl, {'ARESELP-predicted reflections' 'ice–bed reflection'}, 'FontSize', 20, 'Location', 'southeast', 'NumColumns', 2)
			case 2
				pl			= line(NaN, NaN, 'Color', 'b', 'LineWidth', 2);
				legend(pl, 'traced reflections', 'FontSize', 20, 'Location', 'southeast')				
		end
		if (ii == 5)
			axis xy
			xlabel('Distance (km)')
		else
			axis ij
			ax(ii).XTickLabel = {};
		end
		ylabel(ylabels{ii})
    	axis([0 (1e-3 * diff(data_cat.dist_lin(ind_range))) y_range(ii, :)])
		clim([db_min db_max])
		if (ii == 5)
			text(2, (ax(ii).YLim(2) - (0.09 * diff(ax(ii).YLim))), ['(' letters(ii) ')'], 'FontSize', 20, 'FontWeight', 'bold', 'EdgeColor', 'k', 'BackgroundColor', 'w')
		else
			text(2, (ax(ii).YLim(1) + (0.09 * diff(ax(ii).YLim))), ['(' letters(ii) ')'], 'FontSize', 20, 'FontWeight', 'bold', 'EdgeColor', 'k', 'BackgroundColor', 'w')
		end
    	ax(ii).FontSize		= 20;
		ax(ii).LineWidth    = 1;
		ax(ii).XGrid		= 'on';
		ax(ii).GridColor	= [1 1 1];
		ax(ii).GridAlpha	= 0.5;
		ax(ii).GridLineWidth= 1;
		ax(ii).GridLineStyle= '--';
		ax(ii).Layer		= 'top';
		if (ii == 5)
			colorbar('FontSize', 20, 'Limits', [db_min db_max], 'YTick', [db_min db_max], 'YTickLabel', {num2str(db_min) num2str(db_max)}, 'TickLength', 0.025, 'Location', 'north', 'Position', [0.87 0.19 0.1 0.015], 'Color', 'w')
			text(365, 2.4e3, '(dB)', 'FontSize', 20, 'Color', 'w')
		end
    	box on
	end
	axes('Position', [0.89 (0.815 - (0.19 * 2) + 0.01) 0.10 0.15], 'Color', 'w')
	hold on
	image(BM5.x(1:10:end), BM5.y(1:10:end), BM5_mask_color);
	axis xy image
	line(data_cat.x(ind_range(1):25:ind_range(end)), data_cat.y(ind_range(1):25:ind_range(end)), 'Color', 'k', 'LineWidth', 2)
	line(data_cat.x(ind_range(1)), data_cat.y(ind_range(1)), 'Marker', 'o', 'MarkerFaceColor', 'g', 'MarkerSize', 10, 'Color', 'k', 'LineWidth', 1)
	line(data_cat.x(ind_range(end)), data_cat.y(ind_range(end)), 'Marker', '^', 'MarkerFaceColor', 'r', 'MarkerSize', 10, 'Color', 'k', 'LineWidth', 1)
	set(gca, 'XTick', [], 'YTick', [])
	text(5.28e5, -3.11e6, ['(' letters(ii + 1) ')'], 'FontSize', 20, 'FontWeight', 'bold', 'EdgeColor', 'k', 'BackgroundColor', 'w')
	box on

%% DATED RADARGRAM FOR FLOWCHART

	num_color_db			= 2 ^ 10;
	[db_min, db_max]		= deal(-20, 80);	
    figure('Position', [100 100 1800 600], 'Color', 'w')
    colormap(bone(num_color_db))
   	hold on
	imagesc(data_cat_core.dist_lin, elev_vec_core, amp_elev_core, [db_min db_max])
	for jj = 1:pk_cat.num_layer(ind_pk_core)
		ind_good = find(~isnan(depth_pk{ind_pk_core}(jj, :)));
		if (age_type{ind_pk_core}(jj) <= 0.5)
			line(data_cat_core.dist_lin(ind_good), (elev_surf_grd_core(ind_good) - depth_pk{ind_pk_core}(jj, ind_good) + FDM_corr{ind_pk_core}(ind_good)), 'Marker', '.', 'Color', ([227 74 51] ./ 255), 'MarkerSize', 10, 'LineStyle', 'none', 'Tag', num2str(jj))
		elseif (age_type{ind_pk_core}(jj) == 1)
			line(data_cat_core.dist_lin(ind_good), (elev_surf_grd_core(ind_good) - depth_pk{ind_pk_core}(jj, ind_good) + FDM_corr{ind_pk_core}(ind_good)), 'Marker', '.', 'Color', ([253 187 132] ./ 255), 'MarkerSize', 8, 'LineStyle', 'none')
		elseif (age_type{ind_pk_core}(jj) > 1)
			line(data_cat_core.dist_lin(ind_good), (elev_surf_grd_core(ind_good) - depth_pk{ind_pk_core}(jj, ind_good) + FDM_corr{ind_pk_core}(ind_good)), 'Marker', '.', 'Color', ([254 232 200] ./ 255), 'MarkerSize', 6, 'LineStyle', 'none')
		else
			line(data_cat_core.dist_lin(ind_good), (elev_surf_grd_core(ind_good) - depth_pk{ind_pk_core}(jj, ind_good) + FDM_corr{ind_pk_core}(ind_good)), 'Marker', '.', 'Color', 'w', 'MarkerSize', 4, 'LineStyle', 'none')
		end
	end
	line(data_cat_core.dist_lin(ones(1, 2) .* 109638), [200 elev_surf_grd_core(109638)], 'Color', 'm', 'LineWidth', 3)	
	axis xy
	axis([0 data_cat_core.dist_lin(end) -500 elev_vec_core(end)])
	clim([db_min db_max])
	set(gca, 'XTick', [], 'YTick', [], 'XDir', 'reverse')
	grid off
	box on

%% NORMALIZED RADARGRAM FOR FLOWCHART

	num_color_db			= 2 ^ 10;
	[db_min, db_max]		= deal(-200, -100);
    figure('Position', [100 100 600 600], 'Color', 'w')
    colormap(bone(num_color_db))
   	hold on
	imagesc(data_cat_norm.dist_lin, elev_vec_norm, amp_elev_norm, [db_min db_max])
	for jj = 1:pk_cat.num_layer(ind_pk_norm)
		ind_good = find(~isnan(depth_pk{ind_pk_norm}(jj, :)));
		if (age_type{ind_pk_norm}(jj) <= 0.5)
			line(data_cat_norm.dist_lin(ind_good), (elev_surf_grd_norm(ind_good) - depth_pk{ind_pk_norm}(jj, ind_good) + FDM_corr{ind_pk_norm}(ind_good)), 'Marker', '.', 'Color', ([227 74 51] ./ 255), 'MarkerSize', 10, 'LineStyle', 'none')
		elseif (age_type{ind_pk_norm}(jj) == 1)
			line(data_cat_norm.dist_lin(ind_good), (elev_surf_grd_norm(ind_good) - depth_pk{ind_pk_norm}(jj, ind_good) + FDM_corr{ind_pk_norm}(ind_good)), 'Marker', '.', 'Color', ([253 187 132] ./ 255), 'MarkerSize', 8, 'LineStyle', 'none')
		elseif (age_type{ind_pk_norm}(jj) > 1)
			line(data_cat_norm.dist_lin(ind_good), (elev_surf_grd_norm(ind_good) - depth_pk{ind_pk_norm}(jj, ind_good) + FDM_corr{ind_pk_norm}(ind_good)), 'Marker', '.', 'Color', ([254 232 200] ./ 255), 'MarkerSize', 6, 'LineStyle', 'none')
		else
			line(data_cat_norm.dist_lin(ind_good), (elev_surf_grd_norm(ind_good) - depth_pk{ind_pk_norm}(jj, ind_good) + FDM_corr{ind_pk_norm}(ind_good)), 'Marker', '.', 'Color', 'w', 'MarkerSize', 4, 'LineStyle', 'none')
		end
	end
	% for jj = 1:8
	% 	line(data_cat_norm.dist_lin, (elev_surf_grd_norm - ((jj / 10) .* pk_cat.thick{ind_pk_norm})), 'Color', 'w', 'Linewidth', 2, 'LineStyle', '--')
	% end
	axis xy
	axis([0 data_cat_norm.dist_lin(end) elev_vec_norm(1) elev_vec_norm(end)])
	clim([db_min db_max])
	set(gca, 'XTick', [], 'YTick', [], 'XDir', 'reverse')
	grid off
	box on

%% TRACED AND UNTRACED LINES

	color_set				= [253 174 97; 215 25 28; 255 255 191; 171 217 233; 44 123 182; 0 0 0] ./ 255;%{'y' 'r' 'c' 'b' 'm'};
    figure('Position', [100 100 600 1000], 'Color', 'w')
    colormap([([135 206 235] ./ 255); ([159 89 39] ./ 255); 0.9 0.9 0.9])
	subplot('Position', [0.05 0.10 0.89 0.89])
	hold on
	axis equal
	axis([x_min x_max y_min y_max])
	imagesc(BM5.x, BM5.y, BM5.mask_combo_plot)
	for ii = 1:length(x_paral)
    	line(x_paral{ii}, y_paral{ii}, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
	end
	for ii = 1:length(x_merid)
    	line(x_merid{ii}, y_merid{ii}, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
	end
	for ii = 1:num_coast
    	line((1e3 .* x_coast{ii}), (1e3 .* y_coast{ii}), 'Color', 'k', 'linewidth', 1)
	end
	for ii = setdiff(1:21, ind_campaign_ignore)
		for jj = 1:num_segment(ii)
			line(x{ii}{jj}(1:20:end), y{ii}{jj}(1:20:end), 'LineWidth', 2, 'Color', color_set(1, :))
		end
	end
	for ii = setdiff(22:num_campaign, ind_campaign_ignore)
		for jj = 1:num_segment(ii)
			line(x{ii}{jj}(1:20:end), y{ii}{jj}(1:20:end), 'LineWidth', 2, 'Color', color_set(2, :))
		end
	end
	for ii = 1:length(DTU)
		line(DTU(ii).X, DTU(ii).Y, 'LineWidth', 2, 'Color', color_set(3, :))
	end	
	for ii = setdiff(1:num_AWI, [1 7])
		line(x_AWI{ii}, y_AWI{ii}, 'Marker', '.', 'MarkerSize', 4, 'Color', color_set(4, :))
	end
	for ii = [1 7]
		line(x_AWI{ii}, y_AWI{ii}, 'Marker', '.', 'MarkerSize', 4, 'Color', color_set(5, :))
	end
	for ii = 1:num_hiawatha
		line(x_hiawatha{ii}, y_hiawatha{ii}, 'Marker', '.', 'MarkerSize', 4, 'Color', color_set(5, :))
	end	
	for ii = 1:length(M19)
		line(M19(ii).X, M19(ii).Y, 'LineWidth', 2, 'Color', color_set(6, :))
	end
	pl						= gobjects(1, size(color_set, 1));
	for ii = 1:size(color_set, 1)
		pl(ii)				= line(NaN, NaN, 'LineWidth', 3, 'Color', color_set(ii, :));
	end
   	line(x_core, y_core, 'Color', 'k', 'Marker', '^', 'MarkerSize', 10, 'MarkerFaceColor', 'm', 'LineStyle', 'none', 'LineWidth', 1)
	set(gca, 'FontSize', 20, 'XTick', [], 'YTick', [])
	clim([0 2])
	fill((1e3 .* [325 425 425 325]), (1e3 .* [-3230 -3230 -3260 -3260]), 'k')
	fill((1e3 .* [425 525 525 425]), (1e3 .* [-3230 -3230 -3260 -3260]), 'w', 'EdgeColor', 'k')
	text(300e3, -3180e3, '0', 'Color', 'k', 'FontSize', 18)
	text(520e3, -3180e3, '200 km', 'Color', 'k', 'FontSize', 18)
	text(-542e3, -3238e3, '60\circN', 'Color', 'k', 'FontSize', 18, 'Rotation', -10)
	text(-330e3, -3250e3, '50\circW', 'Color', 'k', 'FontSize', 18, 'Rotation', 82)
	text(552e3, -2750e3, '65\circN', 'Color', 'k', 'FontSize', 18, 'Rotation', 15)
	text(750e3, -2975e3, '30\circW', 'Color', 'k', 'FontSize', 18, 'Rotation', -75)
	legend(pl, {'NASA 1993–2013 (v1)' 'NASA 2014–2019 (+v2)' 'DTU 1971–1979' 'AWI 1996–2013' 'AWI 2016–2018', 'drainage basins'}, ...
		   'FontSize', 20, 'Location', 'southoutside', 'NumColumns', 2, 'Position', [0.085 0.01 0.82 0.08], 'Color', [0.95 0.95 0.95])
	grid off
	box on

%% ISOCHRONE DEPTHS
	
	[~, ind_age_iso]		= intersect(age_iso, [11.7 29 57 115]);
	num_color				= 10;
    color_fig				= [([135 206 235] ./ 255); ([159 89 39] ./ 255); 0.9 0.9 0.9; flipud(parula(num_color)); cool(num_color)];
    figure('Position', [100 100 1800 1300], 'Color', 'w')
    colormap(color_fig)	
    range_all				= [500 2500; 0 500];
	range_tmp				= NaN(1, (num_color + 1));
	for ii = 1:2
		range_tmp(ii, :)	= linspace(range_all(ii, 1), range_all(ii, 2), (num_color + 1));
	end
	ax						= gobjects(2, 4);
	for ii = 1:2
		for jj = 1:4
    		ax(ii, jj)      = subplot('Position', [(0.01 + ((jj - 1) * 0.235)) (0.52 - (0.50 * (ii - 1))) 0.23 0.45]);		
    		hold on
    		axis equal
    		axis([x_min x_max y_min y_max])
    		imagesc(BM5.x, BM5.y, BM5.mask_combo_plot, [0 2])
			switch ii
				case 1
					plot_tmp= 2 + discretize(squeeze(depth_iso_smooth(:, :, ind_age_iso(jj))), [-Inf range_tmp(ii, 2:(end - 1)) Inf]);
				case 2
					plot_tmp= 2 + num_color + discretize(squeeze(depth_iso_uncert_tot_smooth(:, :, ind_age_iso(jj))), [-Inf range_tmp(ii, 2:(end - 1)) Inf]);
			end
			imagesc(xx_grd(1, :), flipud(yy_grd(:, 1)), plot_tmp, 'AlphaData', ~isnan(plot_tmp))
    		for kk = 1:length(x_paral)
        		line(x_paral{kk}, y_paral{kk}, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
    		end
    		for kk = 1:length(x_merid)
        		line(x_merid{kk}, y_merid{kk}, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
    		end
			for kk = 1:num_coast
        		line((1e3 .* x_coast{kk}), (1e3 .* y_coast{kk}), 'Color', 'k', 'linewidth', 1)
			end
			for kk = 1:length(M19)
				line(M19(kk).X, M19(kk).Y, 'LineWidth', 2, 'Color', 'k')
			end
			line(x_core, y_core, 'Color', 'k', 'Marker', '^', 'MarkerSize', 10, 'MarkerFaceColor', 'm', 'LineStyle', 'none', 'LineWidth', 1)
    		[ax(ii, jj).FontSize, ax(ii, jj).XTick, ax(ii, jj).YTick] ...
							= deal(20, [], []);
			clim([0 size(color_fig, 1)])
    		fill((1e3 .* [325 425 425 325]), (1e3 .* [-3230 -3230 -3260 -3260]), 'k')
    		fill((1e3 .* [425 525 525 425]), (1e3 .* [-3230 -3230 -3260 -3260]), 'w', 'EdgeColor', 'k')
    		text(300e3, -3180e3, '0', 'Color', 'k', 'FontSize', 18)
    		text(520e3, -3180e3, '200 km', 'Color', 'k', 'FontSize', 18)
			if (ii == 1)
				title([num2str(age_iso(ind_age_iso(jj))) ' ka'], 'Color', 'k', 'FontSize', 24, 'FontWeight', 'bold')
			end
			text(-575e3, -800e3, ['(' letters(jj + ((ii - 1) * 4)) ')'], 'Color', 'k', 'FontSize', 24, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', 'k')
    		text(-625e3, -3205e3, '60\circN', 'Color', 'k', 'FontSize', 18, 'Rotation', -10)
    		text(-340e3, -3270e3, '50\circW', 'Color', 'k', 'FontSize', 18, 'Rotation', 82)
    		text(525e3, -2755e3, '65\circN', 'Color', 'k', 'FontSize', 18, 'Rotation', 15)
    		text(725e3, -2910e3, '30\circW', 'Color', 'k', 'FontSize', 18, 'Rotation', -75)
			if (jj == length(ind_age_iso))
				tick_str	= cell(1, (num_color + 1));
				for kk = 1:(num_color + 1)
            		tick_str{kk} ...
							= num2str(range_tmp(ii, kk));
				end
				if (ii == 1)
					tick_str{1} ...
							= ['≤' tick_str{1}];
				end
				tick_str{end} ...
							= ['≥' tick_str{end}];
				colorbar('FontSize', 20, 'Limits', ([3 (num_color + 3)] + (num_color * (ii - 1))), 'YTick', ((3:(num_color + 3)) + (num_color * (ii - 1))), 'YTickLabel', tick_str, 'TickLength', 0.03, 'Position', [0.95 (0.52 - (0.50 * (ii - 1))) 0.01 0.45])
				switch ii
					case 1
						text((x_max + 10e3), (y_max + 130e3), 'depth (m)', 'FontSize', 20)
					case 2
						text((x_max + 10e3), (y_max + 130e3), 'uncertainty (m)', 'FontSize', 20)
				end
			end
    		grid off
    		box on
		end
	end

%% DEPTH-NORMALIZED AGES

	[~, ind_depth_norm]		= intersect(depth_norm, 20:20:80);
	num_color				= 10;
    color_fig				= [([135 206 235] ./ 255); ([159 89 39] ./ 255); 0.9 0.9 0.9; flipud(parula(num_color)); cool(num_color)];
    figure('Position', [100 100 1800 1300], 'Color', 'w')
    colormap(color_fig)	
    range_all				= [0 5 0 15 0 40 10 90;
							   0 1 0 2  1 6 5 15];
	ax						= gobjects(2, 4);
	for ii = 1:2
		for jj = 1:4
    		ax(ii, jj)      = subplot('Position', [(0.01 + ((jj - 1) * 0.235)) (0.52 - (0.50 * (ii - 1))) 0.23 0.45]);
    		hold on
    		axis equal
    		axis([x_min x_max y_min y_max])
    		imagesc(BM5.x, BM5.y, BM5.mask_combo_plot, [0 2])
			range_tmp		= linspace(range_all(ii, (1 + (2 * (jj - 1)))), range_all(ii, (2 + (2 * (jj - 1)))), (num_color + 1));			
			switch ii
				case 1
					plot_tmp= 2 + discretize(squeeze(1e-3 .* age_norm_smooth(:, :, ind_depth_norm(jj))), [-Inf range_tmp(2:(end - 1)) Inf]);
				case 2
					plot_tmp= 2 + num_color + discretize(squeeze(1e-3 .* age_norm_uncert_tot_smooth(:, :, ind_depth_norm(jj))), [-Inf range_tmp(2:(end - 1)) Inf]);
			end
			imagesc(xx_grd(1, :), flipud(yy_grd(:, 1)), plot_tmp, 'AlphaData', ~isnan(plot_tmp))
    		for kk = 1:length(x_paral)
        		line(x_paral{kk}, y_paral{kk}, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
    		end
    		for kk = 1:length(x_merid)
        		line(x_merid{kk}, y_merid{kk}, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
    		end
			for kk = 1:num_coast
        		line((1e3 .* x_coast{kk}), (1e3 .* y_coast{kk}), 'Color', 'k', 'linewidth', 1)
			end
			for kk = 1:length(M19)
				line(M19(kk).X, M19(kk).Y, 'LineWidth', 2, 'Color', 'k')
			end
			switch ii
				case 1
					line(x_core, y_core, 'Color', 'k', 'Marker', '^', 'MarkerSize', 10, 'MarkerFaceColor', 'm', 'LineStyle', 'none', 'LineWidth', 1)
				case 2
					line(x_core, y_core, 'Color', 'k', 'Marker', '^', 'MarkerSize', 10, 'MarkerFaceColor', 'w', 'LineStyle', 'none', 'LineWidth', 1)
			end
    		[ax(ii, jj).FontSize, ax(ii, jj).XTick, ax(ii, jj).YTick] ...
							= deal(20, [], []);
			clim([0 size(color_fig, 1)])
    		fill((1e3 .* [325 425 425 325]), (1e3 .* [-3230 -3230 -3260 -3260]), 'k')
    		fill((1e3 .* [425 525 525 425]), (1e3 .* [-3230 -3230 -3260 -3260]), 'w', 'EdgeColor', 'k')
    		text(300e3, -3180e3, '0', 'Color', 'k', 'FontSize', 18)
    		text(520e3, -3180e3, '200 km', 'Color', 'k', 'FontSize', 18)
			if (ii == 1)
				title([num2str(depth_norm(ind_depth_norm(jj))) '% depth'], 'Color', 'k', 'FontSize', 24, 'FontWeight', 'bold')
			end
			text(-575e3, -800e3, ['(' letters(jj + ((ii - 1) * 4)) ')'], 'Color', 'k', 'FontSize', 24, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', 'k')
    		text(-625e3, -3205e3, '60\circN', 'Color', 'k', 'FontSize', 18, 'Rotation', -10)
    		text(-340e3, -3270e3, '50\circW', 'Color', 'k', 'FontSize', 18, 'Rotation', 82)
    		text(525e3, -2755e3, '65\circN', 'Color', 'k', 'FontSize', 18, 'Rotation', 15)
    		text(725e3, -2910e3, '30\circW', 'Color', 'k', 'FontSize', 18, 'Rotation', -75)
			tick_str			= cell(1, (num_color + 1));
			for kk = 1:(num_color + 1)
            	tick_str{kk}	= num2str(range_tmp(kk));
			end
			if (range_tmp(1) ~= 0)
				tick_str{1}		= ['≤' tick_str{1}];
			end
			tick_str{end}		= ['≥' tick_str{end}];
			colorbar('FontSize', 20, 'Limits', ([3 (num_color + 3)] + (num_color * (ii - 1))), 'YTick', ((3:(num_color + 3)) + (num_color * (ii - 1))), 'YTickLabel', tick_str, 'TickLength', 0.03, 'Position', [(0.23 + (0.235 * (jj - 1))) (0.52 - (0.50 * (ii - 1))) 0.01 0.45])
			if (jj == 4)
				switch ii
					case 1
						text((x_max + 10e3), (y_max + 130e3), 'age (ka)', 'FontSize', 20)
					case 2
						text((x_max + 10e3), (y_max + 130e3), 'uncertainty (ka)', 'FontSize', 20)
				end
			end
    		grid off
    		box on
		end
	end

%% DEPTH AND AGE HISTOGRAMS

	depth_norm_pk_tmp		= cell(1, pk_cat.num_file_pk);
	for ii = 1:pk_cat.num_file_pk
		depth_norm_pk_tmp{ii} ...
							= NaN(pk_cat.num_layer(ii), pk_cat.num_decim(ii));
		for jj = 1:pk_cat.num_decim(ii)
			depth_norm_pk_tmp{ii}(:, jj) ...
							= mean(depth_norm_pk{ii}(:, pk_cat.ind_decim{ii}(jj):(pk_cat.ind_decim{ii}(jj + 1) - 1)), 2);
		end
		depth_norm_pk_tmp{ii} ...
							= [depth_norm_pk_tmp{ii}(:)]';
	end
	depth_norm_pk_tmp		= 1e2 .* [depth_norm_pk_tmp{:}];
	depth_norm_pk_tmp		= depth_norm_pk_tmp(~isnan(depth_norm_pk_tmp));
	age_tmp					= cellfun(@transpose, age, 'UniformOutput', false);
	age_tmp					= 1e-3 .* [age_tmp{:}];
	age_tmp					= age_tmp(~isnan(age_tmp));

%%

	figure('Position', [10 10 600 1000], 'Color', 'w')
	subplot('Position', [0.12 0.56 0.84 0.42])
	histogram(depth_norm_pk_tmp, 0:5:100, 'Normalization', 'percentage', 'LineWidth', 1.5, 'FaceColor', 'b', 'FaceAlpha', 0.5)
	line([10 10], [0 10], 'Color', 'm', 'LineStyle', '--', 'LineWidth', 3)
	line([80 80], [0 10], 'Color', 'm', 'LineStyle', '--', 'LineWidth', 3)
	axis([0 100 0 10])
	set(gca, 'FontSize', 20, 'FontWeight', 'bold')
	xlabel('Thickness-normalized depth (%)')
	ylabel('Fraction of all traced reflections (%)')
	text(2, 9.5, '(a)', 'FontSize', 24, 'FontWeight', 'bold', 'Color', 'k')
	text(77, 9.8, 'gridded range', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'm', 'Rotation', 270)	
	box on
	subplot('Position', [0.12 0.06 0.84 0.42])
	histogram(age_tmp, 0:5:150, 'Normalization', 'percentage', 'LineWidth', 1.5, 'FaceColor', 'r', 'FaceAlpha', 0.5)
	axis([0 150 0 60])
	set(gca, 'FontSize', 20, 'FontWeight', 'bold')
	xlabel('Reflection age (ka)')
	ylabel('Fraction of all dated reflections (%)')
	text(3, 57, '(b)', 'FontSize', 24, 'FontWeight', 'bold', 'Color', 'k')
	box on
	axes('Position', [0.45 0.20 0.45 0.25], 'Color', 'none')
	histogram(age_tmp, 15:5:150, 'Normalization', 'percentage', 'LineWidth', 1.5, 'FaceColor', 'r', 'FaceAlpha', 0.5)
	axis([15 150 0 2])	
	set(gca, 'FontSize', 20, 'FontWeight', 'bold', 'XTick', [15 50 100 150])
	text(18, 1.85, '(c)', 'FontSize', 24, 'FontWeight', 'bold', 'Color', 'k')
	box on

%% V2-V1

	[~, ind_age_iso]		= intersect(age_iso, [11.7 29 57 115]);
	[~, ind_age_iso_v1]		= intersect(age_grd2_v1.age_iso, (1e3 .* [11.7 29 57 115]));
	[~, ind_depth_norm]		= intersect(depth_norm, 20:20:80);
	[~, ind_depth_norm_v1]	= intersect(round(1e2 .* age_grd2_v1.depth_norm), 20:20:80);
	num_color				= 10;
    color_fig				= [([135 206 235] ./ 255); ([159 89 39] ./ 255); 0.9 0.9 0.9; redblue(num_color, 0.5)];
    figure('Position', [100 100 1800 1300], 'Color', 'w')
    colormap(color_fig)
    range_all				= [-250 250; -5 5];
	range_tmp				= NaN(1, (num_color + 1));
	for ii = 1:2
		range_tmp(ii, :)	= linspace(range_all(ii, 1), range_all(ii, 2), (num_color + 1));
	end
	ax						= gobjects(2, 4);
	for ii = 1:2
		for jj = 1:4
    		ax(ii, jj)      = subplot('Position', [(0.01 + ((jj - 1) * 0.235)) (0.52 - (0.50 * (ii - 1))) 0.23 0.45]);
    		hold on
    		axis equal
    		axis([x_min x_max y_min y_max])
    		imagesc(BM5.x, BM5.y, BM5.mask_combo_plot, [0 2])
			switch ii
				case 1
					plot_tmp= 2 + discretize((squeeze(depth_iso_smooth(:, :, ind_age_iso(jj))) - interp2((1e3 .* age_grd2_v1.x_grd), (1e3 .* flipud(age_grd2_v1.y_grd)), squeeze(age_grd2_v1.depth_iso2(:, :, ind_age_iso_v1(jj))), xx_grd, yy_grd)), [-Inf range_tmp(ii, 2:(end - 1)) Inf]);
				case 2
					plot_tmp= 2 + discretize((1e-3 .* (squeeze(age_norm_smooth(:, :, ind_depth_norm(jj))) - interp2((1e3 .* age_grd2_v1.x_grd), (1e3 .* flipud(age_grd2_v1.y_grd)), squeeze(age_grd2_v1.age_norm2(:, :, ind_depth_norm_v1(jj))), xx_grd, yy_grd))), [-Inf range_tmp(ii, 2:(end - 1)) Inf]);
			end
			imagesc(xx_grd(1, :), flipud(yy_grd(:, 1)), plot_tmp, 'AlphaData', ~isnan(plot_tmp))
    		for kk = 1:length(x_paral)
        		line(x_paral{kk}, y_paral{kk}, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
    		end
    		for kk = 1:length(x_merid)
        		line(x_merid{kk}, y_merid{kk}, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
    		end
			for kk = 1:num_coast
        		line((1e3 .* x_coast{kk}), (1e3 .* y_coast{kk}), 'Color', 'k', 'linewidth', 1)
			end
			for kk = 1:length(M19)
				line(M19(kk).X, M19(kk).Y, 'LineWidth', 2, 'Color', 'k')
			end
			line(x_core, y_core, 'Color', 'k', 'Marker', '^', 'MarkerSize', 10, 'MarkerFaceColor', 'm', 'LineStyle', 'none', 'LineWidth', 1)
    		[ax(ii, jj).FontSize, ax(ii, jj).XTick, ax(ii, jj).YTick] ...
							= deal(20, [], []);
			clim([0 size(color_fig, 1)])
    		fill((1e3 .* [325 425 425 325]), (1e3 .* [-3230 -3230 -3260 -3260]), 'k')
    		fill((1e3 .* [425 525 525 425]), (1e3 .* [-3230 -3230 -3260 -3260]), 'w', 'EdgeColor', 'k')
    		text(300e3, -3180e3, '0', 'Color', 'k', 'FontSize', 18)
    		text(520e3, -3180e3, '200 km', 'Color', 'k', 'FontSize', 18)
			switch ii
				case 1
					title([num2str(age_iso(ind_age_iso(jj))) ' ka'], 'Color', 'k', 'FontSize', 24, 'FontWeight', 'bold')
				case 2
					title([num2str(depth_norm(ind_depth_norm(jj))) '% depth'], 'Color', 'k', 'FontSize', 24, 'FontWeight', 'bold')					
			end
			text(-575e3, -800e3, ['(' letters(jj + ((ii - 1) * 4)) ')'], 'Color', 'k', 'FontSize', 24, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', 'k')
    		text(-625e3, -3205e3, '60\circN', 'Color', 'k', 'FontSize', 18, 'Rotation', -10)
    		text(-340e3, -3270e3, '50\circW', 'Color', 'k', 'FontSize', 18, 'Rotation', 82)
    		text(525e3, -2755e3, '65\circN', 'Color', 'k', 'FontSize', 18, 'Rotation', 15)
    		text(725e3, -2910e3, '30\circW', 'Color', 'k', 'FontSize', 18, 'Rotation', -75)
			if (jj == length(ind_age_iso))
				tick_str	= cell(1, (num_color + 1));
				for kk = 1:(num_color + 1)
            		tick_str{kk} ...
							= num2str(range_tmp(ii, kk));
				end
				if (ii == 1)
					tick_str{1} ...
							= ['≤' tick_str{1}];
				end
				tick_str{end} ...
							= ['≥' tick_str{end}];
				colorbar('FontSize', 20, 'Limits', [3 (num_color + 3)], 'YTick', (3:(num_color + 3)), 'YTickLabel', tick_str, 'TickLength', 0.03, 'Position', [0.95 (0.52 - (0.50 * (ii - 1))) 0.01 0.45])
				switch ii
					case 1
						text((x_max - 170e3), (y_max + 130e3), 'depth difference (m)', 'FontSize', 20)
					case 2
						text((x_max - 130e3), (y_max + 130e3), 'age difference (ka)', 'FontSize', 20)
				end
			end
    		grid off
    		box on
		end
	end

%% LAYER TYPE

	num_color				= 5;
    color_fig				= [([135 206 235] ./ 255); ([159 89 39] ./ 255); 0.9 0.9 0.9; parula(6)];
    figure('Position', [100 100 1600 650], 'Color', 'w')
    colormap(color_fig)
	range_tmp				= [0 25];
	ax						= gobjects(1, 4);
	titles					= {'core dated' 'matched to core-dated' 'interpolated age' 'undated'};
	for ii = 1:4
		ax(ii)				= subplot('Position', [(0.01 + ((ii - 1) * 0.235)) 0.02 0.23 0.91]);
		hold on
		axis equal
		axis([x_min x_max y_min y_max])
		imagesc(BM5.x, BM5.y, BM5.mask_combo_plot, [0 2])
		for jj = 1:length(x_paral)
    		line(x_paral{jj}, y_paral{jj}, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
		end
		for jj = 1:length(x_merid)
    		line(x_merid{jj}, y_merid{jj}, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
		end
		for jj = 1:num_coast
    		line((1e3 .* x_coast{jj}), (1e3 .* y_coast{jj}), 'Color', 'k', 'linewidth', 1)
		end
		for jj = 1:length(M19)
			line(M19(jj).X, M19(jj).Y, 'LineWidth', 2, 'Color', 'k')
		end
		line(x_core, y_core, 'Color', 'm', 'Marker', '^', 'MarkerSize', 10, 'LineStyle', 'none', 'LineWidth', 2)
		for jj = setdiff(1:num_campaign, ind_campaign_ignore)
			for kk = 1:num_segment(jj)
				line(x{jj}{kk}(1:20:end), y{jj}{kk}(1:20:end), 'LineWidth', 1, 'Color', 'k')
			end
		end
		[x_cat, y_cat, num_layer_cat] ...
							= deal([]);
		for jj = 1:pk_cat.num_file_pk
			ind_good		= find(num_layer_type{jj}(ii, :) > 0);
			[x_cat, y_cat, num_layer_cat] ...
							= deal([x_cat pk_cat.x{jj}(pk_cat.ind_decim_mid{jj}(ind_good))], [y_cat pk_cat.y{jj}(pk_cat.ind_decim_mid{jj}(ind_good))], [num_layer_cat num_layer_type{jj}(ii, ind_good)]);
		end
		[~, ind_sort]		= sort(num_layer_cat);
		scatter(x_cat(ind_sort), y_cat(ind_sort), 10, color_fig(interp1(((0:5:25) + (diff(0:5:30) ./ 2)), 4:9, num_layer_cat(ind_sort), 'nearest', 'extrap'), :), 'filled')
		[ax(ii).FontSize, ax(ii).XTick, ax(ii).YTick] ...
							= deal(20, [], []);
		clim([0 size(color_fig, 1)])
		fill((1e3 .* [325 425 425 325]), (1e3 .* [-3230 -3230 -3260 -3260]), 'k')
		fill((1e3 .* [425 525 525 425]), (1e3 .* [-3230 -3230 -3260 -3260]), 'w', 'EdgeColor', 'k')
		text(300e3, -3180e3, '0', 'Color', 'k', 'FontSize', 18)
		text(520e3, -3180e3, '200 km', 'Color', 'k', 'FontSize', 18)
		title(['(' letters(ii) ') ' titles{ii}], 'Color', 'k', 'FontSize', 24, 'FontWeight', 'bold')
		text(-625e3, -3205e3, '60\circN', 'Color', 'k', 'FontSize', 18, 'Rotation', -10)
		text(-340e3, -3270e3, '50\circW', 'Color', 'k', 'FontSize', 18, 'Rotation', 82)
		text(525e3, -2755e3, '65\circN', 'Color', 'k', 'FontSize', 18, 'Rotation', 15)
		text(725e3, -2910e3, '30\circW', 'Color', 'k', 'FontSize', 18, 'Rotation', -75)
		if (ii == 4)
			colorbar('FontSize', 20, 'Limits', [3 9], 'YTick', 3:9, 'YTickLabel', {'1' '5' '10' '15' '20' '25' '>30'}, 'TickLength', 0.0175)
			text((x_max - 170e3), (y_max + 130e3), '# reflections', 'FontSize', 20)
		end
		grid off
		box on
	end

%%
end