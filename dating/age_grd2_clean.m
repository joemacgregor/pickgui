% AGE_GRD2_CLEAN Clean up Python-generated isochrone and normalized ages.
% 
% Joe MacGregor
% Last updated: 18 November 2024

clear

do_save						= false;
plotting					= false;

dir_mat						= '/Users/jamacgre/OneDrive - NASA/research/matlab/pickgui_v2/mat/';

[x_min, x_max, y_min, y_max]= deal(-632e3, 846e3, -3344e3, -670e3);

%%

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

%%

load([dir_mat 'depth_iso_krige.mat'], 'age_iso', 'depth_iso', 'depth_iso_uncert_tot', 'idx_inside_trim_iso_nohull', 'xx_grd', 'yy_grd')

num_age_iso					= length(age_iso);

% mask of convex hull
depth_iso_grd				= cell(1, num_age_iso);
for ii = 1:num_age_iso
	depth_iso_grd{ii}		= zeros(size(depth_iso(:, :, ii)));
	depth_iso_grd{ii}(~isnan(depth_iso(:, :, ii))) ...
							= 1;
	tmp						= flipud(rot90(reshape(idx_inside_trim_iso_nohull(:, ii), size(depth_iso, 2), size(depth_iso, 1))));
	depth_iso_grd{ii}(tmp)	= 2;
end

% re-3d
depth_iso_grd_merge			= NaN(size(depth_iso));
for ii = 1:num_age_iso
	depth_iso_grd_merge(:, :, ii) ...
							= depth_iso_grd{ii};
end

[depth_iso, depth_iso_grd_merge, depth_iso_uncert_tot] ...
							= deal(flipud(depth_iso), flipud(depth_iso_grd_merge), flipud(depth_iso_uncert_tot));

% one kilometer grid to rule them all
[x_grd, y_grd]				= deal(xx_grd(1, :)', flipud(yy_grd(:, 1)));
clear depth_iso_std depth_iso_uncert

BM5.mask_gris_grd			= interp2(BM5.x, BM5.y, BM5.mask_gris, xx_grd, yy_grd, 'nearest');

% windows for smoothing
win_sz_iso					= 11 .* ones(1, num_age_iso);
win_sz_iso(7:end)			= [13 15 15 15]; % larger windowing for older/deeper isochrones

% NaN out off-ice locations
[depth_iso(~BM5.mask_gris_grd(:, :, ones(1, 1, length(age_iso)))), depth_iso_uncert_tot(~BM5.mask_gris_grd(:, :, ones(1, 1, length(age_iso))))] ...
							= deal(NaN); 

% smooth
[depth_iso_smooth, depth_iso_uncert_tot_smooth] ...
							= deal(NaN(size(depth_iso)));
for ii = 1:num_age_iso
	depth_iso_smooth(:, :, ii) ...
							= smoothdata2(depth_iso(:, :, ii), 'gaussian', win_sz_iso(ii), 'includemissing');
	depth_iso_uncert_tot_smooth(:, :, ii) ...
							= smoothdata2(depth_iso_uncert_tot(:, :, ii), 'gaussian', win_sz_iso(ii), 'includemissing');
end

% remove original convex hull
[depth_iso_smooth(depth_iso_grd_merge < 2), depth_iso_uncert_tot_smooth(depth_iso_grd_merge < 2)] ...
							= deal(NaN);

% fill in interior gaps
for ii = 1:num_age_iso
	tmp1					= depth_iso_smooth(:, :, ii);
	tmp1(isnan(tmp1))		= 0;
	mask_hole				= logical(smallholefill((1 .* logical(tmp1)), 500, 1, '==') - logical(tmp1));
	if ~isempty(find(mask_hole, 1))
		[ind_i, ind_j]		= find(mask_hole);
		tmp2				= depth_iso_smooth(:, :, ii);
		tmp2(isnan(tmp2) & ~mask_hole) ...
							= Inf;
		tmp2				= fillmissing2(tmp2, 'natural');
		tmp2(isinf(tmp2))	= NaN;
		depth_iso_smooth(:, :, ii) ...
							= tmp2;
		tmp3				= depth_iso_uncert_tot_smooth(:, :, ii);
		tmp3(isnan(tmp3) & ~mask_hole) ...
							= Inf;
		tmp3				= fillmissing2(tmp3, 'natural');
		tmp3(isinf(tmp3))	= NaN;
		depth_iso_uncert_tot_smooth(:, :, ii) ...
							= tmp3;
	end
end

disp('isochrone overturns...')

ind_bad_iso					= [];

for ii = 1:(num_age_iso - 1)
	for jj = (ii + 1):num_age_iso
		if ~isempty(find((diff(depth_iso_smooth(:, :, [ii jj]), 1, 3) < 0), 1))
			disp([ii jj (1e2 * (length(find((diff(depth_iso_smooth(:, :, [ii jj]), 1, 3) < 0))) / length(find(BM5.mask_gris_grd))))])
			[i_tmp, j_tmp]	= ind2sub(size(xx_grd), find((diff(depth_iso_smooth(:, :, [ii jj]), 1, 3) < 0)));
			ind_bad_iso		= [ind_bad_iso; ii(ones(length(i_tmp), 1)) jj(ones(length(j_tmp), 1)) i_tmp j_tmp]; %#ok<AGROW>
		end
	end
end

% try to fix overturns using filling
for ii = 1:num_age_iso
	ind_bad_curr			= ind_bad_iso((ind_bad_iso(:, 1) == ii), [3 4]);
	if ~isempty(ind_bad_curr)
		mask_bad			= false(size(depth_iso_smooth, 1), size(depth_iso_smooth, 2));
		for jj = 1:size(ind_bad_curr, 1)
			mask_bad(ind_bad_curr(jj, 1), ind_bad_curr(jj, 2)) ...
							= true;
		end
		tmp2				= depth_iso_smooth(:, :, ii);
		tmp2(isnan(tmp2) & ~mask_bad) ...
							= Inf;
		tmp2(mask_bad)		= NaN;
		tmp2				= fillmissing2(tmp2, 'natural');
		tmp2(isinf(tmp2))	= NaN;
		depth_iso_smooth(:, :, ii) ...
							= tmp2;
		tmp3				= depth_iso_uncert_tot_smooth(:, :, ii);
		tmp3(isnan(tmp3) & ~mask_bad) ...
							= Inf;
		tmp2(mask_bad)		= NaN;
		tmp3				= fillmissing2(tmp3, 'natural');
		tmp3(isinf(tmp3))	= NaN;
		depth_iso_uncert_tot_smooth(:, :, ii) ...
							= tmp3;
	end
end

disp('isochrone overturns post...')

ind_bad_iso					= [];

for ii = 1:(num_age_iso - 1)
	for jj = (ii + 1):num_age_iso
		if ~isempty(find((diff(depth_iso_smooth(:, :, [ii jj]), 1, 3) < 0), 1))
			disp([ii jj (1e2 * (length(find((diff(depth_iso_smooth(:, :, [ii jj]), 1, 3) < 0))) / length(find(BM5.mask_gris_grd))))])
			[i_tmp, j_tmp]	= ind2sub(size(xx_grd), find((diff(depth_iso_smooth(:, :, [ii jj]), 1, 3) < 0)));
			ind_bad_iso		= [ind_bad_iso; ii(ones(length(i_tmp), 1)) jj(ones(length(j_tmp), 1)) i_tmp j_tmp]; %#ok<AGROW>
		end
	end
end

% NaN out grid cell of deepest isochrone with issue
[~, ind_unique]				= unique(ind_bad_iso(:, [3 4]), 'rows', 'stable');
ind_bad_iso					= ind_bad_iso(ind_unique, :);
for ii = 1:size(ind_bad_iso, 1)
	[depth_iso_smooth(ind_bad_iso(ii, 3), ind_bad_iso(ii, 4), ind_bad_iso(ii, 2)), depth_iso_uncert_tot_smooth(ind_bad_iso(ii, 3), ind_bad_iso(ii, 4), ind_bad_iso(ii, 2))] ...
							= deal(NaN);
end

%%

load([dir_mat 'age_norm_krige.mat'], 'age_norm', 'age_norm_uncert_tot', 'depth_norm', 'idx_inside_trim_norm_nohull')

ind_depth_norm_good			= find((depth_norm >= 10) & (depth_norm <= 80));
[age_norm, age_norm_uncert_tot, depth_norm, idx_inside_trim_norm_nohull] ...
							= deal(age_norm(:, :, ind_depth_norm_good), age_norm_uncert_tot(:, :, ind_depth_norm_good), depth_norm(ind_depth_norm_good), idx_inside_trim_norm_nohull(:, ind_depth_norm_good));
num_depth_norm				= length(depth_norm);

% mask of convex hull
age_norm_grd				= cell(1, num_depth_norm);
for ii = 1:num_depth_norm
	age_norm_grd{ii}		= zeros(size(age_norm(:, :, ii)));
	age_norm_grd{ii}(~isnan(age_norm(:, :, ii))) ...
							= 1;
	tmp						= flipud(rot90(reshape(idx_inside_trim_norm_nohull(:, ii), size(age_norm, 2), size(age_norm, 1))));
	age_norm_grd{ii}(tmp)	= 2;
end

% re-3d
age_norm_grd_merge			= NaN(size(age_norm));
for ii = 1:num_depth_norm
	age_norm_grd_merge(:, :, ii) ...
							= age_norm_grd{ii};
end

[age_norm, age_norm_grd_merge, age_norm_uncert_tot] ...
							= deal(flipud(age_norm), flipud(age_norm_grd_merge), flipud(age_norm_uncert_tot));

win_sz_norm					= 11 .* ones(1, num_depth_norm);

[age_norm(~BM5.mask_gris_grd(:, :, ones(1, 1, length(depth_norm)))), age_norm_uncert_tot(~BM5.mask_gris_grd(:, :, ones(1, 1, length(depth_norm))))] ...
							= deal(NaN); % NaN out off-ice locations

[age_norm_smooth, age_norm_uncert_tot_smooth] ...
							= deal(NaN(size(age_norm)));
for ii = 1:num_depth_norm
	age_norm_smooth(:, :, ii) ...
							= smoothdata2(age_norm(:, :, ii), 'gaussian', win_sz_norm(ii), 'includemissing');
	age_norm_uncert_tot_smooth(:, :, ii) ...
							= smoothdata2(age_norm_uncert_tot(:, :, ii), 'gaussian', win_sz_norm(ii), 'includemissing');
end

% remove original convex hull
[age_norm_smooth(age_norm_grd_merge < 2), age_norm_uncert_tot_smooth(age_norm_grd_merge < 2)] ...
							= deal(NaN);

% fill in interior gaps
for ii = 1:num_depth_norm
	tmp1					= age_norm_smooth(:, :, ii);
	tmp1(isnan(tmp1))		= 0;
	mask_hole				= logical(smallholefill((1 .* logical(tmp1)), 500, 1, '==') - logical(tmp1));
	if ~isempty(find(mask_hole, 1))
		[ind_i, ind_j]		= find(mask_hole);
		tmp2				= age_norm_smooth(:, :, ii);
		tmp2(isnan(tmp2) & ~mask_hole) ...
							= Inf;
		tmp2				= fillmissing2(tmp2, 'natural');
		tmp2(isinf(tmp2))	= NaN;
		age_norm_smooth(:, :, ii) ...
							= tmp2;
		tmp3				= age_norm_uncert_tot_smooth(:, :, ii);
		tmp3(isnan(tmp3) & ~mask_hole) ...
							= Inf;
		tmp3				= fillmissing2(tmp3, 'natural');
		tmp3(isinf(tmp3))	= NaN;
		age_norm_uncert_tot_smooth(:, :, ii) ...
							= tmp3;
	end
end

disp('age overturns...')

ind_bad_norm				= [];

for ii = 1:(num_depth_norm - 1)
	for jj = (ii + 1):num_depth_norm
		if ~isempty(find((diff(age_norm_smooth(:, :, [ii jj]), 1, 3) < 0), 1))
			disp([ii jj (1e2 * (length(find((diff(age_norm_smooth(:, :, [ii jj]), 1, 3) < 0))) / length(find(BM5.mask_gris_grd))))])
			[i_tmp, j_tmp]	= ind2sub(size(xx_grd), find((diff(age_norm_smooth(:, :, [ii jj]), 1, 3) < 0)));
			ind_bad_norm	= [ind_bad_norm; ii(ones(length(i_tmp), 1)) jj(ones(length(j_tmp), 1)) i_tmp j_tmp]; %#ok<AGROW>
		end
	end
end

% try to fix overturns using filling
for ii = 1:num_depth_norm
	ind_bad_curr			= ind_bad_norm((ind_bad_norm(:, 1) == ii), [3 4]);
	if ~isempty(ind_bad_curr)
		mask_bad			= false(size(age_norm_smooth, 1), size(age_norm_smooth, 2));
		for jj = 1:size(ind_bad_curr, 1)
			mask_bad(ind_bad_curr(jj, 1), ind_bad_curr(jj, 2)) ...
							= true;
		end
		tmp2				= age_norm_smooth(:, :, ii);
		tmp2(isnan(tmp2) & ~mask_bad) ...
							= Inf;
		tmp2(mask_bad)		= NaN;
		tmp2				= fillmissing2(tmp2, 'natural');
		tmp2(isinf(tmp2))	= NaN;
		age_norm_smooth(:, :, ii) ...
							= tmp2;
		tmp3				= age_norm_uncert_tot_smooth(:, :, ii);
		tmp3(isnan(tmp3) & ~mask_bad) ...
							= Inf;
		tmp2(mask_bad)		= NaN;
		tmp3				= fillmissing2(tmp3, 'natural');
		tmp3(isinf(tmp3))	= NaN;
		age_norm_uncert_tot_smooth(:, :, ii) ...
							= tmp3;
	end
end

disp('age volume overturns post...')

ind_bad_norm				= [];

for ii = 1:(num_depth_norm - 1)
	for jj = (ii + 1):num_depth_norm
		if ~isempty(find((diff(age_norm_smooth(:, :, [ii jj]), 1, 3) < 0), 1))
			disp([ii jj (1e2 * (length(find((diff(age_norm_smooth(:, :, [ii jj]), 1, 3) < 0))) / length(find(BM5.mask_gris_grd))))])
			[i_tmp, j_tmp]	= ind2sub(size(xx_grd), find((diff(age_norm_smooth(:, :, [ii jj]), 1, 3) < 0)));
			ind_bad_norm	= [ind_bad_norm; ii(ones(length(i_tmp), 1)) jj(ones(length(j_tmp), 1)) i_tmp j_tmp]; %#ok<AGROW>
		end
	end
end

% NaN out grid cell of deepest age with issue
[~, ind_unique]				= unique(ind_bad_norm(:, [3 4]), 'rows', 'stable');
ind_bad_norm					= ind_bad_norm(ind_unique, :);
for ii = 1:size(ind_bad_norm, 1)
	[age_norm_smooth(ind_bad_norm(ii, 3), ind_bad_norm(ii, 4), ind_bad_norm(ii, 2)), age_norm_uncert_tot_smooth(ind_bad_norm(ii, 3), ind_bad_norm(ii, 4), ind_bad_norm(ii, 2))] ...
							= deal(NaN);
end


%%

%%

if do_save
	save([dir_mat 'age_grd2_clean.mat'], 'age_iso', 'age_norm_smooth', 'age_norm_uncert_tot_smooth', 'depth_iso_smooth', 'depth_iso_uncert_tot_smooth', 'depth_norm', 'num_age_iso', 'num_depth_norm', 'x_grd', 'y_grd')
end

%%
if plotting
%%
    set(0, 'DefaultFigureWindowStyle', 'default') %#ok<UNRCH>

%%
    set(0, 'DefaultFigureWindowStyle', 'docked')

%% 

	figure
	colormap(jet)
	for ii = 1:num_age_iso
		subplot(2, 5, ii)
		imagesc(x_grd, y_grd, depth_iso_smooth(:, :, ii), 'AlphaData', ~isnan(depth_iso_smooth(:, :, ii)))
		% line(x_grd(ind_bad_iso(:, 4)), y_grd(ind_bad_iso(:, 3)), 'Color', 'm', 'Marker', '.', 'MarkerSize', 6, 'LineStyle', 'none')
		clim([500 2000])
		colorbar
		axis ij equal
		axis([x_min x_max y_min y_max])
		title([num2str(age_iso(ii)) ' ka'])
	end

%%

	figure
	colormap(jet)
	for ii = 1:num_depth_norm
		subplot(2, 4, ii)
		imagesc(x_grd, y_grd, age_norm_smooth(:, :, ii), 'AlphaData', ~isnan(age_norm_smooth(:, :, ii)))
		% line(x_grd(ind_bad_iso(:, 4)), y_grd(ind_bad_iso(:, 3)), 'Color', 'm', 'Marker', '.', 'MarkerSize', 6, 'LineStyle', 'none')
		clim([5e3 115e3])
		colorbar
		axis ij equal
		axis([x_min x_max y_min y_max])
		title([num2str(depth_norm(ii)) '%'])
	end

%%	
	figure
	colormap(jet)
	for ii = 1:num_age_iso
		subplot(2, 5, ii)
		imagesc(x_grd, y_grd, (1e2 .* (depth_iso_uncert_tot_smooth(:, :, ii) ./ depth_iso_smooth(:, :, ii))), 'AlphaData', ~isnan(depth_iso_uncert_tot_smooth(:, :, ii)))
		clim([0 20])
		colorbar
		axis ij equal
	end

end