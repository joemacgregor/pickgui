% NCWRITE_AGE_GRD2 Make NetCDF database of gridded Greenland radiostratigraphy v2.
% 
% Joe MacGregor
% Last updated: 31 July 2024

clear
% close all

path_mat					= '/Users/jamacgre/OneDrive - NASA/research/matlab/pickgui_v2/mat/';
path_nc						= '/Users/jamacgre/OneDrive - NASA/research/matlab/pickgui_v2/misc/';
file_nc                     = 'Greenland_isochrone_grid_v2.nc';
nc_ver                      = 0.1;
val_nodata                  = NaN;

load([path_mat 'depth_iso_krige.mat'], 'age_iso', 'depth_iso', 'depth_iso_std', 'xx_grd', 'yy_grd')
[x_grd, y_grd]				= deal(xx_grd(1, :)', flipud(yy_grd(:, 1)));

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

BM5.mask_gris_grd			= interp2(BM5.x, BM5.y, BM5.mask_gris, xx_grd, yy_grd, 'nearest');

depth_iso					= flipud(depth_iso);
depth_iso_std				= flipud(depth_iso_std);

[depth_iso(~BM5.mask_gris_grd(:, :, ones(1, 1, length(age_iso)))), depth_iso_std(~BM5.mask_gris_grd(:, :, ones(1, 1, length(age_iso))))] ...
							= deal(NaN); % NaN out off-ice locations

[depth_iso_smooth, depth_iso_std_smooth] ...
							= deal(NaN(size(depth_iso)));
win_sz						= 11 .* ones(1, length(age_iso));
win_sz(7:end)				= [13 15 15 15]; % smaller windowing for older/deeper isochrones

for ii = 1:size(depth_iso, 3)
	depth_iso_smooth(:, :, ii) ...
							= smoothdata2(depth_iso(:, :, ii), 'gaussian', win_sz(ii), 'includemissing');
	depth_iso_std_smooth(:, :, ii) ...
							= smoothdata2(depth_iso_std(:, :, ii), 'gaussian', win_sz(ii), 'includemissing');
end

% for ii = 1:length(age_iso)
% 	figure
% 	colormap(jet)
% 	subplot(1, 2, 1)
% 	imagesc(x_grd, y_grd, depth_iso(:, :, ii), 'AlphaData', ~isnan(depth_iso(:, :, ii)))
% 	clim([500 2000])
% 	colorbar
% 	axis xy equal
% 	subplot(1, 2, 2)
% 	imagesc(x_grd, y_grd, depth_iso_smooth(:, :, ii), 'AlphaData', ~isnan(depth_iso_smooth(:, :, ii)))
% 	clim([500 2000])
% 	colorbar
% 	axis xy equal
% end

mapping_bm                  = ncinfo('/Users/jamacgre/OneDrive - NASA/research/data/greenland/BedMachine/BedMachineGreenland-v5.nc');
mapping_bm                  = mapping_bm.Variables(1);

% rotate grids
depth_iso_smooth			= rot90(depth_iso_smooth, -1);
depth_iso_std_smooth		= rot90(depth_iso_std_smooth, -1);

%%

% move old file if present
if exist([path_nc file_nc], 'file')
    movefile([path_nc file_nc], [path_nc file_nc(1:(end - 3)) '_old.nc'])
end

% variable writing to NetCDF file
nccreate([path_nc file_nc], 'mapping', 'Datatype', 'char', 'Format', 'netcdf4')
nccreate([path_nc file_nc], 'x', 'Dimensions', {'x' length(x_grd)}, 'Datatype', 'double')
ncwrite([path_nc file_nc], 'x', x_grd)
nccreate([path_nc file_nc], 'y', 'Dimensions', {'y' length(y_grd)}, 'Datatype', 'double')
ncwrite([path_nc file_nc], 'y', y_grd)
nccreate([path_nc file_nc], 'age', 'Dimensions', {'age' length(age_iso)}, 'Datatype', 'double')
ncwrite([path_nc file_nc], 'age', age_iso)

nc_vars                     = {'depth'	'depth_iso_smooth'		'isochrone depth, m';
                               'std'	'depth_iso_std_smooth'	'kriging uncertainty in isochrone depth, m'};
for ii = 1:size(nc_vars, 1)
    nccreate([path_nc file_nc], nc_vars{ii, 1}, 'Dimensions', {'x' length(x_grd) 'y' length(y_grd) 'age' length(age_iso)}, 'Datatype', 'double', 'FillValue', val_nodata)
    ncwrite([path_nc file_nc], nc_vars{ii, 1}, eval(nc_vars{ii, 2}))
end

% global attributes
ncwriteatt([path_nc file_nc], '/', 'Conventions', 'CF-1.7')
ncwriteatt([path_nc file_nc], '/', 'Title', 'Greenland gridded isochrone depths, v2')
ncwriteatt([path_nc file_nc], '/', 'Author', 'Joseph MacGregor')
ncwriteatt([path_nc file_nc], '/', 'Version', [string(datetime('now')) ' (v' num2str(nc_ver) ')'])
ncwriteatt([path_nc file_nc], '/', 'Data_citation', 'MacGregor et al., in prep., A revised and expanded deep radiostratigraphy of the Greenland Ice Sheet from airborne radar sounding surveys between 1993–2019')
ncwriteatt([path_nc file_nc], '/', 'nx', length(x_grd))
ncwriteatt([path_nc file_nc], '/', 'ny', length(y_grd))
ncwriteatt([path_nc file_nc], '/', 'nt', length(age_iso))
ncwriteatt([path_nc file_nc], '/', 'Projection', 'EPSG:3413 NSIDC Sea Ice Polar Stereographic North, centered on 70N, 45W')
ncwriteatt([path_nc file_nc], '/', 'proj4', '+init=epsg:3413')
ncwriteatt([path_nc file_nc], '/', 'xmin', min(x_grd))
ncwriteatt([path_nc file_nc], '/', 'ymax', max(y_grd))
ncwriteatt([path_nc file_nc], '/', 'spacing', 5e3)
ncwriteatt([path_nc file_nc], '/', 'no_data', val_nodata)

% variable-specific attributes
for ii = 2:length(mapping_bm.Attributes) % skip geoid (1)
    ncwriteatt([path_nc file_nc], 'mapping', mapping_bm.Attributes(ii).Name, mapping_bm.Attributes(ii).Value)
end
ncwriteatt([path_nc file_nc], 'mapping', 'long_name', 'mapping')
ncwriteatt([path_nc file_nc], 'mapping', 'spatial_ref', 'PROJCS["WGS 84 / NSIDC Sea Ice Polar Stereographic North",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",70],PARAMETER["central_meridian",-45],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["X",EAST],AXIS["Y",NORTH],AUTHORITY["EPSG","3413"]]')
ncwriteatt([path_nc file_nc], 'mapping', 'crs_wkt', 'PROJCS["WGS 84 / NSIDC Sea Ice Polar Stereographic North",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",70],PARAMETER["central_meridian",-45],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["X",EAST],AXIS["Y",NORTH],AUTHORITY["EPSG","3413"]]')
ncwriteatt([path_nc file_nc], 'mapping', 'proj4text', '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')
ncwriteatt([path_nc file_nc], 'mapping', 'GeoTransform', [-632e3 5e3 0 -670e3 0 -5e3])
ncwriteatt([path_nc file_nc], 'x', 'long_name', 'Cartesian x-coordinate')
ncwriteatt([path_nc file_nc], 'x', 'standard_name', 'projection_x_coordinate')
ncwriteatt([path_nc file_nc], 'x', 'units', 'meter')
ncwriteatt([path_nc file_nc], 'y', 'long_name', 'Cartesian y-coordinate')
ncwriteatt([path_nc file_nc], 'y', 'standard_name', 'projection_y_coordinate')
ncwriteatt([path_nc file_nc], 'y', 'units', 'meter')

for ii = 1:size(nc_vars, 1)
    ncwriteatt([path_nc file_nc], nc_vars{ii, 1}, 'grid_mapping', 'mapping')
    ncwriteatt([path_nc file_nc], nc_vars{ii, 1}, 'long_name', nc_vars{ii, 1})
    ncwriteatt([path_nc file_nc], nc_vars{ii, 1}, 'description', nc_vars{ii, 3})
end