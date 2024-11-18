% NCWRITE_AGE_GRD2 Make NetCDF database of gridded Greenland radiostratigraphy v2.
% 
% Joe MacGregor
% Last updated: 13 November 2024

clear
% close all

path_mat					= '/Users/jamacgre/OneDrive - NASA/research/matlab/pickgui_v2/mat/';
path_nc						= '/Users/jamacgre/OneDrive - NASA/research/matlab/pickgui_v2/misc/';
file_nc                     = 'Greenland_isochrone_grid_v2.nc';
nc_ver                      = 0.5;
val_nodata                  = NaN;

load([path_mat 'age_grd2_clean.mat'], 'age_iso', 'age_norm_smooth', 'age_norm_uncert_tot_smooth', 'depth_iso_smooth', 'depth_iso_uncert_tot_smooth', 'depth_norm', 'num_age_iso', 'num_depth_norm', 'x_grd', 'y_grd')

% rotate grids
depth_iso					= rot90(depth_iso_smooth, -1);
depth_iso_uncert			= rot90(depth_iso_uncert_tot_smooth, -1);
age_norm					= rot90((1e-3 .* age_norm_smooth), -1);
age_norm_uncert				= rot90((1e-3 .* age_norm_uncert_tot_smooth), -1);

%%

mapping_bm                  = ncinfo('/Users/jamacgre/OneDrive - NASA/research/data/greenland/BedMachine/BedMachineGreenland-v5.nc');
mapping_bm                  = mapping_bm.Variables(1);

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
nccreate([path_nc file_nc], 'age_iso', 'Dimensions', {'age' num_age_iso}, 'Datatype', 'double')
ncwrite([path_nc file_nc], 'age_iso', age_iso)
nccreate([path_nc file_nc], 'depth_norm', 'Dimensions', {'depth' num_depth_norm}, 'Datatype', 'double')
ncwrite([path_nc file_nc], 'depth_norm', depth_norm)

nc_vars                     = {'depth'			'depth_iso'			'isochrone depth, m';
                               'depthstd'		'depth_iso_uncert'	'total uncertainty in isochrone depth, m'
							   'age'			'age_norm'			'age at normalized depth, ka';
                               'agestd'			'age_norm_uncert'	'total uncertainty in age at normalized depth, ka'};
for ii = 1:size(nc_vars, 1)
	if any(contains(nc_vars{ii, 2}, 'depth'))
		nccreate([path_nc file_nc], nc_vars{ii, 1}, 'Dimensions', {'x' length(x_grd) 'y' length(y_grd) 'age' num_age_iso}, 'Datatype', 'double', 'FillValue', val_nodata)
	else
		nccreate([path_nc file_nc], nc_vars{ii, 1}, 'Dimensions', {'x' length(x_grd) 'y' length(y_grd) 'depth' num_depth_norm}, 'Datatype', 'double', 'FillValue', val_nodata)
	end
    ncwrite([path_nc file_nc], nc_vars{ii, 1}, eval(nc_vars{ii, 2}))
end

% global attributes
ncwriteatt([path_nc file_nc], '/', 'Conventions', 'CF-1.7')
ncwriteatt([path_nc file_nc], '/', 'Title', 'Greenland gridded isochrone depths and ages, v2')
ncwriteatt([path_nc file_nc], '/', 'Author', 'Joseph MacGregor')
ncwriteatt([path_nc file_nc], '/', 'Version', [string(datetime('now')) ' (v' num2str(nc_ver) ')'])
ncwriteatt([path_nc file_nc], '/', 'Data_citation', 'MacGregor et al., in prep., A revised and expanded deep radiostratigraphy of the Greenland Ice Sheet from airborne radar sounding surveys between 1993–2019')
ncwriteatt([path_nc file_nc], '/', 'nx', length(x_grd))
ncwriteatt([path_nc file_nc], '/', 'ny', length(y_grd))
ncwriteatt([path_nc file_nc], '/', 'na', length(age_iso))
ncwriteatt([path_nc file_nc], '/', 'nd', length(depth_norm))
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