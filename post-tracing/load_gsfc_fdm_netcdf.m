% LOAD_GSFC_FDM_NETCDF Load GSFC FDM NetCDF output from Medley et al. (2022, TC).
% 
% Joe MacGregor (NASA)
% Last updated: 17 July 2024

clear

do_load                     = true;
do_save                     = true;
plotting                    = false;
%%
if do_load
%%    
    % NetCDF path/file
    file_nc                 = '/Users/jamacgre/OneDrive - NASA/research/data/greenland/gsfc_fdm_v1_2_1_gris.nc'; % NetCDF file downloaded on 13 June 2013
    
    disp(['Loading ' file_nc '...'])
    
    % get variable names
    nc_info                 = ncinfo(file_nc);
    num_var                 = length(nc_info.Variables);
    num_var_req				= 4; % only need x, y, time, FAC which are the first four
	
    % load in variables using their original names
    [units, var_long]       = deal(cell(1, 4));
	for ii = 1:num_var_req
        var_long{ii}        = nc_info.Variables(ii).Name;
        eval([var_long{ii} ' = ncread(file_nc, var_long{ii});'])
        if (ii > 3)
            units{ii}       = nc_info.Variables(ii).Attributes(4).Value; % dimensionalized units
		else
            units{ii}       = nc_info.Variables(ii).Attributes(3).Value; % dimensionalized units
        end
	end
	x						= x'; % make row vector
	time_frac_yr			= time; % rename to make initial units more explicit
	clear time
	
	time_dt					= datetime((datenum('1-Jan-1980') + 2.5) + (5 .* (0:(length(time_frac_yr) - 1)))', 'ConvertFrom', 'datenum'); %#ok<DATNM> % regenerate using approximate start date of time dataset and stated interval (5 days)
	
	% limit to period of interest to NASA radar sounding
	time_lim				= datetime([19930101 20200101], 'ConvertFrom', 'yyyymmdd'); % years that bound surveys
	ind_time_good			= find((time_dt >= time_lim(1)) & (time_dt < time_lim(2)));
	[time_dt, FAC]			= deal(time_dt(ind_time_good), FAC(:, :, ind_time_good));
	
	% limit to region of non-NaN values
	ind_x_good				= find(sum(sum(~isnan(FAC), 3)), 1):find(sum(sum(~isnan(FAC), 3)), 1, 'last');
	ind_y_good				= find(sum(sum(~isnan(FAC), 3), 2), 1):find(sum(sum(~isnan(FAC), 3), 2), 1, 'last');
	[x, y, FAC]				= deal(x(ind_x_good), y(ind_y_good), FAC(ind_y_good, ind_x_good, :));
	
    % save output in MATLAB format
    if do_save
        disp('Saving mat/gsfc_fdm_121.mat...')
        save('/Users/jamacgre/OneDrive - NASA/research/matlab/pickgui_v2/mat/gsfc_fdm_121.mat', '-v7.3', 'FAC', 'nc_info', 'time_dt', 'units', 'x', 'y')
        disp('Saved mat/gsfc_fdm_121.mat.')
	end    
end

%%
if plotting
%%
    set(0, 'DefaultFigureWindowStyle', 'docked')
    
%%
    set(0, 'DefaultFigureWindowStyle', 'default')
    
%%
end