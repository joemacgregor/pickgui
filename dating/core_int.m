% CORE_INT Determine which radar segments intersect ice-core sites.
% 
%   CORE_INT determines whether each segment crosses near any of the
%   deep ice-core sites and records this position if so.
% 
% Joe MacGregor (NASA)
% Last updated: 12 August 2024

clear

do_core                     = true;
do_save                     = true;

load('/Users/jamacgre/OneDrive - NASA/research/matlab/pickgui_v2/mat/xyz_all.mat', 'campaign', 'lat', 'lon', 'num_campaign', 'num_segment', 'segment', 'x', 'y')
pk_cat						= load('/Users/jamacgre/OneDrive - NASA/research/matlab/pickgui_v2/mat/pk_cat.mat', 'file_pk', 'ind_pk_merge', 'lat', 'lon', 'num_file_pk', 'num_trace');

%% concatenate x/y positions for each campaign

if do_core
    
    disp('Comparing core sites ...')
    
    rad_threshold           = 5e3; % m
    
    name_core               = {'Camp Century' 'DYE-3' 'EGRIP' 'GISP2' 'GRIP' 'NEEM' 'NorthGRIP'};
    name_core_short         = {'century' 'dye3' 'egrip' 'gisp2' 'grip' 'neem' 'ngrip'};
    lat_core                = [77.18 65.18 75.63 72.6 72.58  77.45  75.10];
    lon_core                = [-61.13 -43.82 -35.99 -38.50 -37.64 -51.06 -42.32];
    num_core                = length(name_core);
    
    [x_core, y_core]        = projfwd(projcrs(3413), lat_core, lon_core);
    
    int_core                = cell(1, num_campaign);
    
    % matrices of intersections between original data / merged picks files and ice-core sites
    int_core_mat			= deal([]);
    
    for ii = 1:num_campaign
        
        disp(['Campaign ' campaign{ii} ' (' num2str(ii) ' / ' num2str(num_campaign) ')...'])
        
        int_core{ii}        = cell(1, num_segment(ii));
        
        for jj = 1:num_segment(ii)
            
			if (isempty(lat{ii}{jj}) || (length(unique(lat{ii}{jj})) < length(lat{ii}{jj})))
				continue
			end

            disp([segment{ii}{jj} ' (' num2str(jj) ' / ' num2str(num_segment(ii)) ')...'])
            
            % core intersections based on x/y for raw data frames
			for kk = 1:num_core
                [tmp1, tmp2]= min(distance([lat_core(kk(ones(length(lat{ii}{jj}), 1)))' lon_core(kk(ones(length(lat{ii}{jj}), 1)))'], [lat{ii}{jj}' lon{ii}{jj}'], wgs84Ellipsoid));
                if (tmp1 < rad_threshold)
                    int_core{ii}{jj} ...
                            = [int_core{ii}{jj}; tmp1 tmp2 kk x{ii}{jj}(tmp2) y{ii}{jj}(tmp2)];
                    int_core_mat ...
                            = [int_core_mat; ii jj tmp1 tmp2 kk x{ii}{jj}(tmp2) y{ii}{jj}(tmp2)]; %#ok<AGROW>
                    disp(['...intersects within ' sprintf('%1.2f', (1e-3 * tmp1)) ' km of ' name_core{kk} '...'])
                end
			end
        end
    end
%%
	int_core_cat			= [];
    % core intersections for merged picks files
	for ii = 1:pk_cat.num_file_pk
        if isempty(pk_cat.lat{kk})
            continue
        end
        for jj = 1:num_core
            [tmp1, tmp2] ...
                    = min(distance([lat_core(jj(ones(pk_cat.num_trace(ii), 1)))' lon_core(jj(ones(pk_cat.num_trace(ii), 1)))'], [pk_cat.lat{ii}' pk_cat.lon{ii}'], wgs84Ellipsoid));
            if (tmp1 < rad_threshold)
                int_core_cat ...
                    = [int_core_cat; ii jj tmp2]; %#ok<AGROW>
				disp([pk_cat.file_pk{ii}(6:(end - 7)) ' intersects within ' sprintf('%1.2f', (1e-3 * tmp1)) ' km of ' name_core{jj} '...'])
            end
        end
	end
%%
    if do_save
        save('/Users/jamacgre/OneDrive - NASA/research/matlab/pickgui_v2/mat/core_int.mat', '-v7.3', ...
			'int_core', 'int_core_mat', 'int_core_cat', 'lat_core', 'lon_core', 'name_core', 'name_core_short', 'segment', 'num_core', 'num_segment', 'num_campaign', 'rad_threshold', 'x_core', 'y_core')
        disp('Done comparing transects to core sites and saved results in mat/core_int.mat.')
    end
%%    
else
	load('/Users/jamacgre/OneDrive - NASA/research/matlab/pickgui_v2/mat/core_int.mat')
    disp(['Loaded core intersections from ' dir_save '.'])
end