% CORE_INT Determine which radar transects intersect ice-core sites.
% 
%   CORE_INT determines whether each radar transect crosses near any of the
%   deep ice-core sites and records this position if so.
% 
% Joe MacGregor (UTIG)
% Last updated: 02/16/15

clear

radar_type                  = 'deep';

do_core                     = true;
do_save                     = false;
plotting                    = false;

switch radar_type
    case 'accum'
        load mat/xy_all_accum name_trans name_year num_year num_trans x y
        load mat/merge_all_accum ind_fence x_pk y_pk
    case 'deep'
        load mat/xy_all name_trans name_year num_year num_trans x y
        load mat/merge_all ind_fence x_pk y_pk
end
load mat/gimp_proj gimp_proj

wgs84                       = wgs84Ellipsoid;

%% concatenate x/y positions for each campaign

if do_core
    
    disp('Comparing core sites ...')
    
    rad_threshold           = 3; % km
    
    name_core               = {'Camp Century' 'DYE-3' 'GISP2' 'GRIP' 'NEEM' 'NorthGRIP'};
    name_core_short         = {'century' 'dye3' 'gisp2' 'grip' 'neem' 'ngrip'};
    lat_core                = [77.18 65.18 72.6 72.58  77.45  75.10];
    lon_core                = [-61.13 -43.82 -38.5 -37.64 -51.06 -42.32];
    num_core                = length(name_core);
    
    [x_core, y_core]        = projfwd(gimp_proj, lat_core, lon_core);
    [x_core, y_core]        = deal((1e-3 .* x_core), (1e-3 .* y_core));
    
    int_core                = cell(1, num_year);
    
    % matrices of intersections between original data / merged picks files and ice-core sites   
    [int_core_mat, int_core_merge] ...
                            = deal([]);
    
    for ii = 1:num_year
        
        disp(['Campaign ' name_year{ii} ' (' num2str(ii) ' / ' num2str(num_year) ')...'])
        
        int_core{ii}        = cell(1, num_trans(ii));
        
        for jj = 1:num_trans(ii)
            
            disp([name_trans{ii}{jj} ' (' num2str(jj) ' / ' num2str(num_trans(ii)) ')...'])
            
            [lat_curr, lon_curr] ...
                            = projinv(gimp_proj, (1e3 .* x{ii}{jj}), (1e3 .* y{ii}{jj}));
            
            for kk = 1:num_core
                [tmp1, tmp2]= min(1e-3 .* distance([lat_core(kk(ones(length(x{ii}{jj}), 1)))' lon_core(kk(ones(length(x{ii}{jj}), 1)))'], [lat_curr' lon_curr'], wgs84Ellipsoid));
                if (tmp1 < rad_threshold)
                    int_core{ii}{jj} ...
                            = [int_core{ii}{jj}; tmp1 tmp2 kk x{ii}{jj}(tmp2) y{ii}{jj}(tmp2)];
                    int_core_mat ...
                            = [int_core_mat; ii jj tmp1 tmp2 kk x{ii}{jj}(tmp2) y{ii}{jj}(tmp2)]; %#ok<AGROW>
                    disp(['...intersects within ' sprintf('%1.2f', tmp1) ' km of ' name_core{kk} '...'])
                end
            end
            
            for kk = 1:length(ind_fence{ii}{jj})
                [lat_merge_curr, lon_merge_curr] ...
                            = projinv(gimp_proj, (1e3 .* x_pk{ii}{jj}{kk}), (1e3 .* y_pk{ii}{jj}{kk}));
                for ll = 1:num_core
                    [tmp1, tmp2] ...
                            = min(1e-3 .* distance([lat_core(ll(ones(length(x_pk{ii}{jj}{kk}), 1)))' lon_core(ll(ones(length(x_pk{ii}{jj}{kk}), 1)))'], [lat_merge_curr' lon_merge_curr'], wgs84Ellipsoid));
                    if (tmp1 < rad_threshold)
                        int_core_merge ...
                            = [int_core_merge; ii jj kk ll tmp2]; %#ok<AGROW>
                        disp(['...intersects within ' sprintf('%1.2f', tmp1) ' km of ' name_core{ll} '...'])
                    end
                end
            end
        end
    end
    
    if do_save
        switch radar_type
            case 'accum'
                save mat/core_int_accum -v7.3 int_core int_core_mat int_core_merge lat_core lon_core name_core name_core_short name_trans num_core num_trans num_year rad_threshold x_core y_core
            case 'deep'
                save mat/core_int -v7.3 int_core int_core_mat int_core_merge lat_core lon_core name_core name_core_short name_trans num_core num_trans num_year rad_threshold x_core y_core
        end
        disp('Done comparing transects to core sites and saved results in mat/.')
    end
    
else
    switch radar_type
        case 'accum'
            load mat/core_int_accum
        case 'deep'
            load mat/core_int
    end
    disp(['Loaded core intersections from ' dir_save '.'])
end
%%
% prepare core intersection results for easy loading into excel
a                           = name_year(int_core_mat(:, 1))'; % name of campaign
b                           = cell(size(int_core_mat, 1), 1);
for ii = 1:size(int_core_mat, 1)
    b{ii}                   = name_trans{int_core_mat(ii, 1)}{int_core_mat(ii, 2)}; % name of transect
end
c                           = int_core_mat(:, 3); % closest approach to core (km)
d                           = int_core_mat(:, 4); % index in original file for intersection
e                           = name_core(int_core_mat(:, 5))'; % name of core with intersection
f                           = int_core_mat(:, 6); % x value of intersection (km)
g                           = int_core_mat(:, 7); % y value of intersection (km)

%%
if plotting
%%
    set(0, 'DefaultFigureWindowStyle', 'default') %#ok<UNRCH>
    
%%
    set(0, 'DefaultFigureWindowStyle', 'docked')

%%    
    [name_trans_core, name_core_int] ...
                            = deal(cell(size(int_core_mat, 1), 1));
    for ii = 1:size(int_core_mat, 1)
        name_trans_core{ii} = name_trans{int_core_mat(ii, 1)}{int_core_mat(ii, 2)};
        name_core_int{ii}   = name_core{int_core_mat(ii, 5)};
    end
    figure('name', 'Core intersections')
    uitable('units', 'normalized', 'position', [0.1 0.1 0.8 0.8], 'fontsize', 20, ...
            'columnname', {'Name' 'Core' 'Year' 'Transect' 'Closest approach (km)' 'Index' 'X (km)' 'Y (km)'}, 'data', [name_trans_core name_core_int mat2cell(int_core_mat(:, [1:4 6:7]), ones(size(int_core_mat, 1), 1), ones(1, 6))]);

%%
end