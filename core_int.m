% CORE_INT Determine which radar transects intersect ice-core sites.
% 
%   CORE_INT determines whether each radar transect crosses near any of the
%   deep ice-core sites and records this position if so.
% 
% Joe MacGregor (UTIG)
% Last updated: 01/31/13

clear

do_core                     = true;
do_save                     = true;
plotting                    = false;
decim                       = 25;
dir_save                    = 'mat/';

load mat/xy_all name_trans name_year num_year num_trans x_gimp y_gimp
load mat/gimp_90m gimp_info
load mat/xy_gl_gimp x_gl_gimp y_gl_gimp

%% concatenate x/y positions for each campaign

if do_core
    
    disp('Comparing core sites ...')
    
    rad_threshold           = 2.5; % km
    
    name_core               = {'Camp Century' 'Dye 3' 'GISP2' 'GRIP' 'NEEM' 'NGRIP'};
    name_core_short         = {'century' 'dye3' 'gisp2' 'grip' 'neem' 'ngrip'};    
    lat_core                = [77.18 65.18 72.6 72.58  77.45  75.10];
    lon_core                = [-61.13 -43.82 -38.5 -37.64 -51.06 -42.32];
    num_core                = length(name_core);
    
%     wgs84                   = almanac('earth', 'wgs84', 'meters');
%     ps_struct               = defaultm('ups');
%     [ps_struct.geoid, ps_struct.mapparallels, ps_struct.falsenorthing, ps_struct.falseeasting, ps_struct.origin] ...
%                             = deal(wgs84, 70, 0, 0, [90 -45 0]);
    
%     [x_core, y_core]        = mfwdtran(ps_struct, lat_core, lon_core);
    [x_core_gimp, y_core_gimp] ...
                            = projfwd(gimp_info, lat_core, lon_core);
    [x_core_gimp, y_core_gimp] ...
                            = deal((1e-3 .* x_core_gimp), (1e-3 .* y_core_gimp));
    
    int_core                = cell(1, num_year);
    int_core_mat            = [];
    
    for ii = 1:num_year
        
        disp(['Campaign ' name_year{ii} ' (' num2str(ii) ' / ' num2str(num_year) ')...'])
        
        int_core{ii}        = cell(1, num_trans(ii));
        
        for jj = 1:num_trans(ii)
            
            disp([name_trans{ii}{jj} ' (' num2str(jj) ' / ' num2str(num_trans(ii)) ')...'])
            
            for kk = 1:num_core
                [tmp1, tmp2]= min(sqrt(((x_gimp{ii}{jj} - x_core_gimp(kk)) .^ 2) + ((y_gimp{ii}{jj} - y_core_gimp(kk)) .^ 2)));
                if (tmp1 < rad_threshold)
                    int_core{ii}{jj} ...
                            = [int_core{ii}{jj}; tmp1 tmp2 kk x_gimp{ii}{jj}(tmp2) y_gimp{ii}{jj}(tmp2)];
                    int_core_mat ...
                            = [int_core_mat; ii jj tmp1 tmp2 kk x_gimp{ii}{jj}(tmp2) y_gimp{ii}{jj}(tmp2)]; %#ok<AGROW>
                    disp(['...intersects within ' sprintf('%1.2f', tmp1) ' km of ' name_core{kk} '...'])
                end
            end
        end
    end
    
    save([dir_save 'core_int'], '-v7.3', 'int_core', 'int_core_mat', 'lat_core', 'lon_core', 'name_core', 'name_core_short', 'name_trans', 'num_core', 'num_trans', 'num_year', 'rad_threshold', 'x_core_gimp', 'y_core_gimp')
    disp(['Done comparing transects to core sites and saved results in ' dir_save '.'])
    
else
    load([dir_save 'core_int']) %#ok<UNRCH>
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
%% all transects and intersections
    figure('position', [503 3 784 1103])
    hold on
    colors                  = colormap(jet(num_year));
    plot((1e-3 .* x_gl_gimp), (1e-3 .* y_gl_gimp), 'k', 'linewidth', 1)
    pk_test                 = false(1, num_year);
    for ii = 1:num_year
        for jj = 1:num_trans(ii)
            if ~isempty(int_core{ii}{jj})
                plot(x_gimp{ii}{jj}(1:decim:end), y_gimp{ii}{jj}(1:decim:end), '.', 'markersize', 12, 'color', colors(ii, :))
                if ~pk_test(ii)
                    pk_test(ii) ...
                            = true;
                end
            end
        end
    end
    for ii = 1:num_core
        plot((x_core_gimp(ii) + (rad_threshold .* cos(0:0.01:(2 * pi)))), (y_core_gimp(ii) + (rad_threshold .* sin(0:0.01:(2 * pi)))), 'k--', 'linewidth', 2)
        plot(x_core_gimp(ii), y_core_gimp(ii), 's', 'color', [0.7 0.7 0.7], 'markersize', 12, 'markerfacecolor', 'k');
    end
    p                       = zeros(1, num_year);
    for ii = 1:num_year
        if pk_test(ii)
            p(ii)           = plot(-1000, -1000, 'color', colors(ii, :), 'linewidth', 4);
        end
    end
    p                       = p(pk_test);
    p(end + 1)              = plot(x_core_gimp(1), y_core_gimp(1), 's', 'color', [0.7 0.7 0.7], 'markersize', 12, 'markerfacecolor', 'k');
    p(end + 1)              = plot(x_core_gimp(1), y_core_gimp(1), 'k--', 'linewidth', 2);
    set(gca, 'fontsize', 20)
    xlabel('Polar stereographic X (km)')
    ylabel('Polar stereographic Y (km)')
    title('Core-intersecting transects')
    axis equal tight
    legend(p, [name_year(pk_test) 'core site' [num2str(rad_threshold) '-km radius']], 'location', 'southwest', 'interpreter', 'none')
    grid on
    box on
%%
end