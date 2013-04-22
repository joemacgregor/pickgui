% RADARINT Determine intersections of radar transects.
% 
%   RADARINT determines the positions of intersections between each
%   transect and itself, all other transects in the campaign, and finally
%   all other transects in all other campaigns. It first loads the
%   positions from original data frames and concatenates them for each
%   transect.
% 
% Joe MacGregor (UTIG)
% Last updated: 04/14/13

clear

dir_save                    = 'mat/';
do_xy                       = false;
do_self                     = false;
do_int                      = false;
dist_density                = 1; % maximum intersection density in km
angle_threshold             = 10; % minimum intersection angle in degrees
decim                       = 25;
plotting                    = false;

%% concatenate x/y positions for each campaign

if do_xy
    
    name_year               = dir; %#ok<*UNRCH>
    ind_dir                 = false(1, length(name_year));
    for ii = 1:length(name_year)
        if name_year(ii).isdir
            ind_dir(ii)     = true;
        end
    end
    
    % extract names and ditch the ones we don't need
    name_year               = {name_year(ind_dir).name};
    name_year               = name_year(~strcmp(name_year, '.'));
    name_year               = name_year(~strcmp(name_year, '..'));
    name_year               = name_year(~strcmp(name_year, 'mat'));
    name_year               = name_year(~strcmp(name_year, 'AGAP'));
    name_year               = name_year(~strcmp(name_year, 'for_ldeo'));
    name_year               = name_year(~strcmp(name_year, 'Alex_Phi'));
    num_year                = length(name_year);
    
    [x, y, dist, dist_diff, num_frame, name_trans] ...
                            = deal(cell(1, num_year));
    num_trans               = zeros(1, num_year);
    
    wgs84                   = almanac('earth', 'wgs84', 'meters');
    ps_struct               = defaultm('ups');
    [ps_struct.geoid, ps_struct.mapparallels, ps_struct.falsenorthing, ps_struct.falseeasting, ps_struct.origin] ...
                            = deal(wgs84, 70, 0, 0, [90 -45 0]);
    
    disp('Extracting x/y positions from data frames...')
    
    for ii = 1:num_year
        
        disp(['Campaign ' name_year{ii} ' (' num2str(ii) ' / ' num2str(num_year) ')...'])
        
        name_trans{ii}      = dir([name_year{ii} '/data/']);
        ind_dir             = false(1, length(name_trans{ii}));
        for jj = 1:length(name_trans{ii})
            if name_trans{ii}(jj).isdir
                ind_dir(jj) = true;
            end
        end
        
        name_trans{ii}      = {name_trans{ii}(ind_dir).name};
        name_trans{ii}      = name_trans{ii}(~strcmp(name_trans{ii}, '.'));
        name_trans{ii}      = name_trans{ii}(~strcmp(name_trans{ii}, '..'));
        
        num_trans(ii)       = length(name_trans{ii});
        num_frame{ii}       = zeros(1, num_trans(ii));
        [x{ii}, y{ii}, dist{ii}, dist_diff{ii}] ...
                            = deal(cell(1, num_trans(ii)));
        
        for jj = 1:num_trans(ii)
            
            disp([ name_trans{ii}{jj} ' (' num2str(jj) ' / ' num2str(num_trans(ii)) ')...'])
            
            dir_curr        = dir([name_year{ii} '/data/' name_trans{ii}{jj} '/*.mat']);
            dir_curr        = {dir_curr.name};
            num_frame{ii}(jj) ...
                            = length(dir_curr);
            
            for kk = 1:num_frame{ii}(jj)
                
                tmp         = load([name_year{ii} '/data/' name_trans{ii}{jj} '/' dir_curr{kk}], 'Latitude', 'Longitude', 'GPS_time');
                [lat_curr, lon_curr, time_curr] ...
                            = deal(tmp.Latitude, tmp.Longitude, tmp.GPS_time);
                if any(isnan(lat_curr))
                    [lon_curr, time_curr] ...
                            = deal(lon_curr(~isnan(lat_curr)), time_curr(~isnan(lat_curr)));
                    lat_curr= lat_curr(~isnan(lat_curr));
                end
                
                % ignore overlapping data
                if (kk > 1)
                    if any(time_curr < time_old(end))
                        ind_tmp ...
                            = find(time_curr > time_old(end));
                        if ~isempty(ind_tmp)
                            [lat_curr, lon_curr] ...
                            = deal(lat_curr(ind_tmp), lon_curr(ind_tmp));
                        else
                            continue % no points in this frame that don't overlap with previous one
                        end
                    end
                end
                
                [x_curr, y_curr] ...
                            = mfwdtran(ps_struct, lat_curr, lon_curr);
                [x{ii}{jj}, y{ii}{jj}] ...
                            = deal([x{ii}{jj} (1e-3 .* single(x_curr))], [y{ii}{jj} (1e-3 .* single(y_curr))]); % m to km
                
                time_old    = time_curr;
                
            end
            
            dist{ii}{jj}    = cumsum([0 sqrt((diff(x{ii}{jj}) .^ 2) + (diff(y{ii}{jj}) .^ 2))]);
            dist_diff{ii}{jj} ...
                            = diff(dist{ii}{jj});
            
        end
    end
    
    save([dir_save 'xy_all'], 'x', 'y', 'dist', 'dist_diff', 'num_frame', 'num_trans', 'num_year', 'name_year', 'name_trans')
    disp(['Done extracting x/y positions and saved in ' dir_save '.'])
    
else
    load([dir_save 'xy_all'])
    disp(['Loaded x/y positions from ' dir_save '.'])
end

%% calculate self intersections

if do_self
    
    disp('Finding transect intersections with themselves...')
    
    dist_threshold          = 0.2; % km, distance threshold within which an intersection is considered
%     [int_self, num_self]    = deal(cell(1, num_year));
    
    load mat/int_self

    for ii = 19%:num_year
        
        disp(['Campaign ' name_year{ii} ' (' num2str(ii) ' / ' num2str(num_year) ')...'])
        
        int_self{ii}        = cell(1, num_trans(ii));
        num_self{ii}        = zeros(1, num_trans(ii));
        
        for jj = 1:num_trans(ii)
            
            disp([ name_trans{ii}{jj} ' (' num2str(jj) ' / ' num2str(num_trans(ii)) ')...'])
            
            num_curr        = length(x{ii}{jj}); % number of points in scatter plot
            decim_vec       = 1:ceil(dist{ii}{jj}(end) ./ (100 * mean(dist_diff{ii}{jj}))):num_curr; % decimate a bit
            num_decim       = length(decim_vec); % turn decim_vec into scalar, number of indices within decim_vec
            
            for kk = 1:num_curr
                
                ind_tmp     = interp1(decim_vec, 1:num_decim, kk, 'nearest', 'extrap'); % find where kk is nearest to decim_vec
                
                % transect section before kk
                if (ind_tmp > 2)
                    
                    [dist_close, ind_close]         = min(sqrt(((x{ii}{jj}(decim_vec(1:(ind_tmp - 1))) - x{ii}{jj}(kk)) .^ 2) + ...
                                                               ((y{ii}{jj}(decim_vec(1:(ind_tmp - 1))) - y{ii}{jj}(kk)) .^ 2))); % finding the closest index below ind_tmp and the closest distance
                    
                    if (dist_close < sqrt(((x{ii}{jj}(decim_vec(ind_tmp - 1)) - x{ii}{jj}(decim_vec(ind_tmp))) ^ 2) + ...
                                          ((y{ii}{jj}(decim_vec(ind_tmp - 1)) - y{ii}{jj}(decim_vec(ind_tmp))) ^ 2))) % pass test if dist less than decim_vec for x/y
                        
                        if ((ind_close < (ind_tmp - 1)) && (ind_close > 1)) % only looking at indices that have 2 or more index gap for ind_tmp and greater than first index
                            
                            [dist_best, ind_best]   = min(sqrt(((x{ii}{jj}(decim_vec(ind_close - 1):decim_vec(ind_close + 1)) - x{ii}{jj}(kk)) .^ 2) + ...
                                                               ((y{ii}{jj}(decim_vec(ind_close - 1):decim_vec(ind_close + 1)) - y{ii}{jj}(kk)) .^ 2)));
                            ind_best                = ind_best + decim_vec(ind_close - 1) - 1; % bring ind_best back to starting pt at ind_close instead of (ind_close-1)
                            
                        elseif ((ind_close == (ind_tmp - 1)) && (ind_close > 1)) % if last index in num_decim
                            
                            [dist_best, ind_best]   = min(sqrt(((x{ii}{jj}(decim_vec(ind_close - 1):decim_vec(ind_close)) - x{ii}{jj}(kk)) .^ 2) + ...
                                                               ((y{ii}{jj}(decim_vec(ind_close - 1):decim_vec(ind_close)) - y{ii}{jj}(kk)) .^ 2)));
                            ind_best                = ind_best + decim_vec(ind_close - 1) - 1;
                            
                        elseif (ind_close <= 1)
                            
                            [dist_best, ind_best]   = deal(NaN);
                            
                        end
                        
                        if (dist_best < dist_threshold) % if dist_best < dist_threshold then it is a point of intersection, MATCH!
                            int_self{ii}{jj}        = [int_self{ii}{jj}; kk ind_best dist_best]; %#ok<*SAGROW> % 1st and 2nd column are index mathes, 3rd column = dist_best
                        end
                    end
                end
                
                % transect section after kk
                if (ind_tmp < num_decim)
                    
                    [dist_close, ind_close]         = min((sqrt(((x{ii}{jj}(decim_vec((ind_tmp + 1):num_decim)) - x{ii}{jj}(kk)) .^ 2) + ((y{ii}{jj}(decim_vec((ind_tmp + 1):num_decim)) - y{ii}{jj}(kk)) .^ 2))));
                    ind_close                       = ind_close + ind_tmp; % adjust for subset, finding closest index above ind_tmp
                    
                    if (dist_close < sqrt(((x{ii}{jj}(decim_vec(ind_tmp + 1)) - x{ii}{jj}(decim_vec(ind_tmp))) ^ 2) + ((y{ii}{jj}(decim_vec(ind_tmp + 1)) - y{ii}{jj}(decim_vec(ind_tmp))) ^ 2)))
                        
                        if ((ind_close > (ind_tmp + 1)) && (ind_close < (num_decim - 1))) % indices in between first and 2nd to last
                            
                            [dist_best, ind_best]   = min(sqrt(((x{ii}{jj}(decim_vec(ind_close - 1):decim_vec(ind_close + 1)) - x{ii}{jj}(kk)) .^ 2) + ...
                                                               ((y{ii}{jj}(decim_vec(ind_close - 1):decim_vec(ind_close + 1)) - y{ii}{jj}(kk)) .^ 2)));
                            ind_best                = ind_best + decim_vec(ind_close - 1) - 1; % bring ind_best back to starting pt at ind_close instead of (ind_close-1)
                        
                        elseif ((ind_close == (ind_tmp + 1)) && (ind_close < (num_decim - 1))) % first indices possible (2) and before the 2nd to last num_decim
                            
                            [dist_best, ind_best]   = min(sqrt(((x{ii}{jj}(decim_vec(ind_close):decim_vec(ind_close + 1)) - x{ii}{jj}(kk)) .^ 2) + ...
                                                               ((y{ii}{jj}(decim_vec(ind_close):decim_vec(ind_close + 1)) - y{ii}{jj}(kk)) .^ 2)));
                            ind_best                = ind_best + decim_vec(ind_close) - 1; % bring ind_best back to starting pt at ind_close instead of (ind_close+1)
                            
                        end
                        
                        if (dist_best < dist_threshold) % if dist_best < dist_threshold then it is an index match point, MATCH!
                            int_self{ii}{jj}        = [int_self{ii}{jj}; kk ind_best dist_best]; % 1st and 2nd column are index mathes, 3rd column = dist_best
                        end
                    end
                end
            end
            
            % next do some intersection trimming
            if isempty(int_self{ii}{jj})
                
                continue
                
            else
                
                % trim to unique intersections
                for kk = 1:size(int_self{ii}{jj}, 1)
                    int_self{ii}{jj}(kk, 1:2)       = sort(int_self{ii}{jj}(kk, 1:2));
                end
                int_self{ii}{jj}                    = unique(int_self{ii}{jj}, 'rows');
                
                % remove similar intersections
                ind_diff                            = diff(int_self{ii}{jj}(:, 2)); % find the differences between the intersections for all rows in first column
                if any(ind_diff < 3) % if any rows that are alike within 3 indices, then remove all but the best point
                    
                    ind2remove                      = [];
                    kk                              = find((ind_diff < 3), 1); % finds the first alike intersections that are less than 3 when compared
                    
                    while (ind_diff(kk) < 3) % as long as ind_diff < 3, add to kk until all alike intersections are counted
                        
                        ll                          = kk;
                        while (ind_diff(ll) < 3);
                            ll                      = ll + 1; % building ll until end of alikes
                            if (ll > length(ind_diff)); % if they're are no more alike intersections then stop
                                break
                            end
                        end
                        
                        [~, ind_min]                = min(int_self{ii}{jj}(kk:ll, 3)); % find the minimum value for the alike group of intersections
                        ind_min                     = ind_min + kk - 1; % adjust for subset
                        ind2remove                  = [ind2remove; setdiff(kk:ll, ind_min)']; % creates an array of indices other than the minimum, then will later be removed with setdiff
                        kk                          = ll + 1; % jj will jump to next section of close intersections (skips 2 b/c the next index can not be similar)
                        
                        if (kk > length(ind_diff)) % break if there are no more alike intersections
                            break
                        end
                    end
                    
                    int_self{ii}{jj}                = int_self{ii}{jj}(setdiff(1:size(int_self{ii}{jj}, 1), ind2remove), :); % all but the minimum and non-alike values are removed from int_self
                    
                end
                
                num_self{ii}(jj)                    = size(int_self{ii}{jj}, 1);
                
            end
        end
    end
    
    save([dir_save 'int_self'], 'int_self', 'num_self')
    disp(['Done calculating self intersections and saved in ' dir_save '.'])
    
else
    load([dir_save 'int_self'])
    disp(['Loaded self intersections from ' dir_save '.'])
end

%% Calculate intersections between each transect all others (all others in that campaign and all other years)

if do_int
    
    disp('Finding all other intersections...')
    
    num_int                 = deal(cell(1, num_year));
    int_all                 = [];
    
    for ii = 1:num_year
        
        disp(['Campaign ' name_year{ii} ' (' num2str(ii) ' / ' num2str(num_year) ')...'])
        
        num_int{ii}         = zeros(1, num_trans(ii));
        
        for jj = ii:num_year
            
            disp(['...compared with ' name_year{jj} '...'])
            
            tmp             = size(int_all, 1);
            
            if (ii ~= jj) % transects from different years
                
                for kk = 1:num_trans(ii)
                    
                    disp(['...transect ' name_trans{ii}{kk} ' compared with...' ])
                    
                    ind_tmp1                            = 1:decim:length(x{ii}{kk});
                    ind_tmp1(end)                       = length(x{ii}{kk});
                    
                    for ll = 1:num_trans(jj)
                        
                        disp(['...' name_trans{jj}{ll} '...'])
                        
                        ind_tmp2                        = 1:decim:length(x{jj}{ll});
                        ind_tmp2(end)                   = length(x{jj}{ll});
                        
                        [~, ~, ind_int1, ind_int2]      = intersecti(x{ii}{kk}(ind_tmp1), y{ii}{kk}(ind_tmp1), x{jj}{ll}(ind_tmp2), y{jj}{ll}(ind_tmp2)); % intersection of two vectors using linear interpolation
                        
                        if ~isempty(ind_int1)
                            
                            for mm = 1:length(ind_int1)
                                
                                if (ind_int1(mm) == 1)
                                    ind_tmp3                = 1:ind_tmp1(2);
                                elseif (ind_int1(mm) == length(ind_tmp1))
                                    ind_tmp3                = ind_tmp1(end - 1):ind_tmp1(end);
                                else
                                    ind_tmp3                = ind_tmp1(ind_int1(mm) - 1):ind_tmp1(ind_int1(mm) + 1);
                                end
                                
                                if (ind_int2(mm) == 1)
                                    ind_tmp4                = 1:ind_tmp2(2);
                                elseif (ind_int2(mm) == length(ind_tmp2))
                                    ind_tmp4                = ind_tmp2(end - 1):ind_tmp2(end);
                                else
                                    ind_tmp4                = ind_tmp2(ind_int2(mm) - 1):ind_tmp2(ind_int2(mm) + 1);
                                end
                                
                                [~, ~, ind_int3, ind_int4]  = intersecti(x{ii}{kk}(ind_tmp3), y{ii}{kk}(ind_tmp3), x{jj}{ll}(ind_tmp4), y{jj}{ll}(ind_tmp4)); % intersection of two vectors using linear interpolation
                                
                                if ~isempty(ind_int3) % only when there is an intersection will we make a matrix for ind_match_trans
                                    [ind_int3, ind_int4]    = deal((ind_int3 + ind_tmp3(1) - 1), (ind_int4 + ind_tmp4(1) - 1));
                                    if ~isempty(find((diff(dist{ii}{kk}(ind_int3)) < dist_density), 1)) % ensure density of consecutive intersections is not too great
                                        tmp1                = find((diff(dist{ii}{kk}(ind_int3)) >= dist_density));
                                        [ind_int3, ind_int4]= deal(ind_int3(tmp1), ind_int4(tmp1));
                                    end
                                    if ~isempty(find((diff(dist{jj}{ll}(ind_int4)) < dist_density), 1))
                                        tmp1                = find(diff(dist{jj}{ll}(ind_int4)) >= dist_density);
                                        [ind_int3, ind_int4]= deal(ind_int3(tmp1), ind_int4(tmp1));
                                    end
                                    if ~isempty(ind_int3)
                                        int_all             = [int_all; repmat([ii kk], length(ind_int3), 1) ind_int3 NaN(length(ind_int3), 1) repmat([jj ll], length(ind_int3), 1) ind_int4 NaN(length(ind_int3), 1)];
                                    end
                                end
                            end
                        end
                    end
                end
                
            else % transects from the same year
                
                for kk = 1:(num_trans(ii) - 1) % beginning transect within year
                    
                    disp(['...transect ' name_trans{ii}{kk} ' compared with...' ])
                    
                    ind_tmp1                            = 1:decim:length(x{ii}{kk});
                    ind_tmp1(end)                       = length(x{ii}{kk});
                    
                    for ll = (kk + 1):num_trans(ii) % look at all of the following transects (no self)
                        
                        disp(['...' name_trans{jj}{ll} '...']);
                        
                        ind_tmp2                        = 1:decim:length(x{jj}{ll});
                        ind_tmp2(end)                   = length(x{jj}{ll});
                        
                        [~, ~, ind_int1, ind_int2]      = intersecti(x{ii}{kk}(ind_tmp1), y{ii}{kk}(ind_tmp1), x{jj}{ll}(ind_tmp2), y{jj}{ll}(ind_tmp2)); % intersection of two vectors using linear interpolation
                        
                        if ~isempty(ind_int1)
                            
                            for mm = 1:length(ind_int1)
                                
                                if (ind_int1(mm) == 1)
                                    ind_tmp3                = 1:ind_tmp1(2);
                                elseif (ind_int1(mm) == length(ind_tmp1))
                                    ind_tmp3                = ind_tmp1(end - 1):ind_tmp1(end);
                                else
                                    ind_tmp3                = ind_tmp1(ind_int1(mm) - 1):ind_tmp1(ind_int1(mm) + 1);
                                end
                                
                                if (ind_int2(mm) == 1)
                                    ind_tmp4                = 1:ind_tmp2(2);
                                elseif (ind_int2(mm) == length(ind_tmp2))
                                    ind_tmp4                = ind_tmp2(end - 1):ind_tmp2(end);
                                else
                                    ind_tmp4                = ind_tmp2(ind_int2(mm) - 1):ind_tmp2(ind_int2(mm) + 1);
                                end
                                
                                [~, ~, ind_int3, ind_int4]  = intersecti(x{ii}{kk}(ind_tmp3), y{ii}{kk}(ind_tmp3), x{jj}{ll}(ind_tmp4), y{jj}{ll}(ind_tmp4));
                                
                                if ~isempty(ind_int3) % only when there is an intersection will we make a matrix for ind_match_trans
                                    [ind_int3, ind_int4]    = deal((ind_int3 + ind_tmp3(1) - 1), (ind_int4 + ind_tmp4(1) - 1));
                                    if ~isempty(find((abs(diff(dist{ii}{kk}(ind_int3))) < dist_density), 1)) % ensure density of consecutive intersections is not too great
                                        tmp1                = find((abs(diff(dist{ii}{kk}(ind_int3))) >= dist_density));
                                        [ind_int3, ind_int4]= deal(ind_int3(tmp1), ind_int4(tmp1));
                                    end
                                    if ~isempty(find((abs(diff(dist{jj}{ll}(ind_int4))) < dist_density), 1))
                                        tmp1                = find(abs(diff(dist{jj}{ll}(ind_int4))) >= dist_density);
                                        [ind_int3, ind_int4]= deal(ind_int3(tmp1), ind_int4(tmp1));
                                    end
                                    if ~isempty(ind_int3)
                                        int_all             = [int_all; repmat([ii kk], length(ind_int3), 1) ind_int3 NaN(length(ind_int3), 1) repmat([jj ll], length(ind_int3), 1) ind_int4 NaN(length(ind_int3), 1)];
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            ind2remove      = [];
            
            % remove intersections due to transect jumps or small intersection angles or more than one intersection per km
            for kk = (tmp + 1):size(int_all, 1)
                
                ind_tmp1    = (int_all(kk, 3) - 1):(int_all(kk, 3) + 1);
                ind_tmp1    = ind_tmp1((ind_tmp1 > 0) & (ind_tmp1 <= length(dist_diff{ii}{int_all(kk, 2)})));
                
                if any(dist_diff{ii}{int_all(kk, 2)}(ind_tmp1) > (2 * median(dist_diff{ii}{int_all(kk, 2)})))
                    ind2remove ...
                            = [ind2remove kk];
                end
                
                ind_tmp2    = (int_all(kk, 7) - 1):(int_all(kk, 7) + 1);
                ind_tmp2    = ind_tmp2((ind_tmp2 > 0) & (ind_tmp2 <= length(dist_diff{jj}{int_all(kk, 6)})));
                
                if any(dist_diff{jj}{int_all(kk, 6)}(ind_tmp2) > (2 * median(dist_diff{jj}{int_all(kk, 6)})))
                    ind2remove ...
                            = [ind2remove kk];
                end
                
                int_all(kk, [4 8]) ...
                            = [atand(diff(y{ii}{int_all(kk, 2)}(ind_tmp1([1 end]))) ./ diff(x{ii}{int_all(kk, 2)}(ind_tmp1([1 end])))) ...
                               atand(diff(y{jj}{int_all(kk, 6)}(ind_tmp2([1 end]))) ./ diff(x{jj}{int_all(kk, 6)}(ind_tmp2([1 end]))))];
                
                % transect intersection angle must be greater than angle_threshold
                if (abs(diff(int_all(kk, [4 8]))) < angle_threshold)
                    ind2remove ...
                            = [ind2remove kk];
                end
            end
            
            ind2remove      = unique(ind2remove);
            
            if ~isempty(ind2remove)
                int_all     = int_all(setdiff(1:size(int_all, 1), ind2remove), :); % the transect intersection where data actually exist
            end
            
            num_int{ii}(jj) = length((tmp + 1):size(int_all, 1));
            
        end
    end
    
    save([dir_save 'int_all'], 'int_all')
    disp(['Done calculating all intersections and saved in ' dir_save '.'])
    
else
    load([dir_save 'int_all'])
    disp(['Loaded all intersections from ' dir_save '.'])
end

%%
% populate fence log
fence                       = cell(size(int_all, 1), 4);
for ii = 1:size(int_all, 1)
    fence{ii, 1}            = name_year{int_all(ii, 1)};
    fence{ii, 2}            = name_trans{int_all(ii, 1)}{int_all(ii, 2)};
    fence{ii, 3}            = name_year{int_all(ii, 5)};
    fence{ii, 4}            = name_trans{int_all(ii, 5)}{int_all(ii, 6)};
end

[~, ind_unique]             = unique(int_all(:, [1 2 5 6]), 'rows');
fence                       = fence(ind_unique, :);

%%
if plotting
%% greenland outline
    h                       = figure(1);
    worldmap('greenland'); % displays empty map
    gl_outline              = shaperead('landareas', 'UseGeoCoords', true, 'Selector', {@(name) strcmp(name, 'Greenland'), 'Name'}); % create Greenland outline
    wgs84                   = almanac('earth', 'wgs84', 'meters'); % parameters for Earth
    ps_struct               = defaultm('ups'); %  initializes a map projection structure
    [ps_struct.geoid, ps_struct.mapparallels, ps_struct.falsenorthing, ps_struct.falseeasting, ps_struct.origin] ...
                            = deal(wgs84, 70, 0, 0, [90 -45 0]); % assigns ps_struct to variables
    [x_gl, y_gl]            = mfwdtran(ps_struct, gl_outline.Lat, gl_outline.Lon); % convert lat/lon to polar stereographic x/y
    [x_gl, y_gl]            = deal((1e-3 .* x_gl), (1e-3 .* y_gl)); % m to km
    close(h)
    
%% all transects and intersections
    figure
    hold on
    colors                  = colormap(jet(num_year));
    plot(x_gl, y_gl, 'k', 'linewidth', 1)
    for ii = 1%:num_year
        for jj = 1:num_trans(ii)
            plot(x{ii}{jj}(1:decim:end), y{ii}{jj}(1:decim:end), '.', 'markersize', 12, 'color', colors(ii, :))
%             for kk = 1:num_self{ii}(jj)
%                 plot(x{ii}{jj}(int_self{ii}{jj}(kk, 1)), y{ii}{jj}(int_self{ii}{jj}(kk, 1)), 'm.', 'markersize', 20)
%                 plot(x{ii}{jj}(int_self{ii}{jj}(kk, 2)), y{ii}{jj}(int_self{ii}{jj}(kk, 2)), 'r.', 'markersize', 20)
%             end
            plot(x{ii}{jj}(int_all(((int_all(:, 1) == ii) & (int_all(:, 2) == jj)), 3)), y{ii}{jj}(int_all(((int_all(:, 1) == ii) & (int_all(:, 2) == jj)), 3)), 'ks', 'markersize', 10, 'markerfacecolor', 'm')
            plot(x{ii}{jj}(int_all(((int_all(:, 4) == ii) & (int_all(:, 5) == jj)), 6)), y{ii}{jj}(int_all(((int_all(:, 4) == ii) & (int_all(:, 5) == jj)), 6)), 'k^', 'markersize', 10, 'markerfacecolor', 'm')
        end
    end
    set(gca, 'fontsize', 20)
    xlabel('Polar stereographic X (km)')
    ylabel('Polar stereographic Y (km)')
    axis equal tight
    grid on
    box on
    
%%
end;