% DATE_LAYERS Date layers using ice-core depth/age scales, lists of
% matching layers and 1D/2D interpolation/extrapolation or quasi-Nye dating
% of overlapping dated layers.
% 
% Joe MacGregor (NASA)
% Last updated: 9 August 2024

clear

dir_mat						= '/Users/jamacgre/OneDrive - NASA/research/matlab/pickgui_v2/mat/';
do_snr						= false;
do_date						= true;
do_age_check				= true;
do_match					= [true true];
do_interp					= true;
interp_type					= 'quasi Nye';
do_save						= true;
do_grd1						= true;
do_grd2						= false;
do_nye_norm					= false;

% variables of interest
depth_uncert				= 5; % depth uncertainty representing unknown firn correction, m
age_uncert_rel_max			= 0.25; % maximum relative age uncertainty
age_uncert_frac_max			= 50e-2; % maximum fraction of age uncertainty at which point to declare an age overturn
dist_int_max				= 250; % m, +/- range to extract core-intersecting layer depths
thick_diff_max				= 0.2; % fraction of ice thickness within which overlapping layer must lie
layer_diff_max				= 1; % fraction of layer thickness within which overlapping layer must lie
num_depth_norm				= 0;%25; % number of normalized depths
age_iso						= 1e3 .* [3 8 9 11.7 12.8 14.7 19 29 57 115]'; % isochrone ages, a
num_age_iso					= length(age_iso); % number of isochrones to calculate
thick_diff_max_iso			= 0.2; % fraction of ice thickness within which overlapping layer must lie (for age_norm1)
layer_diff_max_iso			= 1; % fraction of layer thickness within which overlapping layer must lie (for age_norm1 and depth_iso1)
age_diff_max_iso			= 5e3; % age range within which layer must exist (for depth_iso1)
snr_ref						= 10 ^ (5 / 10); % linear ratio (5 is dB), reference SNR value if unknown
depth_shallow_min			= 0.03; % minimum fraction of ice thickness to attempt shallow fitting
depth_shallow_max			= 0.2; % maximum fraction of ice thickness to attempt shallow fitting
num_date_loop_max			= 10; % maximum number of dating loops
age_max						= 1.5e5; % maximum permitted reflector age, a
campaign_ord				= [22 18 20 28 21 29 30 6 7 8 4 9 5 4 14 16 17 12 11 15 19 3 1 2 23 25]; % order in which to analyze campaigns based on their overall data quality
num_gauss					= [ones(1, 15) (2 .* ones(1, 8)) 1]; % number of Gaussians to use to fit Nye differences of normalized ages

if strcmp(interp_type, 'quasi Nye')
    qn_tol                  = 1e-2; % residual tolerance, m
    qn_iter_max             = 100; % maximum number of iterations in search of best-fit strain rate
else
    [qn_tol, qn_iter_max]   = deal([]); %#ok<*UNRCH>
end

% load concatenated data and core intersection data
load([dir_mat 'xyz_all.mat'], 'campaign')
load([dir_mat 'core_int.mat'], 'int_core_cat', 'num_core', 'name_core_short')
int_core_cat				= sortrows(int_core_cat, 1:3, {'ascend', 'descend', 'ascend'}); % try to get mostly newer cores to date first
load([dir_mat 'layer_bin_clean.mat'], 'id_match', 'layer_bin')
load([dir_mat 'pk_cat.mat'], 'depth', 'dist', 'file_pk', 'ind_campaign', 'ind_decim', 'ind_decim_mid', 'ind_trace_layer', 'num_decim', 'num_layer', 'num_file_pk', 'num_trace', 'thick', 'thick_decim', 'x', 'y')
[x_pk, y_pk]				= deal(x, y);
load([dir_mat 'range_resolution.mat'], 'range_resolution')
FDM							= load([dir_mat 'gsfc_fdm_121.mat'], 'FAC', 'time_dt', 'x', 'y');
FDM.num_time_dt				= length(FDM.time_dt);
[FDM.x, FDM.y]				= deal(double(FDM.x), double(FDM.y));

core_avail					= true(1, num_core); % cores to include in dating, following order in core names, now all available with part of EGRIP

% re-order int_core_cat by core_ord
int_core_pk_unique			= unique(int_core_cat(:, 1));
for ii = 1:length(int_core_pk_unique)
	if (length(find(int_core_cat(:, 1) == int_core_pk_unique(ii))) > 1)
		ind_curr			= find(int_core_cat(:, 1) == int_core_pk_unique(ii));
		int_curr			= int_core_cat(ind_curr, 2);

	end
end

% simplify filenames and unsparsify depth (much simpler), and correct for firn
depth_pk					= cell(1, num_file_pk);
for ii = 1:num_file_pk
	file_pk{ii}				= file_pk{ii}(6:(end - 7));	
	depth_tmp				= full(depth{ii});
	depth_tmp(~depth_tmp)	= NaN;
	depth_pk{ii}			= NaN(num_layer(ii), num_trace(ii));
	depth_pk{ii}(:, ind_trace_layer{ii}) ...
							= depth_tmp;
	ind_FDM_curr			= interp1(FDM.time_dt, 1:FDM.num_time_dt, datetime(str2double(file_pk{ii}(1:8)), 'ConvertFrom', 'yyyymmdd'), 'nearest'); % time index for flight in FDM time vector
	depth_pk{ii}			= depth_pk{ii} + interp2(FDM.x, FDM.y, squeeze(FDM.FAC(:, :, ind_FDM_curr)), x_pk{ii}(ones(num_layer(ii), 1), :), y_pk{ii}(ones(num_layer(ii), 1), :), 'linear', 0); % add correction for firn air content (FAC)
end
clear depth

% % start core pool if possible
% if license('checkout', 'distrib_computing_toolbox')
%     pool_check              = gcp('nocreate');
%     if isempty(pool_check)
%         try
%             pool            = parpool('local', 4); % start 4 cores (no need to get crazy here)
%         catch
%             pool            = parpool('local');
%         end
%     end
%     num_pool                = pool.NumWorkers; % number of workers the pool (will be less than 4 if less than that available)
%     parallel_check          = true; % flag for parallelization
% else
    num_pool                = 0;
    parallel_check          = false;
% end

% load core depth-age scales
core                        = cell(1, num_core);
for ii = find(core_avail)
    core{ii}                = load([dir_mat 'depth_age/' name_core_short{ii} '_depth_age']);
    if isnan(core{ii}.age(end))
        [core{ii}.age, core{ii}.age_uncert, core{ii}.depth] ...
                            = deal(core{ii}.age(1:(end - 1)), core{ii}.age_uncert(1:(end - 1)), core{ii}.depth(1:(end - 1)));
    end
end

% index of NorthGRIP core, which is the backup core for age uncertainty
ind_ngrip                   = find(strcmp(name_core_short, 'ngrip'));

% Greenland-centered, EPSG:3413-projected 1-km grid limits
[x_min, x_max, y_min, y_max]= deal(-632e3, 846e3, -3344e3, -670e3);
[x_grd, y_grd]				= meshgrid(x_min:1e3:x_max, y_min:1e3:y_max);

% BedMachine v5
BM5                         = struct;
BM5.x                       = double(ncread('/Users/jamacgre/OneDrive - NASA/research/data/greenland/BedMachine/BedMachineGreenland-v5.nc', 'x'))'; % projected x, m
BM5.y                       = double(flipud(ncread('/Users/jamacgre/OneDrive - NASA/research/data/greenland/BedMachine/BedMachineGreenland-v5.nc', 'y'))); % projected y, m
BM5.thick					= rot90(ncread('/Users/jamacgre/OneDrive - NASA/research/data/greenland/BedMachine/BedMachineGreenland-v5.nc', 'thickness')); % ice thickness, m
thick_grd					= interp2(BM5.x, BM5.y, double(BM5.thick), x_grd, y_grd);
clear BM5

% load MAR accumulation data (https://arcticdata.io/catalog/view/doi%3A10.18739%2FA28G8FJ7F)
if (do_date || do_grd1 || (do_grd2 && do_nye_norm))
    MAR						= load([dir_mat 'greenland_mar_311_accum.mat'], 'accum', 'x', 'y');
    accum_grd               = MAR.accum(((MAR.y >= y_min) & (MAR.y <= y_max)), ((MAR.x >= x_min) & (MAR.x <= x_max)));
    clear MAR
end

% layer signal-to-noise ratio
if do_snr
    load([dir_mat 'snr_all.mat'], 'snr_all')
    % loop through segments to check for unavailable SNRs and assign them if necessary
    for ii = 1:num_file_pk
        if isempty(snr_all{ii})
            continue
        end
        snr_all{ii}(snr_all{ii} <= 0) ...
							= NaN;
	    snr_all{ii}(~isnan(snr_all{ii})) ...
							= 10 .^ (snr_all{ii}(~isnan(snr_all{ii})) ./ 10);
    end
else % assign all SNRs to the reference value (snr_ref) otherwise
    snr_all                 = cell(1, num_file_pk);
    for ii = 1:num_file_pk
		snr_all{ii}			= snr_ref .* ones(num_layer(ii), num_core);
    end
end

%%
if do_date
%%
    % initialize age-related cells
    date_counter            = 1;
    [age, age_core, age_old, age_ord, age_n, age_range, age_type, age_uncert] ...
                            = deal(cell(1, num_file_pk));
    for ii = 1:num_file_pk
        age_core{ii}		= NaN(num_layer(ii), num_core);
		[age{ii}, age_old{ii}, age_ord{ii}, age_n{ii}, age_range{ii}, age_type{ii}, age_uncert{ii}] ...
                            = deal(NaN(num_layer(ii), 1));
    end
    
	age_id_match			= id_match;
	
    % initialize reference age list
    num_age_ref				= length(layer_bin);
    age_ref					= NaN(num_age_ref, (num_core + 6));
    
    disp('Dating core-intersecting layers in core-intersecting segments...')
    
    % loop through merged picks file that intersect cores with available depth age scales
    for ii = find(core_avail(int_core_cat(:, 2)))
        
		[jj, kk, ll]		= deal(int_core_cat(ii, 1), int_core_cat(ii, 2), int_core_cat(ii, 3)); % jj: picks file index; kk: core index; ll: core intersection trace
		
        disp([file_pk{jj} '...' name_core_short{kk} '...'])
		
        % initialize uncertainty variables
        [age_uncert_depth, age_uncert_interp, age_uncert_radar] ...
							= deal(NaN(num_layer(jj), num_core));
		
        depth_curr			= mean(depth_pk{jj}(:, interp1(dist{jj}, 1:num_trace(jj), (dist{jj}(ll) - dist_int_max), 'nearest', 'extrap'):...
												   interp1(dist{jj}, 1:num_trace(jj), (dist{jj}(ll) + dist_int_max), 'nearest', 'extrap')), 2, 'omitnan'); % depths near core
        
        for mm = find(~isnan(depth_curr))' % loop through each of the merged file's layers, only interpolate age if layer present at/near core intersection
            
            age_core_curr	= interp1(core{kk}.depth, core{kk}.age, depth_curr(mm), 'spline'); % interpolated age at core intersection
			
			% age-overturned
			if ~age_check(age{jj}, age_uncert{jj}, mm, age_core_curr, age_uncert_frac_max, depth_pk{jj}, file_pk{jj})
				continue
			end
			
            % average if segment intersects the same ice core multiple times
            if isnan(age_core{jj}(mm, kk))
                age_core{jj}(mm, kk) ...
							= age_core_curr;
            else
                age_core{jj}(mm, kk) ...
							= mean([age_core_curr age_core{jj}(mm, kk)], 'omitnan');
            end
            
            % get signal-to-noise ratio
            if isnan(snr_all{jj}(mm, kk))
                snr_curr	= mean(snr_all{jj}(:, kk), 'omitnan');
                if isnan(snr_curr)
                    snr_curr= snr_ref;
                end
            else
                snr_curr	= snr_all{jj}(mm, kk);
            end
            
            % radar-induced uncertainty (range resolution / sqrt(snr))
            age_uncert_radar_curr ...
							= 0.5 * sum(abs(interp1(core{kk}.depth, core{kk}.age, ...
										(depth_curr(mm) + ([range_resolution(ind_campaign(jj)) -range_resolution(ind_campaign(jj))] ./ sqrt(snr_curr))), 'linear') - repmat(age_core{jj}(mm, kk), 1, 2)));
            if isnan(age_uncert_radar(mm, kk)) % same averaging if for uncertainty
                age_uncert_radar(mm, kk) ...
							= age_uncert_radar_curr;
            else
                age_uncert_radar(mm, kk) ...
							= mean([age_uncert_radar_curr age_uncert_radar(mm, kk)], 'omitnan');
            end
            
            % depth-induced uncertainty
            age_uncert_depth_curr ...
							= 0.5 * sum(abs(interp1(core{kk}.depth, core{kk}.age, (depth_curr(mm) + [depth_uncert -depth_uncert]), 'linear') - repmat(age_core{jj}(mm, kk), 1, 2)));
			if isnan(age_uncert_depth(mm, kk))
                age_uncert_depth(mm, kk) ...
							= age_uncert_depth_curr;
            else
                age_uncert_depth(mm, kk) ...
							= mean([age_uncert_depth_curr age_uncert_depth(mm, kk)], 'omitnan');
			end
        end
		
        % core-reported age uncertainty
        age_uncert_interp_curr ...
						= interp1(core{kk}.depth, core{kk}.age_uncert, depth_curr(~isnan(depth_curr)), 'linear');
        if isnan(age_uncert_radar(mm, kk))
            age_uncert_interp(~isnan(depth_curr), kk) ...
						= age_uncert_interp_curr;
        else
            age_uncert_interp(~isnan(depth_curr), kk) ...
						= mean([age_uncert_interp_curr age_uncert_interp(~isnan(depth_curr), kk)], 2, 'omitnan');
        end
        
        % order in which layers were dated
        age_ord{jj}(~isnan(depth_curr)) ...
							= date_counter;
        date_counter		= date_counter + 1;
        
        % assign core-intersecting ages, uncertainties, range, number of layers used (1) and type (0 for layers at closest trace in core-intersecting segments)
        age{jj}				= mean(age_core{jj}, 2, 'omitnan');
        age_uncert{jj}		= sqrt(mean((age_uncert_interp .^ 2), 2, 'omitnan') + mean((age_uncert_radar .^ 2), 2, 'omitnan') + mean((age_uncert_depth .^ 2), 2, 'omitnan')); % RSS of core age uncertainty, range resolution and depth uncertainty
        age_n{jj}(~isnan(age{jj})) ...
							= 1;
        age_range{jj}		= range(age_core{jj}, 2);
        age_type{jj}(~isnan(age{jj})) ...
							= 0;
		
        disp([num2str(length(find(~isnan(age{jj})))) '/' num2str(num_layer(jj)) ' layers dated.'])
    end
    
%%
    
    if do_match(1)
        
        disp('Assigning ages to matched core-intersecting layers...')
        
        [age_ref, age_core, age, age_uncert, age_range, age_n, age_type, age_ord, date_counter] ...
                            = match_age(1:num_age_ref, layer_bin, num_core, age_ref, age_core, age, age_uncert, age_range, age_n, age_type, age_ord, date_counter, true, 0, [0.5 1], do_age_check, age_uncert_frac_max, file_pk, depth_pk);
		
        disp('...done assigning ages to matched core-intersecting layers.')
    end
    
%%  
    
    % determine how many layers are currently dated and which segments
    [num_layer_dated_all, num_layer_dated_old] ...
                            = deal(0);
    for ii = 1:num_file_pk
        num_layer_dated_all = num_layer_dated_all + length(find(~isnan(age{ii})));
    end
    
    % initialize dating loop number
    num_date_loop           = 1;
    
    % keep going until dating does not improve anymore
    while ((num_layer_dated_all > num_layer_dated_old) && (num_date_loop <= num_date_loop_max))
        
        % don't bother if dating interpolation not requested
        if ~do_interp
            break
        end
        
        disp(['Dating loop #' num2str(num_date_loop) '...'])
        
        % loop through campaigns in preferred order based on quality
        for ii = campaign_ord
            
            if (ii == campaign_ord(1))
                disp('Dating overlapping layers in all segments...')
            end
            
            disp(['Campaign: ' campaign{ii} '...'])
            
            % loop through segments
            for jj = find(ind_campaign == ii)
                
                % move onto to next sub-segment if all or no layers were dated, or not enough dates, or no change from previous loop
                if (~any(isnan(age{jj})) || all(isnan(age{jj})) || (length(find(~isnan(age{jj}))) < 2) || isempty(setdiff(find(~isnan(age{jj})), find(~isnan(age_old{jj})))))
                    continue
                else
                    disp([file_pk{jj} '...'])
                end
                
                % Nye strain rate
                if strcmp(interp_type, 'quasi Nye')
                    accum_curr ...
                            = interp2(x_grd, y_grd, accum_grd, x_pk{jj}, y_pk{jj}, 'linear', NaN); % m/yr
                    strain_rate_curr ...
                            = accum_curr ./ thick{jj}; % start with Nye strain rate as initial guess
                    ind_nothick ...
                            = find(isnan(strain_rate_curr) | isinf(strain_rate_curr));
                    if ~isempty(ind_nothick)
                        strain_rate_curr(ind_nothick) ...
                            = accum_curr(ind_nothick) ./ interp2(x_grd, y_grd, thick_grd, x_pk{jj}(ind_nothick), y_pk{jj}(ind_nothick), 'linear', NaN);
                    end
                else
                    strain_rate_curr ...
                            = [];
                end
                
                % remember ages before attempting dating this time
                age_old{jj} = age{jj};
                
                % initial list of undated reflectors
                [ind_layer_undated, ind_layer_undated_ref] ...
							= deal(find(isnan(age{jj}))');
                
                % dummy list of undated reflectors to start first while loop
                ind_layer_undated_old ...
							= 1:num_layer(jj);
                
                disp([num2str(length(ind_layer_undated)) '/' num2str(num_layer(jj)) ' undated...'])
                
                % repeat dating until no more progress is made
                while (length(ind_layer_undated) < length(ind_layer_undated_old))
                    
                    ind_layer_ignore ...
							= [];
                    
                    % progressively assign ages to layers in core-intersecting segments that did not intersect an ice core
                    while any(isnan(age{jj}(setdiff(ind_layer_undated, ind_layer_ignore))))
                        [age{jj}, age_ord{jj}, age_n{jj}, age_range{jj}, age_type{jj}, age_uncert{jj}, date_counter, ind_layer_ignore] ...
							= date_interp(depth_pk{jj}, thick{jj}, age{jj}, age_ord{jj}, age_n{jj}, age_range{jj}, age_type{jj}, age_uncert{jj}, num_layer(jj), ind_layer_undated, ind_layer_ignore, thick_diff_max, ...
										  layer_diff_max, file_pk{jj}, interp_type, strain_rate_curr, parallel_check, num_pool, qn_tol, qn_iter_max, age_max, do_age_check, age_uncert_frac_max, age_uncert_rel_max, date_counter, (num_date_loop + 1));
                    end
                    
                    % updated undated layer set
                    ind_layer_undated_old ...
							= ind_layer_undated;
                    ind_layer_undated ...
							= find(isnan(age{jj}))';
                end
                
                % shallow quasi-Nye using surface as a last resort, if possible
                for ll = find(isnan(age{jj}))'
                    [age{jj}, age_ord{jj}, age_n{jj}, age_range{jj}, age_type{jj}, age_uncert{jj}, date_counter] ...
							= date_interp_shallow(depth_pk{jj}, thick{jj}, age{jj}, age_ord{jj}, age_n{jj}, age_range{jj}, age_type{jj}, age_uncert{jj}, interp_type, date_counter, ll, age_uncert_rel_max, depth_shallow_min, ...
												  depth_shallow_max, age_max, do_age_check, age_uncert_frac_max, file_pk{jj}, (num_date_loop + 1.25));
                end
                
                % display number of additional layers that were dated
                disp([num2str(length(ind_layer_undated_ref) - length(ind_layer_undated)) ' dated...'])
            end
        end
        
        disp('...done dating overlapping layers.')
        
        if do_match(2)
            
            disp('Assigning ages to matched layers...')
            
            [age_ref, age_core, age, age_uncert, age_range, age_n, age_type, age_ord, date_counter] ...
                            = match_age(find(isnan(age_ref(:, (num_core + 1))))', layer_bin, num_core, age_ref, age_core, age, age_uncert, age_range, age_n, age_type, age_ord, date_counter, false, (num_date_loop + 1.25), ...
										(num_date_loop + [1.5 1.75]), do_age_check, age_uncert_frac_max, file_pk, depth_pk);
            
            disp('...done assigning ages to matched layers.')
            
        end
        
        % reassign old number of dated layers
        num_layer_dated_old = num_layer_dated_all;
        
        % new number of dated layers
        num_layer_dated_all = 0;
        for ii = 1:num_file_pk
			num_layer_dated_all ...
                            = num_layer_dated_all + length(find(~isnan(age{ii})));
        end
        
        disp([num2str(num_layer_dated_all - num_layer_dated_old) ' dated layers added in this dating loop...'])
        if ((num_layer_dated_all - num_layer_dated_old) <= 0)
            disp('Dating loop shutting down...')
        end
        
        % increment dating loop counter
        num_date_loop       = num_date_loop + 1;
        
    end
    
    if do_save
        disp('Saving layer ages...')
		save([dir_mat 'date_all.mat'], '-v7.3', 'age', 'age_core', 'age_ref', 'age_id_match', 'age_ord', 'age_n', 'age_range', 'age_type', 'age_uncert', 'age_uncert_rel_max', 'dist_int_max', 'layer_bin')
    end
    
%%
elseif (do_grd1 || do_grd2)
%%
    disp('Loading layer ages...')
    load([dir_mat 'date_all.mat'], 'age', 'age_uncert')
    disp(['Loaded layer ages from ' dir_mat 'date_all.mat'])
%%
end

%%
if do_grd1
    
    disp('1-D inter/extrapolation of age at normalized depths and depth at specified ages...')
    
    % decimation
    if (num_depth_norm > 0)
        depth_norm          = ((1 / num_depth_norm):(1 / num_depth_norm):1)'; % normalized/relative depth vector
    else
        depth_norm          = 0;
    end
    
    % initialize along-track variables
	if (num_depth_norm > 0)
		[age_diff1, age_norm1, age_uncert_norm1, depth_norm_grd1, dist_grd1] ...
                            = deal(cell(1, num_file_pk));
	end
	if (num_age_iso > 0)
		[depth_iso1, depth_uncert_iso1] ...
                            = deal(cell(1, num_file_pk));
	end
    
    for ii = 1:num_file_pk
        
        % dated layer indices
        ind_layer_dated		= find(~isnan(age{ii}) & ~isnan(age_uncert{ii}));
        
        % not enough dated layers in this segment
        if (length(ind_layer_dated) < 2)
            continue
        else
            disp([file_pk{ii} '...'])
        end
        
        % decimated ice thickness and layer depths
        dist_decim			= NaN(1, num_decim(ii));
        depth_decim			= NaN(length(ind_layer_dated), num_decim(ii));
        for jj = 1:num_decim(ii)
            dist_decim(jj)	= mean(dist{ii}(ind_decim{ii}(jj):ind_decim{ii}(jj + 1)), 'omitnan');
            depth_decim(:, jj) ...
							= mean(depth_pk{ii}(ind_layer_dated, ind_decim{ii}(jj):ind_decim{ii}(jj + 1)), 2, 'omitnan');
        end
        
        % decimated Nye strain rate
        if strcmp(interp_type, 'quasi Nye')
            accum_curr		= interp2(x_grd, y_grd, accum_grd, x_pk{ii}, y_pk{ii}, 'linear', NaN);
            accum_decim		= NaN(1, num_decim(ii));
            for jj = 1:num_decim(ii)
                accum_decim(jj) ...
							= mean(accum_curr(ind_decim{ii}(jj):ind_decim{ii}(jj + 1)), 'omitnan');
            end
            strain_rate_decim ...
							= accum_decim ./ thick_decim{ii}; % start with Nye strain rate as initial guess
            strain_rate_decim(isnan(strain_rate_decim) | isinf(strain_rate_decim)) ...
							= mean(strain_rate_decim(~isinf(strain_rate_decim)), 'omitnan');
        end
        
        % thickness-normalized layer depths
        depth_decim_norm	= depth_decim ./ repmat(thick_decim{ii}, length(ind_layer_dated), 1);

		% segment-specific grid dimensions
		if (num_depth_norm > 0)
			[dist_grd1{ii}, depth_norm_grd1{ii}] ...
							= meshgrid(dist_decim, depth_norm);
		end
		
        % dated layers and uncertainties
        [age_dated, age_uncert_dated] ...
							= deal(age{ii}(ind_layer_dated), age_uncert{ii}(ind_layer_dated)); % only use dated layers
        
		% initialize
		if (num_depth_norm > 0)
			[age_norm1{ii}, age_uncert_norm1{ii}] ...
							= deal(NaN(num_depth_norm, num_decim(ii)));
		end
		if (num_age_iso > 0)
			[depth_iso1{ii}, depth_uncert_iso1{ii}] ...
							= deal(NaN(num_age_iso, num_decim(ii)));
		end
		
        % do age and depth uncertainties first
		for jj = find(~isnan(thick_decim{ii}) & (thick_decim{ii} > 0) & (sum(~isnan(depth_decim_norm), 1) > 1))
			
			if (num_depth_norm > 0)
				[depth_norm_unique, ind_depth_norm_unique]	...
							= unique(depth_decim_norm(~isnan(depth_decim_norm(:, jj)), jj));
			
				age_uncert_curr ...
							= age_uncert_dated(~isnan(depth_decim_norm(:, jj)));
				age_uncert_norm1{ii}(:, jj) ...
							= interp1(depth_norm_unique, age_uncert_curr(ind_depth_norm_unique), depth_norm, 'linear', 'extrap');
			end
			
			if (num_age_iso > 0)
				[age_unique, ind_age_unique] ...
							= unique(age_dated);
				depth_uncert_curr ...
							= mean(diff(core{ind_ngrip}.depth(interp1(core{ind_ngrip}.age, 1:length(core{ind_ngrip}.age), [(age_dated - age_uncert_dated) age_dated (age_dated + age_uncert_dated)], 'nearest', 'extrap')), 1, 2), 2, 'omitnan');
				depth_uncert_iso1{ii}(:, jj) ...
							= interp1(age_unique, depth_uncert_curr(ind_age_unique), age_iso, 'linear', 'extrap');
			end
		end
		
        % loop through each decimated trace and calculate 1D vertical age profile at uniform horizontal intervals in along-segment grid
        switch interp_type
            
            case 'linear' % NOT RECOMMENDED BECAUSE LINEAR LESS PHYSICALLY CONSISTENT WITH DEPTH-AGE PROFILES
				
				% skip if NaN/zero thickness or not enough dated layers
                for jj = find(~isnan(thick_decim{ii}) & (thick_decim{ii} > 0) & (sum(~isnan(depth_decim_norm), 1) > 1))
                    
                    % linearly interpolate age at normalized depths
                    if (num_depth_norm > 0)
						[depth_norm_unique, ind_depth_norm_unique] ...
							= unique(depth_decim_norm(~isnan(depth_decim_norm(:, jj)), jj));						
                        age_curr ...
                            = age_dated(~isnan(depth_decim_norm(:, jj)));
                        age_norm1{ii}(:, jj) ...
                            = interp1(depth_norm_unique, age_curr(ind_depth_norm_unique), depth_norm, 'linear', 'extrap');
                    end
                    
                    % linearly interpolate depth at specified ages
                    if (length(tmp1) > 1)
						[age_unique, ind_age_unique] ...
							= unique(age_dated(~isnan(depth_decim(:, jj))));
                        depth_decim_curr ...
                            = depth_decim(~isnan(depth_decim(:, jj)), jj);
                        depth_iso1{ii}(:, jj) ...
                            = interp1(age_unique, depth_decim_curr(ind_age_unique), age_iso, 'linear', 'extrap');
                    end
                    
                    depth_iso1{ii}((depth_iso1{ii}(:, jj) > thick_decim{ii}(jj)), jj) ...
                            = NaN;
                end
                
            case 'quasi Nye'
                
                for jj = find(~isnan(thick_decim{ii}) & (thick_decim{ii} > 0) & (sum(~isnan(depth_decim_norm), 1) > 1))
                    
                    % current normalized layer depths
					if (num_depth_norm > 0)
						[depth_curr, thick_diff_max_curr] ...
							= deal((depth_norm .* thick_decim{ii}(jj)), (thick_diff_max_iso * thick_decim{ii}(jj)));						
					end
                    [age_curr, depth_decim_curr] ...
							= deal(age_dated(~isnan(depth_decim(:, jj))), depth_decim(~isnan(depth_decim(:, jj)), jj));
                    [depth_decim_curr, ind_depth_ord] ...
							= sort(depth_decim_curr);
                    age_curr= age_curr(ind_depth_ord);
                    
                    for kk = 1:num_depth_norm
                        
                        [age_bound, depth_bound]			= deal([]);
                        
						% find two dated layers above and below current one
                        [ind_top, ind_bot]					= deal(find((depth_decim_curr < depth_curr(kk)), 2, 'last'), find((depth_decim_curr > depth_curr(kk)), 2));
                        
                        if (isempty(ind_top) && isempty(ind_bot))
                            continue
                        end
                        
                        if (isempty(ind_top) && (length(ind_bot) == 2)) % only two below
                            if all((depth_decim_curr(ind_bot(1)) - depth_curr(kk)) <= [(layer_diff_max_iso * diff(depth_decim_curr(ind_bot))) thick_diff_max_curr])
                                [age_bound, depth_bound]	= deal(age_curr(ind_bot), depth_decim_curr(ind_bot));
                            end
                        elseif (isempty(ind_bot) && (length(ind_top) == 2)) % only two above
                            if all((depth_curr(kk) - depth_decim_curr(ind_top(2))) <= [(layer_diff_max_iso * diff(depth_decim_curr(ind_top))) thick_diff_max_curr])
                                [age_bound, depth_bound]	= deal(age_curr(ind_top), depth_decim_curr(ind_top));
                            end
                        elseif (~isempty(ind_top) && ~isempty(ind_bot)) % at least one above and below (ideal)
                            [age_bound, depth_bound]		= deal(age_curr([ind_top(end); ind_bot(1)]), depth_decim_curr([ind_top(end); ind_bot(1)]));
                        end
                        
                        % quasi-Nye age
                        if (~isempty(depth_bound) && ~isempty(age_bound) && (diff(depth_bound) > 0) && (diff(age_bound) > 0))
                            age_norm1{ii}(kk, jj)			= date_quasi_nye(depth_bound, age_bound, (-diff(log(1 - (depth_bound ./ thick_decim{ii}(jj)))) / diff(age_bound)), depth_curr(kk), qn_tol, qn_iter_max, 'age');
                            % sanity check
                            if (((depth_curr(kk) > depth_bound(2)) && (age_norm1{ii}(kk, jj) < age_bound(2))) || ((depth_curr(kk) < depth_bound(1)) && (age_norm1{ii}(kk, jj) > age_bound(1))) || ...
                                (((depth_curr(kk) > depth_bound(1)) && (depth_curr(kk) < depth_bound(2))) && ((age_norm1{ii}(kk, jj) < age_bound(1)) || (age_norm1{ii}(kk, jj) > age_bound(2)))))
                                age_norm1{ii}(kk, jj)		= NaN;
                            end
                        end
                    end
                    
                    for kk = 1:num_age_iso
                        
                        [age_bound, depth_bound]			= deal([]);
                        
                        [ind_top, ind_bot]					= deal(find((age_curr < age_iso(kk)), 2, 'last'), find((age_curr > age_iso(kk)), 2));
                        
                        if (isempty(ind_top) && isempty(ind_bot))
                            continue
                        end
                        
                        if (isempty(ind_top) && (length(ind_bot) == 2))
                            if all((age_curr(ind_bot(1)) - age_iso(kk)) <= [(layer_diff_max_iso * diff(age_curr(ind_bot))) age_diff_max_iso])
                                [age_bound, depth_bound]	= deal(age_curr(ind_bot), depth_decim_curr(ind_bot));
                            end
                        elseif (isempty(ind_bot) && (length(ind_top) == 2))
                            if all((age_iso(kk) - age_curr(ind_top(2))) <= [(layer_diff_max_iso * diff(age_curr(ind_top))) age_diff_max_iso])
                                [age_bound, depth_bound]	= deal(age_curr(ind_top), depth_decim_curr(ind_top));
                            end
                        elseif (~isempty(ind_top) && ~isempty(ind_bot))
                            [age_bound, depth_bound]		= deal(age_curr([ind_top(end); ind_bot(1)]), depth_decim_curr([ind_top(end); ind_bot(1)]));
                        end
                        
                        % quasi-Nye isochrone depth
                        if (~isempty(depth_bound) && ~isempty(age_bound) && (diff(depth_bound) > 0) && (diff(age_bound) > 0))
                            depth_iso1{ii}(kk, jj)			= date_quasi_nye(depth_bound, age_bound, (-diff(log(1 - (depth_bound ./ thick_decim{ii}(jj)))) / diff(age_bound)), age_iso(kk), qn_tol, qn_iter_max, 'depth');
                            if (((age_iso(kk) > age_bound(2)) && (depth_iso1{ii}(kk, jj) < depth_bound(2))) || ((age_iso(kk) < age_bound(1)) && (depth_iso1{ii}(kk, jj) > depth_bound(1))) || ...
                                (((age_iso(kk) > age_bound(1)) && (age_iso(kk) < age_bound(2))) && ((depth_iso1{ii}(kk, jj) < depth_bound(1)) || (depth_iso1{ii}(kk, jj) > depth_bound(2)))))
                                depth_iso1{ii}(kk, jj)		= NaN;
                            end
                        end
                    end
                    depth_iso1{ii}((depth_iso1{ii}(:, jj) > thick_decim{ii}(jj)), jj) ...
							= NaN;
                end
        end
        
        % NaN out misbehaving values
		if (num_depth_norm > 0)
			age_norm1{ii}(isinf(age_norm1{ii}) | (age_norm1{ii} < 0) | (age_norm1{ii} > age_max)) ...
							= NaN;
			age_uncert_norm1{ii}(isinf(age_uncert_norm1{ii}) | (age_uncert_norm1{ii} < 0)) ...
							= NaN;
		end
		if (num_age_iso > 0)
			depth_iso1{ii}(isinf(depth_iso1{ii}) | (depth_iso1{ii} < 0)) ...
							= NaN;
			depth_uncert_iso1{ii}(depth_uncert_iso1{ii} < 0) ...
							= NaN;
		end
		
        % reference Nye age field
        if (num_depth_norm > 0)
            age_nye_tmp		= [-(1 ./ strain_rate_decim(ones((num_depth_norm - 1), 1), :)) .* log(1 - depth_norm(1:(end - 1), ones(1, num_decim(ii)))); NaN(1, num_decim(ii))];
            age_diff1{ii}	= (age_norm1{ii} - age_nye_tmp) ./ age_nye_tmp;
        end
    end
    
    disp('...done 1-D gridding layer ages.')
%%
    if do_save
        disp('Saving 1-D age grids...')

		if (num_depth_norm > 0)
			save([dir_mat 'age_grd1_norm.mat'], '-v7.3', 'age_diff1', 'age_norm1', 'age_uncert_norm1', 'depth_norm', 'depth_norm', 'depth_norm_grd1', 'num_depth_norm')
		end

		if (num_age_iso > 0)

			save([dir_mat 'age_grd1_iso.mat'], '-v7.3', 'age_iso', 'depth_iso1', 'depth_uncert_iso1', 'num_age_iso')
			[depth_iso1_cell, depth_uncert_iso1_cell] ...
						= deal(depth_iso1, depth_uncert_iso1);

			[depth_iso1, depth_uncert_iso1] ...
						= deal(cell2struct(depth_iso1, 'depth_iso1'), cell2struct(depth_uncert_iso1, 'depth_uncert_iso1'));
			save([dir_mat 'age_grd1_iso_struct.mat'], 'age_iso', 'depth_iso1', 'depth_uncert_iso1', 'num_age_iso')
			[depth_iso1, depth_uncert_iso1] ...
						= deal(depth_iso1_cell, depth_uncert_iso1_cell);
			
%%			
			[x_pk_decim, y_pk_decim] ...
						= deal(cell(1, num_file_pk));
			for ii = find(cellfun(@length, depth_iso1))
				[x_pk_decim{ii}, y_pk_decim{ii}] ...
						= deal(x_pk{ii}(ind_decim_mid{ii}), y_pk{ii}(ind_decim_mid{ii}));
			end
			[x_pk_decim, y_pk_decim, thick_decim_all] ...
						= deal([x_pk_decim{:}]', [y_pk_decim{:}]', [thick_decim{:}]');
			[depth_iso1_split1, depth_uncert_iso1_split1] ...
						= deal(cell(num_age_iso, num_file_pk));
			for ii = 1:num_age_iso
				for jj = find(cellfun(@length, depth_iso1))
					[depth_iso1_split1{ii, jj}, depth_uncert_iso1_split1{ii, jj}] ...
						= deal(depth_iso1{jj}(ii, :), depth_uncert_iso1{jj}(ii, :));
				end
			end
			[depth_iso1_split2, depth_uncert_iso1_split2] ...
						= deal(cell(1, num_age_iso));
			for ii = 1:num_age_iso
				[depth_iso1_split2{ii}, depth_uncert_iso1_split2] ...
						= deal([depth_iso1_split1{ii, :}]', [depth_uncert_iso1_split1{ii, :}]');
			end
%%			
			for ii = 1:num_age_iso
				ind_good= find(~isnan(depth_iso1_split2{ii}) & ~isnan(thick_decim_all));
				str_tmp = num2str((1e-3 * age_iso(ii)), '%1.1f');
				if strcmp(str_tmp((end - 1):end), '.0')
					str_tmp ...
						= str_tmp(1:(end - 2));
				end
				writetable(table(x_pk_decim(ind_good), y_pk_decim(ind_good), depth_iso1_split2{ii}(ind_good), depth_uncert_iso1_split2(ind_good), thick_decim_all(ind_good), 'VariableNames', {'x' 'y' 'depth' 'uncertainty' 'thickness'}), ...
						   [dir_mat 'depth_iso1_' str_tmp '_ka.csv'])
			end
		end
		
        disp(['Saved 1-D age grids as ' dir_mat 'age_grd1_*.mat.'])
    end
%%
elseif do_grd2
%%
    disp('Loading 1-D age gridding...')
    load([dir_mat 'age_grd1.mat'], 'age_norm1', 'age_uncert_norm1', 'depth_iso1', 'depth_norm', 'depth_uncert_iso1', 'num_age_iso', 'num_depth_norm')
    if do_nye_norm
        load([dir_mat 'age_grd1.mat'], 'age_diff1')
    end
    disp('Loaded 1-D age gridding from age_grd1.mat')
%%
end

if do_grd2
    
    disp('2-D gridding of normalized ages, isochrone depths and their uncertainties...')
    
    % prepare cells for ages/depths to be gridded
    [age_norm2, age_uncert_norm2] ...
                            = deal(NaN(size(x_grd, 1), size(x_grd, 2), num_depth_norm));
    [depth_iso2, depth_uncert_iso2] ...
                            = deal(NaN(size(x_grd, 1), size(x_grd, 2), num_age_iso));
    [thick_cat, x_all_age_cat, x_all_age_uncert_cat, x_all_depth_cat, x_all_depth_uncert_cat, y_all_age_cat, y_all_age_uncert_cat, y_all_depth_cat, y_all_depth_uncert_cat] ...
                            = deal([]);
    [age_all, age_all_min, age_all_max, age_uncert_all, x_all_age, x_all_age_min, x_all_age_max, x_all_age_uncert, y_all_age, y_all_age_min, y_all_age_max, y_all_age_uncert] ...
                            = deal(cell(1, num_depth_norm));
    if do_nye_norm
        age_diff_all        = cell(1, num_depth_norm);
    end
    [depth_all, depth_all_min, depth_all_max, depth_uncert_all, x_all_depth, x_all_depth_min, x_all_depth_max, x_all_depth_uncert, y_all_depth, y_all_depth_min, y_all_depth_max, y_all_depth_uncert] ...
                            = deal(cell(1, num_age_iso));
    
    % loop through each 1-D set and add its values to the mix
    for ii = 1:num_file_pk
        if ~isempty(age_norm1{ii})
            [x_all_age_cat, y_all_age_cat] ...
                    = deal([x_all_age_cat; x_pk{ii}(ind_decim_mid{ii})'], [y_all_age_cat; y_pk{ii}(ind_decim_mid{ii})']);
            for jj = 1:num_depth_norm
                age_all{jj} ...
                    = [age_all{jj}; age_norm1{ii}(jj, :)'];
                if do_nye_norm
                    age_diff_all{jj} ...
                        = [age_diff_all{jj}; age_diff1{ii}(jj, :)'];
                end
                [age_min_tmp, age_max_tmp] ...
                    = deal(age_norm1{ii}(jj, :));
                for kk = find(isnan(age_norm1{ii}(jj, :)))
                    if ~isempty(find(age_norm1{ii}(1:(jj - 1), kk), 1, 'last'))
                        age_min_tmp(jj) ...
                            = age_norm1{ii}(find(age_norm1{ii}(1:(jj - 1), kk), 1, 'last'), kk);
                    end
                    if ~isempty(find(age_norm1{ii}((jj + 1):end, kk), 1))
                        age_max_tmp(jj) ...
                            = age_norm1{ii}(find(age_norm1{ii}((jj + 1):end, kk), 1), kk);
                    end
                end
                [age_all_min{jj}, age_all_max{jj}] ...
                    = deal([age_all_min{jj}; age_min_tmp'], [age_all_max{jj}; age_max_tmp']);
            end
        end
        if ~isempty(age_uncert_norm1{ii})
            [x_all_age_uncert_cat, y_all_age_uncert_cat] ...
                    = deal([x_all_age_uncert_cat; x_pk{ii}(ind_decim_mid{ii})'], [y_all_age_uncert_cat; y_pk{ii}(ind_decim_mid{ii})']);
            for jj = 1:num_depth_norm
                age_uncert_all{jj} ...
                    = [age_uncert_all{jj}; age_uncert_norm1{ii}(jj, :)'];
            end
        end
        if ~isempty(depth_iso1{ii})
            [x_all_depth_cat, y_all_depth_cat, thick_cat] ...
                    = deal([x_all_depth_cat; x_pk{ii}(ind_decim_mid{ii})'], [y_all_depth_cat; y_pk{ii}(ind_decim_mid{ii})'], [thick_cat; double(thick_decim{ii}')]);
            for jj = 1:num_age_iso
                depth_all{jj} ...
                    = [depth_all{jj}; depth_iso1{ii}(jj, :)'];
                [depth_min_tmp, depth_max_tmp] ...
                        = deal(depth_iso1{ii}(jj, :));
                for kk = find(isnan(depth_iso1{ii}(jj, :)))
                    if ~isempty(find(depth_iso1{ii}(1:(jj - 1), kk), 1, 'last'))
                        depth_min_tmp(jj) ...
                            = depth_iso1{ii}(find(depth_iso1{ii}(1:(jj - 1), kk), 1, 'last'), kk);
                    end
                    if ~isempty(find(depth_iso1{ii}((jj + 1):end, kk), 1))
                        depth_max_tmp(jj) ...
                            = depth_iso1{ii}(find(depth_iso1{ii}((jj + 1):end, kk), 1), kk);
                    end
                end
                [depth_all_min{jj}, depth_all_max{jj}] ...
                    = deal([depth_all_min{jj}; depth_min_tmp'], [depth_all_max{jj}; depth_max_tmp']);
            end
        end
        if ~isempty(depth_uncert_iso1{ii})
            [x_all_depth_uncert_cat, y_all_depth_uncert_cat] ...
                    = deal([x_all_depth_uncert_cat; x_pk{ii}(ind_decim_mid{ii})'], [y_all_depth_uncert_cat; y_pk{ii}(ind_decim_mid{ii})']);
            for jj = 1:num_age_iso
                depth_uncert_all{jj} ...
                    = [depth_uncert_all{jj}; depth_uncert_iso1{ii}(jj, :)'];
            end
        end
    end
    
    % concatenate and remove the bad stuff prior to gridding
    for ii = 1:num_depth_norm
        [x_all_age{ii}, x_all_age_min{ii}, x_all_age_max{ii}] ...
                            = deal(x_all_age_cat);
        [y_all_age{ii}, y_all_age_min{ii}, y_all_age_max{ii}] ...
                            = deal(y_all_age_cat);
        ind_good            = find(~isnan(age_all{ii}));
        [x_all_age{ii}, y_all_age{ii}, age_all{ii}] ...
                            = deal(x_all_age{ii}(ind_good), y_all_age{ii}(ind_good), age_all{ii}(ind_good));
        if do_nye_norm
            age_diff_all{ii}= age_diff_all{ii}(ind_good);
        end
        [xy_all, ind_unique]= unique([x_all_age{ii} y_all_age{ii}], 'rows');
        [x_all_age{ii}, y_all_age{ii}, age_all{ii}] ...
                            = deal(xy_all(:, 1), xy_all(:, 2), age_all{ii}(ind_unique));
        if do_nye_norm
            age_diff_all{ii}= age_diff_all{ii}(ind_unique);
        end
        ind_good            = find(~isnan(age_all_min{ii}));
        [x_all_age_min{ii}, y_all_age_min{ii}, age_all_min{ii}] ...
                            = deal(x_all_age_min{ii}(ind_good), y_all_age_min{ii}(ind_good), age_all_min{ii}(ind_good));
        [xy_all, ind_unique]= unique([x_all_age_min{ii} y_all_age_min{ii}], 'rows');
        [x_all_age_min{ii}, y_all_age_min{ii}, age_all_min{ii}] ...
                            = deal(xy_all(:, 1), xy_all(:, 2), age_all_min{ii}(ind_unique));
        ind_good            = find(~isnan(age_all_max{ii}));
        [x_all_age_max{ii}, y_all_age_max{ii}, age_all_max{ii}] ...
                            = deal(x_all_age_max{ii}(ind_good), y_all_age_max{ii}(ind_good), age_all_max{ii}(ind_good));
        [xy_all, ind_unique]= unique([x_all_age_max{ii} y_all_age_max{ii}], 'rows');
        [x_all_age_max{ii}, y_all_age_max{ii}, age_all_max{ii}] ...
                            = deal(xy_all(:, 1), xy_all(:, 2), age_all_max{ii}(ind_unique));
        [x_all_age_uncert{ii}, y_all_age_uncert{ii}] ...
                            = deal(x_all_age_uncert_cat, y_all_age_uncert_cat);
        ind_good            = find(~isnan(age_uncert_all{ii}));
        [x_all_age_uncert{ii}, y_all_age_uncert{ii}, age_uncert_all{ii}] ...
                            = deal(x_all_age_uncert{ii}(ind_good), y_all_age_uncert{ii}(ind_good), age_uncert_all{ii}(ind_good));
        [xy_all, ind_unique]= unique([x_all_age_uncert{ii} y_all_age_uncert{ii}], 'rows');
        [x_all_age_uncert{ii}, y_all_age_uncert{ii}, age_uncert_all{ii}] ...
                            = deal(xy_all(:, 1), xy_all(:, 2), age_uncert_all{ii}(ind_unique));
    end
    
    for ii = 1:num_age_iso
        [x_all_depth{ii}, x_all_depth_min{ii}, x_all_depth_max{ii}] ...
                            = deal(x_all_depth_cat);
        [y_all_depth{ii}, y_all_depth_min{ii}, y_all_depth_max{ii}] ...
                            = deal(y_all_depth_cat);
        [thick_all_curr, thick_all_min_curr, thick_all_max_curr] ...
                            = deal(thick_cat);
        ind_good            = find(~isnan(depth_all{ii}) & ~isnan(thick_all_curr));
        [x_all_depth{ii}, y_all_depth{ii}, depth_all{ii}, thick_all_curr] ...
                            = deal(x_all_depth{ii}(ind_good), y_all_depth{ii}(ind_good), depth_all{ii}(ind_good), thick_all_curr(ind_good));
        [xy_all, ind_unique]= unique([x_all_depth{ii} y_all_depth{ii}], 'rows', 'stable');
        [x_all_depth{ii}, y_all_depth{ii}, depth_all{ii}, thick_all_curr] ...
                            = deal(xy_all(:, 1), xy_all(:, 2), depth_all{ii}(ind_unique), thick_all_curr(ind_unique));
        depth_all{ii}       = depth_all{ii} ./ thick_all_curr;
        ind_good            = find(~isnan(depth_all_min{ii}) & ~isnan(thick_all_min_curr));
        [x_all_depth_min{ii}, y_all_depth_min{ii}, depth_all_min{ii}, thick_all_min_curr] ...
                            = deal(x_all_depth_min{ii}(ind_good), y_all_depth_min{ii}(ind_good), depth_all_min{ii}(ind_good), thick_all_min_curr(ind_good));
        [xy_all, ind_unique]= unique([x_all_depth_min{ii} y_all_depth_min{ii}], 'rows', 'stable');
        [x_all_depth_min{ii}, y_all_depth_min{ii}, depth_all_min{ii}, thick_all_min_curr] ...
                            = deal(xy_all(:, 1), xy_all(:, 2), depth_all_min{ii}(ind_unique), thick_all_min_curr(ind_unique));
        depth_all_min{ii}   = depth_all_min{ii} ./ thick_all_min_curr;
        ind_good            = find(~isnan(depth_all_max{ii}) & ~isnan(thick_all_max_curr));
        [x_all_depth_max{ii}, y_all_depth_max{ii}, depth_all_max{ii}, thick_all_max_curr] ...
                            = deal(x_all_depth_max{ii}(ind_good), y_all_depth_max{ii}(ind_good), depth_all_max{ii}(ind_good), thick_all_max_curr(ind_good));
        [xy_all, ind_unique]= unique([x_all_depth_max{ii} y_all_depth_max{ii}], 'rows', 'stable');
        [x_all_depth_max{ii}, y_all_depth_max{ii}, depth_all_max{ii}, thick_all_max_curr] ...
                            = deal(xy_all(:, 1), xy_all(:, 2), depth_all_max{ii}(ind_unique), thick_all_max_curr(ind_unique));
        depth_all_max{ii}   = depth_all_max{ii} ./ thick_all_max_curr;
        [x_all_depth_uncert{ii}, y_all_depth_uncert{ii}, thick_uncert_curr] ...
                            = deal(x_all_depth_uncert_cat, y_all_depth_uncert_cat, thick_cat);
        ind_good            = find(~isnan(depth_uncert_all{ii}) & ~isnan(thick_uncert_curr));
        [x_all_depth_uncert{ii}, y_all_depth_uncert{ii}, depth_uncert_all{ii}, thick_uncert_curr] ...
                            = deal(x_all_depth_uncert{ii}(ind_good), y_all_depth_uncert{ii}(ind_good), depth_uncert_all{ii}(ind_good), thick_uncert_curr(ind_good));
        [xy_all, ind_unique]= unique([x_all_depth_uncert{ii} y_all_depth_uncert{ii}], 'rows', 'stable');
        [x_all_depth_uncert{ii}, y_all_depth_uncert{ii}, depth_uncert_all{ii}, thick_uncert_curr] ...
                            = deal(xy_all(:, 1), xy_all(:, 2), depth_uncert_all{ii}(ind_unique), thick_uncert_curr(ind_unique));
        depth_uncert_all{ii}= depth_uncert_all{ii} ./ thick_uncert_curr;
    end
    
    clear *_cat
    
    % use Nye normalization and Gaussian fits to remove major anomalies from age_norm1 sets to be gridded
    if do_nye_norm
        for ii = 1:(num_depth_norm - 1)
            fit_curr        = fit((min(age_diff_all{ii}):0.1:max(age_diff_all{ii}))', histcounts(age_diff_all{ii}, min(age_diff_all{ii}):0.1:max(age_diff_all{ii}))', ['gauss' num2str(num_gauss(ii))]); % Gaussian fit to histogram of age differences from Nye model
            switch num_gauss(ii)
                case 1
                    age_diff_range ...
                            = norminv([0.05 0.99], fit_curr.b1, fit_curr.c1); % lower/upper bounds of acceptable Nye differences
                case 2
                    age_diff_range ...
                            = zeros(2);
                    age_diff_range(1, :) ...
                            = norminv([0.05 0.99], fit_curr.b1, fit_curr.c1);
                    age_diff_range(2, :) ...
                            = norminv([0.05 0.99], fit_curr.b2, fit_curr.c2);
                    age_diff_range ...
                            = [min(age_diff_range(:, 1)) max(age_diff_range(:, 2))]; % take min/max of the two fits' lower/upper bounds, respectively
            end
            ind_good        = find((age_diff_all{ii} >= age_diff_range(1)) & (age_diff_all{ii} <= age_diff_range(2))); % indices of ages with Nye differences within spec
            [x_all_age{ii}, y_all_age{ii}, age_all{ii}, age_diff_all{ii}] ...
                            = deal(x_all_age{ii}(ind_good), y_all_age{ii}(ind_good), age_all{ii}(ind_good), age_diff_all{ii}(ind_good));
            [xy_all, ~, ind_common] ...
                            = intersect([x_all_age{ii} y_all_age{ii}], [x_all_age_min{ii} y_all_age_min{ii}], 'rows', 'stable');
            [x_all_age_min{ii}, y_all_age_min{ii}, age_all_min{ii}] ...
                            = deal(xy_all(:, 1), xy_all(:, 2), age_all_min{ii}(ind_common));
            [xy_all, ~, ind_common] ...
                            = intersect([x_all_age{ii} y_all_age{ii}], [x_all_age_max{ii} y_all_age_max{ii}], 'rows', 'stable');
            [x_all_age_max{ii}, y_all_age_max{ii}, age_all_max{ii}] ...
                            = deal(xy_all(:, 1), xy_all(:, 2), age_all_max{ii}(ind_common));
        end
    end
%%    
    if parallel_check
        disp('...fixed-depth layers...')
        parfor ii = 1:num_depth_norm
            disp([num2str(ii) '/' num2str(num_depth_norm) '...'])
            try
                age_interpolant ...
                            = scatteredInterpolant(x_all_age{ii}, y_all_age{ii}, age_all{ii}, 'natural', 'none');
                age_norm2(:, :, ii) ...
                            = age_interpolant(x_grd, y_grd);
                age_uncert_interpolant ...
                            = scatteredInterpolant(x_all_age_uncert{ii}, y_all_age_uncert{ii}, age_uncert_all{ii}, 'natural', 'none');
                age_uncert_norm2(:, :, ii) ...
                            = age_uncert_interpolant(x_grd, y_grd);
            catch
                disp([num2str(depth_norm(ii)) ' isopach failed...'])
            end
        end
        disp('...fixed-age layers...')
        parfor ii = 1:num_age_iso
            disp([num2str(ii) '/' num2str(num_age_iso) '...'])
            try
                depth_interpolant ...
                            = scatteredInterpolant(x_all_depth{ii}, y_all_depth{ii}, depth_all{ii}, 'natural', 'none');
                depth_iso2(:, :, ii) ...
                            = depth_interpolant(x_grd, y_grd) .* thick_grd;
                depth_uncert_interpolant ...
                            = scatteredInterpolant(x_all_depth_uncert{ii}, y_all_depth_uncert{ii}, depth_uncert_all{ii}, 'natural', 'none');
                depth_uncert_iso2(:, :, ii) ...
                            = depth_uncert_interpolant(x_grd, y_grd) .* thick_grd;
            catch
                disp([num2str(age_iso(ii)) '-ka isochrone failed...'])
            end
        end
    else
        disp('...fixed-depth layers...')
        for ii = 1:num_depth_norm
            disp([num2str(ii) '/' num2str(num_depth_norm) '...'])
            try
                age_interpolant ...
                            = scatteredInterpolant(x_all_age{ii}, y_all_age{ii}, age_all{ii}, 'natural', 'none');
                age_norm2(:, :, ii) ...
                            = age_interpolant(x_grd, y_grd);
                age_uncert_interpolant ...
                            = scatteredInterpolant(x_all_age_uncert{ii}, y_all_age_uncert{ii}, age_uncert_all{ii}, 'natural', 'none');
                age_uncert_norm2(:, :, ii) ...
                            = age_uncert_interpolant(x_grd, y_grd);
            catch
                disp([num2str(depth_norm(ii)) ' isopach failed...'])
            end
        end
        disp('...fixed-age layers...')
        for ii = 1:num_age_iso
            disp([num2str(ii) '/' num2str(num_age_iso) '...'])
            try
                depth_interpolant ...
                            = scatteredInterpolant(x_all_depth{ii}, y_all_depth{ii}, depth_all{ii}, 'natural', 'none');
                depth_iso2(:, :, ii) ...
                            = depth_interpolant(x_grd, y_grd) .* thick_grd;
                depth_uncert_interpolant ...
                            = scatteredInterpolant(x_all_depth_uncert{ii}, y_all_depth_uncert{ii}, depth_uncert_all{ii}, 'natural', 'none');
                depth_uncert_iso2(:, :, ii) ...
                            = depth_uncert_interpolant(x_grd, y_grd) .* thick_grd;
            catch
                disp([num2str(age_iso(ii)) '-ka isochrone failed...'])
            end
        end
    end
    
    % trim bad ages/depths
    age_norm2((age_norm2 < 0) | (age_norm2 > age_max)) ...
                            = NaN;
    age_uncert_norm2(age_uncert_norm2 < 0) ...
                            = NaN;
    depth_iso2(depth_iso2 < 0) ...
                            = NaN;
    depth_iso2(depth_iso2 > thick_grd(:, :, ones(1, 1, num_age_iso))) ...
                            = NaN;
    depth_uncert_iso2(depth_uncert_iso2 < 0) ...
                            = NaN;
    
    disp('...done 2-D gridding age.')
%%
    if do_save
        disp('Saving 2-D age grids...')
    	save([dir_mat 'age_grd2.mat'], '-v7.3', 'age_norm2', 'age_uncert_norm2', 'age_iso', 'depth_iso2', 'depth_norm', 'depth_uncert_iso2', 'num_age_iso', 'num_depth_norm', 'x_grd', 'y_grd')
    	disp('Saved 2-D age grids as mat/age_grd2.mat.')
    end
%%
end

if parallel_check
    delete(pool)
end