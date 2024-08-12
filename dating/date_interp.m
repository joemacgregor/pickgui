function [age, age_ord, age_n, age_range, age_type, age_uncert, date_counter, ind_layer_ignore] ...
                            = date_interp(depth, thick, age, age_ord, age_n, age_range, age_type, age_uncert, num_layer, ind_layer_undated, ind_layer_ignore, thick_diff_max, layer_diff_max, file_pk, interp_type, strain_rate_ref, parallel_check, num_pool, tol, iter_max, age_max, do_age_check, age_uncert_frac_max, age_uncert_rel_max, date_counter, age_type_assign)
% DATE_INTERP Date an undated layer.
%
% To be used within DATE_LAYERS only.
%
% Joe MacGregor (NASA)
% Last updated: 22 July 2024

% determine which traces are potentially usable because they have at least two dated layers
ind_trace_usable            = find(sum(~isnan(depth(~isnan(age), :)), 1) > 1);

% stop if insufficient usable traces for overlapping dating
if isempty(ind_trace_usable)
    ind_layer_ignore		= find(isnan(age))';
    return
end

% loop through undated layers to find number of overlapping traces with dated layers
num_overlap                 = zeros(num_layer, 1);
for ii = intersect(find(isnan(age))', setdiff(ind_layer_undated, ind_layer_ignore), 'stable')
    num_overlap(ii)         = length(find(~isnan(depth(ii, ind_trace_usable))));
end

% stop if there are no more undated layers that overlap with dated layers
if isempty(find((num_overlap > 1), 1))
    ind_layer_ignore		= find(isnan(age))';
    return
end

% layer with most overlap
[~, curr_layer]             = max(num_overlap);

% along-transect overlapping indices of layer with most overlap
ind_overlap                 = ind_trace_usable(~isnan(depth(curr_layer, ind_trace_usable)));
num_overlap_max             = length(ind_overlap);

% prepare for loop iteration by slicing all necessary variables
[age_bound, age_uncert_bound, depth_bound] ...
                            = deal(NaN(2, num_overlap_max));

% initially include all overlapping layers and then whittle them down to the bounding layers
[age_bound_tmp, age_uncert_bound_tmp, depth_bound_tmp, depth_curr, ind_overlap_bounded, thick_diff_max_curr, thick_overlap, strain_rate_overlap] ...
                            = deal(age(~isnan(age)), age_uncert(~isnan(age)), depth(~isnan(age), ind_overlap), depth(curr_layer, ind_overlap), false(1, num_overlap_max), (thick(ind_overlap) .* thick_diff_max), thick(ind_overlap), strain_rate_ref(ind_overlap));

% if thickness unavailable, then set thickness criterion to infinity (rare for deep radar, always for accumulation radar)
thick_diff_max_curr(isnan(thick_diff_max_curr)) ...
                            = Inf;

% sort dated layers by increasing depth
[depth_bound_tmp, ind_depth_ord] ...
                            = sort(depth_bound_tmp, 1);
[age_bound_tmp, age_uncert_bound_tmp] ...
                            = deal(age_bound_tmp(:, ones(1, num_overlap_max)), age_uncert_bound_tmp(:, ones(1, num_overlap_max)));
[age_bound_tmp, age_uncert_bound_tmp] ...
                            = deal(age_bound_tmp(ind_depth_ord), age_uncert_bound_tmp(ind_depth_ord));

% trim depths/ages to upper/lower bounds
for ii = 1:num_overlap_max
    
    % upper/lower bounding layer indices
    [ind_top, ind_bot]      = deal(find((depth_bound_tmp(:, ii) < depth_curr(ii)), 2, 'last'), find((depth_bound_tmp(:, ii) > depth_curr(ii)), 2));
    
    % skip if no bounding dated reflectors
    if (isempty(ind_top) && isempty(ind_bot))
        continue
    end
    
    % bin ages and depths that survive criteria if unbounded above or below (not both)
    if (isempty(ind_top) && (length(ind_bot) == 2))
        if all((depth_bound_tmp(ind_bot(1), ii) - depth_curr(ii)) <= [(layer_diff_max * diff(depth_bound_tmp(ind_bot, ii))) thick_diff_max_curr(ii)])
            [age_bound(:, ii), age_uncert_bound(:, ii), depth_bound(:, ii), ind_overlap_bounded(ii)] ...
                            = deal(age_bound_tmp(ind_bot, ii), age_uncert_bound_tmp(ind_bot, ii), depth_bound_tmp(ind_bot, ii), true);
        end
    elseif (isempty(ind_bot) && (length(ind_top) == 2))
        if all((depth_curr(ii) - depth_bound_tmp(ind_top(2), ii)) <= [(layer_diff_max * diff(depth_bound_tmp(ind_top, ii))) thick_diff_max_curr(ii)])
            [age_bound(:, ii), age_uncert_bound(:, ii), depth_bound(:, ii), ind_overlap_bounded(ii)] ...
                            = deal(age_bound_tmp(ind_top, ii), age_uncert_bound_tmp(ind_top, ii), depth_bound_tmp(ind_top, ii), true);
        end
    elseif (~isempty(ind_top) && ~isempty(ind_bot))
        [age_bound(:, ii), age_uncert_bound(:, ii), depth_bound(:, ii), ind_overlap_bounded(ii)] ...
                            = deal(age_bound_tmp([ind_top(end); ind_bot(1)], ii), age_uncert_bound_tmp([ind_top(end); ind_bot(1)], ii), depth_bound_tmp([ind_top(end); ind_bot(1)], ii), true);
    end
end

% skip traces with bad age bounds
if ~isempty(find((diff(age_bound) <= 0), 1))
    ind_overlap_bounded(diff(age_bound) <= 0) ...
                            = false;
    disp([file_pk ' layer #' num2str(curr_layer) ' is sometimes badly age-bounded (' num2str(1e-3 * mean(diff(age_bound(:, (diff(age_bound) <= 0))), 'omitnan')) ' kyr mean top/bottom difference).'])
end

% remove unbounded or poorly dated traces
if isempty(find(ind_overlap_bounded, 1))
    ind_layer_ignore		= [ind_layer_ignore curr_layer];
    return
elseif any(~ind_overlap_bounded)
    [age_bound, age_uncert_bound, depth_bound, ind_overlap, depth_curr, thick_overlap, strain_rate_overlap] ...
                            = deal(age_bound(:, ind_overlap_bounded), age_uncert_bound(:, ind_overlap_bounded), depth_bound(:, ind_overlap_bounded), ind_overlap(ind_overlap_bounded), depth_curr(ind_overlap_bounded), thick_overlap(ind_overlap_bounded), strain_rate_overlap(ind_overlap_bounded));
    num_overlap_max         = length(ind_overlap);
end
age_overlap                 = NaN(1, num_overlap_max);

% date layer using either linear or quasi-Nye dating
switch interp_type
    
    case 'linear'
        
        % parallelized loop through each trace where current layer exists, linearly inter/extrapolating age
        if (parallel_check && (num_overlap_max > (5 * num_pool)))
            parfor ii = 1:num_overlap_max
                age_overlap(ii) ...
                            = interp1(depth_bound(:, ii), age_bound(:, ii), depth_curr(ii), 'linear', 'extrap');
            end
        else
            % non-parallel loop doing same thing
            for ii = 1:num_overlap_max
                age_overlap(ii) ...
                            = interp1(depth_bound(:, ii), age_bound(:, ii), depth_curr(ii), 'linear', 'extrap');
            end
        end
        
    case 'quasi Nye'
        
        % reference/starting strain rate, supplemented with reference/Nye strain rate where necessary
        strain_rate_curr    = -(diff(log(1 - (depth_bound ./ thick_overlap(ones(2, 1), :))))) ./ diff(age_bound);
        strain_rate_curr(isnan(strain_rate_curr) | isinf(strain_rate_curr)) ...
                            = strain_rate_overlap(isnan(strain_rate_curr) | isinf(strain_rate_curr));
        
        if (parallel_check && (num_overlap_max >= (5 * num_pool)))
            parfor ii = 1:num_overlap_max
                age_overlap(ii) ...
                            = date_quasi_nye(depth_bound(:, ii), age_bound(:, ii), strain_rate_curr(ii), depth_curr(ii), tol, iter_max, 'age');
            end
        else
            for ii = 1:num_overlap_max
                age_overlap(ii) ...
                            = date_quasi_nye(depth_bound(:, ii), age_bound(:, ii), strain_rate_curr(ii), depth_curr(ii), tol, iter_max, 'age');
            end
        end
end

% sanity checks
age_overlap(isinf(age_overlap) | (age_overlap < 0) | (age_overlap > age_max) | ((depth_curr > depth_bound(2, :)) & (age_overlap < age_bound(2, :))) | ((depth_curr < depth_bound(1, :)) & (age_overlap > age_bound(1, :))) | ...
            (((depth_curr > depth_bound(1, :)) & (depth_curr < depth_bound(2, :))) & ((age_overlap < age_bound(1, :)) | age_overlap > age_bound(2, :)))) ...
                            = NaN;

% mean
age_overlap_mean            = mean(age_overlap, 2, 'omitnan');

% check if new age causes overturning
if isnan(age_overlap_mean)
    ind_layer_ignore		= [ind_layer_ignore curr_layer];
    return
elseif (do_age_check && ~age_check(age, age_uncert, curr_layer, age_overlap_mean, age_uncert_frac_max, depth, file_pk))
    ind_layer_ignore		= [ind_layer_ignore curr_layer];
    return
end

% weighted uncertainty as the sum of uncertainties of both bounding layers
age_overlap_uncert          = sqrt(var(age_overlap, 'omitnan') + mean(mean((age_uncert_bound .^ 2), 1, 'omitnan'), 2, 'omitnan'));

% limit age assignment to well-constrained layers
if (isempty(find(~isnan(age_overlap), 1)) || ((age_overlap_uncert / age_overlap_mean) > age_uncert_rel_max) || (age_overlap_mean <= 0))
    ind_layer_ignore		= [ind_layer_ignore curr_layer];
    return
end

% breakout age information for newly dated layer
age(curr_layer)             = age_overlap_mean;
age_uncert(curr_layer)      = age_overlap_uncert;
age_ord(curr_layer)         = date_counter;
date_counter                = date_counter + 1;
age_n(curr_layer)           = 1;
age_range(curr_layer)       = range(age_overlap);
age_type(curr_layer)        = age_type_assign;