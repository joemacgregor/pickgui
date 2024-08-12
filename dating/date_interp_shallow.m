function [age, age_ord, age_n, age_range, age_type, age_uncert, date_counter] ...
                            = date_interp_shallow(depth, thick, age, age_ord, age_n, age_range, age_type, age_uncert, interp_type, date_counter, curr_layer, age_uncert_rel_max, depth_shallow_min, depth_shallow_max, age_max, do_age_check, age_uncert_frac_max, file_pk, age_type_assign)
% DATE_INTERP_SHALLOW Date a shallow undated layer.
% 
% To be used within DATE_LAYERS only.
% 
% Joe MacGregor (NASA)
% Last updated: 24 July 2024

% determine which traces are potentially usable because they have at least one dated layer
ind_overlap                 = find((sum(~isnan(depth(~isnan(age), :)), 1) > 0) & ~isnan(depth(curr_layer, :)));

% stop if insufficient usable traces for overlapping dating
if isempty(ind_overlap)
    return
end

% stop if insufficient usable traces for overlapping dating
if isempty(ind_overlap)
    return
end

% along-transect overlapping indices of layer with most overlap
num_overlap_max             = length(ind_overlap);

% prepare for loop iteration by slicing all necessary variables
[age_bound, age_uncert_bound, depth_bot] ...
                            = deal(NaN(1, num_overlap_max));

% initially include all layers and then whittle them down to the bounding layers
[age_bound_tmp, age_uncert_slice_tmp, depth_bound_tmp, depth_curr, ind_overlap_bounded, depth_min, depth_max, thick_curr] ...
                            = deal(age(~isnan(age)), age_uncert(~isnan(age)), depth(~isnan(age), ind_overlap), depth(curr_layer, ind_overlap), false(1, num_overlap_max), (depth_shallow_min * thick(ind_overlap)), (depth_shallow_max * thick(ind_overlap)), thick(ind_overlap));

% sort dated layers by increasing depth
[depth_bound_tmp, ind_depth_ord] ...
                            = sort(depth_bound_tmp, 1);
[age_bound_tmp, age_uncert_slice_tmp] ...
                            = deal(age_bound_tmp(:, ones(1, num_overlap_max)), age_uncert_slice_tmp(:, ones(1, num_overlap_max)));
[age_bound_tmp, age_uncert_slice_tmp] ...
                            = deal(age_bound_tmp(ind_depth_ord), age_uncert_slice_tmp(ind_depth_ord));

% trim depths/ages to upper/lower bounds
for ii = 1:num_overlap_max
    
    [ind_top, ind_bot]      = deal(find((depth_bound_tmp(:, ii) < depth_curr(ii)), 1, 'last'), find((depth_bound_tmp(:, ii) > depth_curr(ii)), 1));
    
    % only want dated reflectors below, not above
    if (~isempty(ind_top) || isempty(ind_bot))
        continue
    end
    
    % bin good ages and depths, following criterion for maximum depth based on thickness
    if (all([depth_curr(ii) depth_bound_tmp(ind_bot, ii)] <= depth_max(ii)) && (depth_curr(ii) >= depth_min(ii)))
        [age_bound(ii), age_uncert_bound(ii), depth_bot(ii), ind_overlap_bounded(ii)] ...
                            = deal(age_bound_tmp(ind_bot, ii), age_uncert_slice_tmp(ind_bot, ii), depth_bound_tmp(ind_bot, ii), true);
    end
end

% remove unbounded or poorly dated traces
if isempty(find(ind_overlap_bounded, 1))
    return
elseif any(~ind_overlap_bounded)
    [age_bound, age_uncert_bound, depth_bot, depth_curr, thick_curr] ...
                            = deal(age_bound(ind_overlap_bounded), age_uncert_bound(ind_overlap_bounded), depth_bot(ind_overlap_bounded), depth_curr(ind_overlap_bounded), thick_curr(ind_overlap_bounded));
end

switch interp_type
    case 'linear'
        age_overlap         = age_bound .* (depth_curr ./ depth_bot);
    case 'quasi Nye'
        % quasi-Nye dating for shallow reflectors
        age_overlap         = age_bound .* (log(1 - (depth_curr ./ thick_curr)) ./ log(1 - (depth_bot ./ thick_curr)));
end

% sanity checks
age_overlap(isinf(age_overlap) | (age_overlap < 0) | (age_overlap > age_max) | (age_overlap > age_bound)) ...
                            = NaN;

% mean overlapping age
age_overlap_mean            = mean(age_overlap, 2, 'omitnan');

% check if mean causes overturning
if (do_age_check && ~age_check(age, age_uncert, curr_layer, age_overlap_mean, age_uncert_frac_max, depth, file_pk))
    return
end

% weighted uncertainty as sum of uncertainties of both bounding layers
age_overlap_uncert          = sqrt(var(age_overlap, 'omitnan') + mean((age_uncert_bound .^ 2), 2, 'omitnan'));

% limit age assignment to well-constrained layers
if (isempty(find(~isnan(age_overlap), 1)) || ((age_overlap_uncert / age_overlap_mean) > age_uncert_rel_max) || (age_overlap_mean <= 0))
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