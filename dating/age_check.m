function age_test           = age_check(age, age_uncert, curr_layer, age_candidate, age_uncert_frac_max, depth, file_pk)
% AGE_CHECK Test whether a layer age would lead to age overturning (decreasing with increasing depth) prior to assigning it.
% 
% To be used within DATE_LAYERS only.
% 
% Joe MacGregor (NASA)
% Last updated: 22 July 2024

% assume candidate age will pass
age_test                    = true;

% skip if only one layer or no dated layers
if (isscalar(age) || isempty(find(~isnan(age), 1)))
    return
end

% assign candidate age for testing
[age(curr_layer), age_uncert(curr_layer)] ...
							= deal(age_candidate, -Inf);

% determine dated layers
ind_layer_good              = find(~isnan(age));

if isscalar(ind_layer_good)
    return
end

% index of current layer
ind_curr_layer              = find(ind_layer_good == curr_layer);

% trim age to dated layers
[age, age_uncert]			= deal(age(ind_layer_good), age_uncert(ind_layer_good));

% sort layers and dated layers by increasing depth only where current layer exists
[depth_sort, ind_depth_sort]= sort(depth(ind_layer_good, ~isnan(depth(curr_layer, :))), 1);

% sort ages and NaN out empty ones
[age_sort, age_uncert_sort] = deal(age(ind_depth_sort), age_uncert(ind_depth_sort));
[age_sort(isnan(depth_sort)), age_uncert_sort(isnan(depth_sort))] ...
							= deal(NaN);

% age difference normalized by uncertainty between layers sorted by depth, column-wise
age_diff_rel				= diff(age_sort) ./ max(cat(3, age_uncert_sort(1:(end - 1), :), age_uncert_sort(2:end, :)), [], 3, 'omitnan');

% review if negative age difference exceeds maximum allowable fraction of existing age uncertainty
% if ~isempty(find((age_diff <= 0), 1))
if ~isempty(find((age_diff_rel <= -age_uncert_frac_max), 1))
	
    [i_bad, j_bad]          = find(age_diff_rel <= -age_uncert_frac_max); % row/column indices of overturns (column vectors)
	
    % check if current layer is part of the overturn
    if ~isempty(find((ind_depth_sort(sub2ind(size(ind_depth_sort), [i_bad(:); min([(i_bad(:) + 1) (size(ind_depth_sort, 1) .* ones(length(j_bad), 1))], [], 2)], repmat(j_bad(:), 2, 1))) == ind_curr_layer), 1)) 
		% figure;imagesc(1e-3 .* diff(age_sort));colorbar				
		% figure;imagesc(1e-3 .* max(cat(3, age_uncert_sort(1:(end - 1), :), age_uncert_sort(2:end, :)), [], 3, 'omitnan'));colorbar		
		% figure;imagesc(1e2 .* age_diff_rel);colorbar
		% disp(min(1e2 .* age_diff_rel, [], 'all'))
        disp(['!!! ' file_pk ' #' num2str(curr_layer) ' would be significantly age-overturned so not dated here.'])
        age_test            = false;
        return
    end
end