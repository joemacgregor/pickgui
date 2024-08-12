function [age_ref, age_core, age, age_uncert, age_range, age_n, age_type, age_ord, date_counter] ...
                            = match_age(ind_layer_curr, layer_bin, num_core, age_ref, age_core, age, age_uncert, age_range, age_n, age_type, age_ord, date_counter, do_core, age_type_max, age_type_assign, do_age_check, age_uncert_frac_max, file_pk, depth)
% MATCH_AGE Match ages of dated and undated layers.
% 
% To be used within DATE_LAYERS only.
% 
% Joe MacGregor (NASA)
% Last updated: 22 July 2024

% loop through layers to consider
for ii = ind_layer_curr
    
    age_cat                 = [];
    id                      = layer_bin{ii};
	
    % loop through all layers matched to current one and test whether they are dated; concatenate them if so
    for jj = 1:size(id, 1)
        if (~isnan(age{id(jj, 1)}(id(jj, 2))) && (age_type{id(jj, 1)}(id(jj, 2)) <= age_type_max))
            if do_core % use core ages if match_age being called right after those were dated
                age_cat     = [age_cat; age_core{id(jj, 1)}(id(jj, 2), :) age_uncert{id(jj, 1)}(id(jj, 2))]; %#ok<AGROW>
            else
                age_cat     = [age_cat; age{id(jj, 1)}(id(jj, 2)) age_uncert{id(jj, 1)}(id(jj, 2))]; %#ok<AGROW>
            end
        end
    end
    
    % pass to next iteration if no dated layers
    if isempty(age_cat)
        continue
    end
    
    % break out concatenated dated layers
    age_ref_tmp				= NaN(1, (num_core + 6));
    if do_core
        age_ref_tmp(1:num_core) ...
                            = mean(age_cat(:, 1:num_core), 1, 'omitnan');
        age_ref_tmp(num_core + 1) ...
                            = mean(age_ref_tmp(1:num_core), 'omitnan');
        age_ref_tmp(num_core + 2) ...
                            = mean(age_cat(:, (num_core + 1)), 'omitnan');
        age_ref_tmp(num_core + 3) ...
                            = range(age_cat(:, (num_core + 1)));
    else
        age_ref_tmp(num_core + 1) ...
                            = mean(age_cat(:, 1), 'omitnan');
        age_ref_tmp(num_core + 2) ...
                            = mean(age_cat(:, 2), 'omitnan');
        age_ref_tmp(num_core + 3) ...
                            = range(age_cat(:, 1));
    end
    age_ref_tmp(num_core + 4) ...
                            = size(age_cat, 1);
    age_ref_tmp(num_core + 5) ...
                            = age_type_assign(1);
    age_ref_tmp(num_core + 6) ...
                            = date_counter; % dating order increment
    
    % assign best age set to master layer
    age_ref(ii, :)			= age_ref_tmp;
    
    % reassign ages of matched layers to that of reference layer, depending on status of the matched layer and if age does not create overturning
    for jj = 1:size(id, 1)
        if (do_age_check && ~age_check(age{id(jj, 1)}, age_uncert{id(jj, 1)}, id(jj, 2), age_ref_tmp(num_core + 1), age_uncert_frac_max, depth{id(jj, 1)}, file_pk{id(jj, 1)}))
            continue
        end
        if isnan(age{id(jj, 1)}(id(jj, 2)))
            age_type{id(jj, 1)}(id(jj, 2)) ...
                            = age_type_assign(2);
        else
            age_type{id(jj, 1)}(id(jj, 2)) ...
                            = age_type_assign(1);
        end
        if do_core
            age_core{id(jj, 1)}(id(jj, 2), :) ...
                            = age_ref_tmp(1:num_core);
        end
        [age{id(jj, 1)}(id(jj, 2)), age_uncert{id(jj, 1)}(id(jj, 2)), age_range{id(jj, 1)}(id(jj, 2)), age_n{id(jj, 1)}(id(jj, 2)), age_ord{id(jj, 1)}(id(jj, 2))] ...
                            = deal(age_ref_tmp(num_core + 1), age_ref_tmp(num_core + 2), age_ref_tmp(num_core + 3), age_ref_tmp(num_core + 4), age_ref_tmp(num_core + 6));
    end
    
    date_counter            = date_counter + 1; % dating order increment
end