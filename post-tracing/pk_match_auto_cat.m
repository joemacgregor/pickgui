% PK_MATCH_AUTO_CAT Automatically determine matching layers within merged radar segments.
% 
% Joe MacGregor (NASA)
% Last updated: 18 July 2024

clear

do_save                     = true;

dir_mat						= '/Users/jamacgre/OneDrive - NASA/research/matlab/pickgui_v2/mat/';
dist_int_max                = 0.5e3; % +/- range to extract local layer depths, m
depth_diff_fact				= 3; % factor by which to multiply range resolution for depth match

load([dir_mat 'xyz_all.mat'], 'num_campaign')
load([dir_mat 'range_resolution.mat'], 'range_resolution')
load([dir_mat 'int_all_cat.mat'], 'int_all')

pk_cat						= load([dir_mat 'pk_cat.mat'], 'depth', 'dist', 'file_pk', 'ind_campaign', 'ind_trace_layer', 'num_file_pk', 'num_layer', 'num_trace');
% simplify filenames and unsparsify depth (simpler)
[depth_pk, file_pk]			= deal(cell(1, pk_cat.num_file_pk));
for ii = 1:pk_cat.num_file_pk
	file_pk{ii}				= pk_cat.file_pk{ii}(6:(end - 7));
	depth_tmp				= full(pk_cat.depth{ii});
	depth_tmp(~depth_tmp)	= NaN;
	depth_pk{ii}			= NaN(pk_cat.num_layer(ii), pk_cat.num_trace(ii));
	depth_pk{ii}(:, pk_cat.ind_trace_layer{ii}) ...
							= depth_tmp;
end
clear depth_tmp
pk_cat						= rmfield(pk_cat, 'depth');

% do not permit matches between campaigns at far end of system from each other, i.e., ICORDSv1 cannot be matched with MCoRDSv2
campaign_match              = true(num_campaign);
campaign_match(1:4, 16:end) = false;
campaign_match(16:end, 1:4) = false;

%% Populate fence log

% initialize layer match list
layer_match_list			= NaN(0, 4);
[depth_diff_match, int_match] ...
							= deal([]);
num_match					= zeros(1, size(int_all, 1));

disp('Auto-populating fenced list...')

% first evaluate each segment intersection for initial set of layer matches
for ii = 1:size(int_all, 1)
	
	[jj, kk]				= deal(int_all(ii, 1), int_all(ii, 4));
	
	if ~campaign_match(pk_cat.ind_campaign(jj), pk_cat.ind_campaign(kk))
        continue
	end
	
	if ~isempty(find((pk_cat.num_trace([jj kk]) == 0), 1))
		disp('EMPTY MERGED FILE')
		continue
	end
	
	ind_curr1				= interp1(pk_cat.dist{jj}, 1:pk_cat.num_trace(jj), (pk_cat.dist{jj}(int_all(ii, 2)) + [-dist_int_max dist_int_max]), 'nearest', 'extrap');
	ind_curr2				= interp1(pk_cat.dist{kk}, 1:pk_cat.num_trace(kk), (pk_cat.dist{kk}(int_all(ii, 5)) + [-dist_int_max dist_int_max]), 'nearest', 'extrap');
		
	% extract depths for current pk_cat traces and unsparsify
	depth_curr1				= mean(depth_pk{jj}(:, ind_curr1(1):ind_curr1(2)), 2, 'omitnan');
	depth_curr2				= mean(depth_pk{kk}(:, ind_curr2(1):ind_curr2(2)), 2, 'omitnan');
	
	% skip if no layers to match
	if (isempty(find(~isnan(depth_curr1), 1)) || isempty(find(~isnan(depth_curr2), 1)))
		continue
	end
	
	% set threshold depth difference to factor of minimum range resolution of layer pair's respective campaigns
	depth_diff_max_test		= depth_diff_fact * max(range_resolution(pk_cat.ind_campaign([jj kk])));
	ind_remain				= find(~isnan(depth_curr2));
	
	for ll = find(~isnan(depth_curr1))'
		
		% layer index in intersecting (second) segment of the pair that matches the primary (first) segment layer depth
		[depth_diff_min, curr_layer] ...
							= min(abs(depth_curr1(ll) - depth_curr2(ind_remain)));
		
		% only keep going if a layer pair meets the current maximum range resolution
		if (depth_diff_min > depth_diff_max_test)
    		continue
		end
		
		layer_match_list	= [layer_match_list; int_all(ii, 1) ll int_all(ii, 4) ind_remain(curr_layer(1))]; %#ok<AGROW>
		depth_diff_match	= [depth_diff_match; depth_diff_min]; %#ok<AGROW>
		int_match			= [int_match; ii]; %#ok<AGROW>
		
		ind_remain(curr_layer(1)) ...
							= []; % don't let another match with same layer occur
		num_match(ii)		= num_match(ii) + 1; % increment number of matches
		
		if isempty(ind_remain)
			break
		end
	end
	
	if num_match(ii)
		disp([num2str(ii) '/' num2str(size(int_all, 1)) ' ... ' pk_cat.file_pk{jj} ' + ' pk_cat.file_pk{kk} ' ... ' num2str(num_match(ii)) ' match(es).'])
	end
end

% trim to unique layer matches and ensure sorted
[layer_match_list, ind_unique] ...
							= unique(layer_match_list, 'rows', 'stable');
[layer_match_list, ind_sort]= sortrows(layer_match_list, [1 3 2 4]);
depth_diff_match			= depth_diff_match(ind_unique(ind_sort));
int_match					= int_match(ind_unique(ind_sort));

%% 

% now verify that each layer match is valid at ALL intersections of that segment pair (easier to do after building the match list)
[depth_diff_match_test, depth_diff_max_ref] ...
							= deal(zeros(size(layer_match_list, 1), 1));

for ii = 1:size(layer_match_list, 1)
	
	if ~mod(ii, 1e3)
		disp(ii)
	end
	
	% current segment pair (file_pk)
	[jj, kk]				= deal(layer_match_list(ii, 1), layer_match_list(ii, 3));
	
	% all intersections of current segment pair
	int_curr				= int_all(((int_all(:, 1) == jj) & (int_all(:, 4) == kk)), [2 5]);
	
	% loop through all the file's intersections and get layer depth difference there
	depth_diff_max_curr		= NaN(1, size(int_curr, 1));
	
	for ll = 1:size(int_curr, 1)
		% get current horizontal indices of intersection and environs
		ind_curr1			= interp1(pk_cat.dist{jj}, 1:pk_cat.num_trace(jj), (pk_cat.dist{jj}(int_curr(ll, 1)) + [-dist_int_max dist_int_max]), 'nearest', 'extrap');
		ind_curr2			= interp1(pk_cat.dist{kk}, 1:pk_cat.num_trace(kk), (pk_cat.dist{kk}(int_curr(ll, 2)) + [-dist_int_max dist_int_max]), 'nearest', 'extrap');
		depth_diff_max_curr(ll)	...
							= abs(mean(depth_pk{jj}(layer_match_list(ii, 2), ind_curr1(1):ind_curr1(2)), 'omitnan') - mean(depth_pk{kk}(layer_match_list(ii, 4), ind_curr2(1):ind_curr2(2)), 'omitnan'));
	end
	% preserve maximum depth difference
	depth_diff_match_test(ii) ...
							= max(depth_diff_max_curr, [], 'omitnan');
	depth_diff_max_ref(ii)	= depth_diff_fact * max(range_resolution(pk_cat.ind_campaign([jj kk])));
end
%%
disp(['Removing ' num2str(length(find(depth_diff_match_test > depth_diff_max_ref))) ' matches due to too high depth difference at other intersections...'])

% remove those with large depth differences
layer_match_list((depth_diff_match_test > depth_diff_max_ref), :) ...
							= [];
depth_diff_match(depth_diff_match_test > depth_diff_max_ref) ...
							= [];
int_match(depth_diff_match_test > depth_diff_max_ref) ...
							= [];

%%

if do_save
    save([dir_mat 'layer_match_list_auto.mat'], '-v7.3', 'layer_match_list', 'depth_diff_match', 'int_match')
	disp('SAVED layer_match_list_auto.mat')
end

% % NO LONGER NEEDED
% 
% layer_match_list_pk			= layer_match_list; % copy and then adjust numbers to be relevant to pk files
% 
% % loop through pk files that lost layers and adjust layer numbers
% for ii = find(pk_cat.num_layer_removed)
% 	disp(pk_cat.file_pk{ii})
% 	ind_curr				= find((layer_match_list_pk(:, 1) == ii) | (layer_match_list_pk(:, 3) == ii)); % current indices for pk file with removed layers
% 	ind_layer_remap			= setdiff(1:(pk_cat.num_layer(ii) + pk_cat.num_layer_removed(ii)), pk_cat.ind_layer_removed{ii}); % remapping to layers in original pk files
% 	layer_match_list_pk(ind_curr(layer_match_list_pk(ind_curr, 1) == ii), 2) ...
% 							= ind_layer_remap(layer_match_list_pk(ind_curr(layer_match_list_pk(ind_curr, 1) == ii), 2)); % remap layer numbers in first column
% 	layer_match_list_pk(ind_curr(layer_match_list_pk(ind_curr, 3) == ii), 4) ...
% 							= ind_layer_remap(layer_match_list_pk(ind_curr(layer_match_list_pk(ind_curr, 3) == ii), 4)); % now second column
% end