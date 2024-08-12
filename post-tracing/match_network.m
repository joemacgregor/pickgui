% MATCH_NETWORK Testing layer match network for overlap sequentially.
% 
% Joe MacGregor (NASA)
% Last updated: 19 July 2024

clear

do_save						= false;
path_cat					= '/Users/jamacgre/OneDrive - NASA/data/gris_strat2/arctic_cat_r1/';

dir_mat						= '/Users/jamacgre/OneDrive - NASA/research/matlab/pickgui_v2/mat/';
load([dir_mat 'xyz_all.mat'], 'campaign', 'num_campaign')
load([dir_mat 'pk_cat.mat'], 'depth', 'dist', 'file_pk', 'ind_trace_layer', 'num_file_pk', 'num_layer', 'num_trace', 'x', 'y')
[x_pk, y_pk]				= deal(x, y);
clear x y
load([dir_mat 'int_all_cat.mat'], 'int_all')
load([dir_mat 'layer_match_list.mat'], 'layer_match_list', 'depth_diff_match', 'int_match')

% BedMachine v5
BM5                         = struct;
BM5.x                       = double(ncread('/Users/jamacgre/OneDrive - NASA/research/data/greenland/BedMachine/BedMachineGreenland-v5.nc', 'x'))'; % projected x, m
BM5.y                       = double(flipud(ncread('/Users/jamacgre/OneDrive - NASA/research/data/greenland/BedMachine/BedMachineGreenland-v5.nc', 'y'))); % projected y, m
BM5.elev_surf               = double(rot90(ncread('/Users/jamacgre/OneDrive - NASA/research/data/greenland/BedMachine/BedMachineGreenland-v5.nc', 'surface'))); % surface elevation (geoid), m
BM5.thick					= double(rot90(ncread('/Users/jamacgre/OneDrive - NASA/research/data/greenland/BedMachine/BedMachineGreenland-v5.nc', 'thickness'))); % ice thickness, m
decim_bm					= 5;
[BM5.elev_surf, BM5.thick, BM5.x, BM5.y] ...
							= deal(BM5.elev_surf(1:decim_bm:end, 1:decim_bm:end), BM5.thick(1:decim_bm:end, 1:decim_bm:end), BM5.x(1:decim_bm:end), BM5.y(1:decim_bm:end));

% simplify filenames and unsparsify depth/elevation (much simpler)
[depth_pk_mean, depth_pk, elev_pk] ...
							= deal(cell(1, num_file_pk));
for ii = 1:num_file_pk
	file_pk{ii}				= file_pk{ii}(6:(end - 7));
	depth_tmp				= full(depth{ii});
	depth_tmp(~depth_tmp)	= NaN;
	depth_pk{ii}			= NaN(num_layer(ii), num_trace(ii));
	depth_pk{ii}(:, ind_trace_layer{ii}) ...
							= depth_tmp;
	depth_pk_mean{ii}		= mean(depth_pk{ii}, 2, 'omitnan');
	elev_pk{ii}				= interp2(BM5.x, BM5.y, BM5.elev_surf, x_pk{ii}, y_pk{ii}) - depth_pk{ii};	
end
clear depth depth_tmp

speed_vacuum                = 299792458; % m/s
permitt_ice                 = 3.15; % dimensionless
speed_ice                   = speed_vacuum / sqrt(permitt_ice); % m/s

ind_bed_cutoff				= 10; % number of samples below which to blank the bed in radargrams
dist_smooth					= 100e3; % smoothing distance for median subtraction, m
num_trace_blank_min			= 5; % minimum number of traces with non-NaN values to consider before blanking row
num_color					= 2 ^ 8; % number of colors for display
num_std_color				= 2; % number of standard deviations for color range
sample_decim				= 2; % number of samples to decimate in fast time (vertical)
dist_decim					= 0.5e3; % decimation interval, m

% layer color range for display
colors_def                  = [0    0       0.75;
                               0    0       1;
                               0    0.25    1;
                               0    0.50    1;
                               0    0.75    1;
                               0    1       1;
                               0.25 1       0.75;
                               0.50 1       0.50;
                               0.75 1       0.25;
                               1    1       0;
                               1    0.75    0;
                               1    0.50    0;
                               1    0.25    0;
                               1    0       0;
                               0.75 0       0];
num_color_def				= size(colors_def, 1);

%%

disp('Assigning layer IDs sequentially...')

num_layer_match				= size(layer_match_list, 1);
layer_match_bin				= NaN(num_layer_match, 2);
% layer_match_test			= true(num_layer_match, 1);

% simplify matching set
layer_unique				= unique([layer_match_list(:, 1:2); layer_match_list(:, 3:4)], 'rows');
layer_match_simple			= NaN(num_layer_match, 2);
for ii = 1:size(layer_unique, 1)
    % if ~mod(ii, 1e3)
    %     disp(ii)
    % end
    layer_match_simple(ismember(layer_match_list(:, 1:2), layer_unique(ii, :), 'rows'), 1) ...
                            = ii;
    layer_match_simple(ismember(layer_match_list(:, 3:4), layer_unique(ii, :), 'rows'), 2) ...
                            = ii;
end

disp('Binning layers via intersections...')

% intial reference id
bin_id_curr					= 1; % starting bin id is 1
diff_threshold				= 0; % depth difference above which to assess overlap, m

for ii = 1:num_layer_match
	
	if ~mod(ii, 1e3)
		disp([num2str(ii) '/' num2str(num_layer_match)])
	end
	
	% if no bin ID assigned to either layer then plan to add a new one
	if all(isnan(layer_match_bin(ii, :)))
		new_bin_flag		= true;
	else
		new_bin_flag		= false;
	end
	
	if new_bin_flag % make a new bin
		
		id_match_curr		= bin_id_curr;
		
		% % build layer set for overlap testing before adding match
		% layer_bin_test      = [layer_match_list(ii, 1:2); layer_match_list(ii, 3:4)];
		
	else % current matching IDs
		
    	id_match_curr		= sort(layer_match_bin(ii, :));
		
		% layer_bin_test      = unique([layer_match_list((layer_match_bin(:, 1) == id_match_curr(1)), 1:2); layer_match_list((layer_match_bin(:, 2) == id_match_curr(1)), 3:4); layer_match_list(ii, 1:2); layer_match_list(ii, 3:4)], 'rows');
		% 
    	% % add higher ID bin's layers to test bin 
		% if (~isnan(id_match_curr(2)) && ~isequal(id_match_curr(1), id_match_curr(2)))
        % 	layer_bin_test	= unique([layer_bin_test; layer_match_list((layer_match_bin(:, 1) == id_match_curr(2)), 1:2); layer_match_list((layer_match_bin(:, 2) == id_match_curr(2)), 3:4)], 'rows');
		% end
	end
	
	% overlap_flag			= false;
	% 
	% % unique segments contained in current test bin
	% seg_unique				= unique(layer_bin_test(:, 1));
	% 
	% % evaluate all segments within bin
	% for jj = 1:length(seg_unique)
	% 
	% 	% indices within layer_bin that from current segment
	% 	ind_layer_risk		= find(layer_bin_test(:, 1) == seg_unique(jj));
	% 
    % 	% skip if only one layer from that segment within bin (no overlap risk)
	% 	if isscalar(ind_layer_risk)
    %     	continue
	% 	end
	% 
	% 	% loop through layers in current segment, searching for layers that overlap horizontally (>diff_threshold)
	% 	for kk = 1:(length(ind_layer_risk) - 1)
	% 		for ll = (kk + 1):length(ind_layer_risk)
    %     		diff_depth_curr ...
    %                 		= abs(diff(depth_pk{seg_unique(jj)}(layer_bin_test(ind_layer_risk([kk ll]), 2), :)));
    %     		if ~isempty(find((diff_depth_curr > diff_threshold), 1)) % matched layers at different depths (>0 permits layers at same depth)
    %         		overlap_flag ...
	% 						= true;
	% 				break
    %     		end
	% 		end
	% 		if overlap_flag
	% 			break
	% 		end
	% 	end
	% 	if overlap_flag
	% 		break
	% 	end
	% end
	
	% handle overlap_flag, if not add a new bin or merge bins
	% if overlap_flag
	% 	layer_match_test(ii)= false;
	% 	continue
	if new_bin_flag % assign new bin
		layer_match_bin(ii, :) ...
							= id_match_curr([1 1]);
		layer_match_bin((layer_match_simple == layer_match_simple(ii, 1)) | (layer_match_simple == layer_match_simple(ii, 2))) ...
							= id_match_curr;
		bin_id_curr			= bin_id_curr + 1; % increment bin for next time
	else % reassign all layers in higher bin to lower bin
		layer_match_bin((layer_match_bin == id_match_curr(2)) | (layer_match_simple == layer_match_simple(ii, 1)) | (layer_match_simple == layer_match_simple(ii, 2))) ...
							= id_match_curr(1);
	end
end

% num_overlap					= length(find(~layer_match_test));

disp('Sorting bins...')

num_bin_total				= max(layer_match_bin(:));

layer_bin					= cell(1, num_bin_total);
for ii = 1:num_bin_total
	layer_bin{ii}			= unique([layer_match_list((layer_match_bin(:, 1) == ii), 1:2); layer_match_list((layer_match_bin(:, 2) == ii), 3:4)], 'rows');
end

% bin size
[num_bin, ~]				= cellfun(@size, layer_bin);

% bin descending order
[~, ind_bin_ord]            = sort(num_bin, 'descend');

% trim bins
layer_bin                   = layer_bin(ind_bin_ord); % reorder bins
layer_bin                   = layer_bin(logical(num_bin(ind_bin_ord)));
num_bin                     = num_bin(ind_bin_ord);
num_bin                     = num_bin(logical(num_bin));
num_bin_total				= length(num_bin);

disp('Evaluating each bin for overlapping layers...')

% cell of logical vectors with overlap evaluation for each transect (match?)
layer_overlap               = cell(1, num_bin_total);

% evaluate all bins
for ii = 1:num_bin_total
    
    % assume no layer overlap within current bin
    layer_overlap{ii}       = zeros(num_bin(ii), 1);
    
    % unique segments contained in current bin
    seg_unique				= unique(layer_bin{ii}(:, 1), 'stable');
	mm						= 1;

    % evaluate all segments within bin
    for jj = 1:size(seg_unique, 1)
        
        % indices within layer_bin that from current segment
        ind_layer_risk        = find(layer_bin{ii}(:, 1) == seg_unique(jj, :));
        
		% skip if scalar
		if isscalar(ind_layer_risk)
			continue
		end
		
		% loop through layers in current segment, searching for layers that overlap horizontally (non-zero depth difference only)
        ind_bin_bad         = false(length(ind_layer_risk), 1);
		
        for kk = 1:(length(ind_layer_risk) - 1)
			ll_curr			= (kk + 1):length(ind_layer_risk);
			diff_depth_test = find(sum((abs(depth_pk{seg_unique(jj)}(layer_bin{ii}(ind_layer_risk(kk .* ones(length(ll_curr), 1)), 2), :) - depth_pk{seg_unique(jj)}(layer_bin{ii}(ind_layer_risk(ll_curr), 2), :)) > diff_threshold), 2));
			if ~isempty(diff_depth_test)
				ind_bin_bad([kk ll_curr(diff_depth_test)]) ...
                        = true;
			end
        end
        
        % record problematic overlapping layers from same segment within bin
		if ~isempty(find(ind_bin_bad, 1))
        	layer_overlap{ii}(ind_layer_risk(ind_bin_bad)) ...
                            = mm;
        	mm              = mm + 1;
			if (mm > size(colors_def, 1))
            	mm          = 1;
			end
		end
    end
end

% number of layers with overlaps in each bin
num_overlap                 = zeros(1, length(layer_bin));
for ii = 1:length(layer_bin)
    num_overlap(ii)         = length(find(layer_overlap{ii}));
end

beep

%% work through each bin where find(num_overlap)>0

curr_bin					= 1;

% get _curr values to shorten a bit
layer_bin_curr				= layer_bin{curr_bin};
num_bin_curr				= num_bin(curr_bin);
layer_overlap_curr			= layer_overlap{curr_bin};

disp(['Current bin (#' num2str(curr_bin) ') has ' num2str(num_bin_curr) ' layers'])

% simplify matching set to only those relevant for layers in curr_bin
% layer_match_list_curr		= NaN(num_layer_match, 2);
ind_layer_unique			= NaN(num_bin_curr, 1);
for ii = 1:num_bin_curr % for each layer in current bin, simplify nomenclature from segment/layer to simplify unique layer id
	ind_layer_unique(ii)	= find(ismember(layer_unique, layer_bin_curr(ii, :), 'rows'));
	% layer_match_list_curr(ismember(layer_unique, layer_bin_curr(ii, :), 'rows'), 1) ...
    %                         = ii;
	% layer_match_list_curr(ismember(layer_match_list(:, 3:4), layer_bin_curr(ii, :), 'rows'), 2) ...
    %                         = ii;
end
layer_match_simple_curr		= layer_match_simple;
layer_match_simple_curr(~ismember(layer_match_simple_curr, ind_layer_unique)) ...
							= NaN;
layer_match_simple_curr2	= layer_match_simple_curr;
for ii = 1:num_bin_curr
	layer_match_simple_curr2(layer_match_simple_curr == ind_layer_unique(ii)) ...
							= ii;
end
layer_match_simple_curr		= layer_match_simple_curr2;
clear layer_match_simple_curr2

ind_layer_match_curr		= find(sum(~isnan(layer_match_simple_curr), 2)); % indices in layer_match_list from whole
layer_match_simple_curr		= layer_match_simple_curr(ind_layer_match_curr, :); % trim unused rows
num_layer_match_curr		= length(ind_layer_match_curr); % number of current matches
depth_diff_match_curr		= depth_diff_match(ind_layer_match_curr);

% sort rows in simple space
[layer_match_simple_curr, ind_sort] ...
							= sortrows(layer_match_simple_curr, [1 2]);
depth_diff_match_curr		= depth_diff_match_curr(ind_sort);
ind_layer_match_curr		= ind_layer_match_curr(ind_sort);

layer_match_wt				= NaN(1, num_layer_match_curr);
for ii = 1:num_layer_match_curr
	layer_match_wt(ii)		= mean([length(find(layer_match_simple_curr(ii, 1) == layer_match_simple_curr))  length(find(layer_match_simple_curr(ii, 2) == layer_match_simple_curr))]); 
end

% sparse matrix for current bin
layer_match_net             = sparse(layer_match_simple_curr(:, 1), layer_match_simple_curr(:, 2), layer_match_wt, num_bin_curr, num_bin_curr);

% depth_mean for current bin's layers
depth_mean_curr				= zeros(num_bin_curr, 1);
for ii = 1:num_bin_curr
	depth_mean_curr(ii)		= depth_pk_mean{layer_bin_curr(ii, 1)}(layer_bin_curr(ii, 2));
end

% ranges for display
depth_mean_range			= [linspace(min(depth_mean_curr), max(depth_mean_curr), (1 + num_color_def))];
depth_diff_range			= linspace(0, max(depth_diff_match_curr), (1 + num_color_def));

% assign node labels, descriptions, sizes
[node_label, node_desc]		= deal(cell(num_bin_curr, 1));
node_size					= NaN(num_bin_curr, 1);
for ii = 1:num_bin_curr
    node_label{ii}			= num2str(layer_bin_curr(ii, :));
    node_desc{ii}			= [num2str(layer_bin_curr(ii, :)) ' / ' file_pk{layer_bin_curr(ii, 1)} ' #' num2str(layer_bin_curr(ii, 2))];
    node_size(ii)			= max([10 (10 .* sqrt(length(find(layer_match_net(ii, :)))))]);
end

% colorize nodes by mean layer depth
node_color_depth			= colors_def(discretize(depth_mean_curr, depth_mean_range), :);

seg_unique_curr				= unique(layer_bin_curr(:, 1));
colors_curr					= repmat(colors_def, ceil(length(seg_unique_curr) / num_color_def), 1);
colors_curr					= colors_curr(1:length(seg_unique_curr), :);
node_color_seg				= NaN(num_bin_curr, 3);
for ii = 1:length(seg_unique_curr)
	node_color_seg(ismember(layer_bin_curr(:, 1), seg_unique_curr(ii), 'rows'), :) ...
							= repmat(colors_curr(ii, :), length(find(ismember(layer_bin_curr(:, 1), seg_unique_curr(ii), 'rows'))), 1);
end

% gray out non-overlapping layers
grays						= flipud(gray(num_color_def));
grays(1, :)					= [0.93 0.93 0.93]; % don't let lightest edge color be white
node_color_depth(~logical(layer_overlap_curr), :) ...
							= repmat([0.5 0.5 0.5], length(find(~logical(layer_overlap_curr))), 1);
node_color_seg(~logical(layer_overlap_curr), :) ...
							= repmat([0.5 0.5 0.5], length(find(~logical(layer_overlap_curr))), 1);

% generate edge labels for easy identification of matches to remove
edge_label					= cell(num_layer_match_curr, 1);
for ii = 1:num_layer_match_curr
% 	edge_label{ii}			= '';
% end
% for ii = find(depth_diff_match_curr > 100)'
	edge_label{ii}			= num2str(layer_match_simple_curr(ii, :));
end

layer_match_graph			= graph(layer_match_net, 'upper');

set(0, 'DefaultFigureWindowStyle', 'docked')
figure;
pg							= plot(layer_match_graph);
layout(pg, 'force', 'UseGravity', true)
% pg.NodeLabel				= 1:num_bin_curr;
pg.NodeLabelColor			= 'k';
pg.NodeLabelMode			= 'manual';
pg.NodeLabel				= '';
% pg.NodeColor				= node_color_depth;
pg.NodeColor				= node_color_seg;
pg.NodeFontSize				= 14;
pg.EdgeColor				= grays(discretize(depth_diff_match_curr, depth_diff_range), :);
pg.LineWidth				= 2 .* discretize(depth_diff_match_curr, depth_diff_range);
pg.EdgeLabel				= edge_label;
pg.EdgeFontWeight			= 'bold';
pg.EdgeFontSize				= 10;
% pg.EdgeLabelColor			= 'm';
pg.MarkerSize				= node_size;
title(['bin #' num2str(curr_bin) ', ' num2str(num_overlap(curr_bin)) '/' num2str(num_bin_curr) ' layers overlap, ' sprintf('%4.0f', depth_mean_range(1)) '-' sprintf('%4.0f', depth_mean_range(end)) ' m depth range'], 'FontWeight', 'bold', 'FontSize', 16)

%% nodes with overlap and their number of matches (and the actual matches) outputted

for ii = find(layer_overlap_curr)'
	disp([num2str(ii) ': ' num2str(length(find(layer_match_simple_curr == ii))) ' : ' num2str(neighbors(layer_match_graph, ii)')])
end

%% list matches of specific nodes and format in easy way to put into i2u

for ii = []
	a = reshape(sort([neighbors(layer_match_graph, ii)'; ii(ones(1, length(find(layer_match_simple_curr == ii))))]), 1, (2 * length(find(layer_match_simple_curr == ii))));
	as = '';
	for ii = 1:length(a)
		as = [as ' ' num2str(a(ii))];
		if ~mod(ii, 2)
			as = [as '; '];
		end
	end
	disp(as)
end

%% match pairs to potentially remove, test by whiting out

i2u							= [];
highlight(pg, i2u(:, 1), i2u(:, 2), 'EdgeColor', 'w', 'EdgeLabelColor', 'w')

%% remove match pairs that cause overlap

ind2unmatch					= ind_layer_match_curr(ismember(layer_match_simple_curr, sortrows(i2u, [1 2]), 'rows') | ismember(layer_match_simple_curr, sortrows(i2u, [2 1]), 'rows'));
layer_match_list(ind2unmatch, :) ...
							= [];
layer_match_simple(ind2unmatch, :) ...
							= [];
depth_diff_match(ind2unmatch) ...
							= [];
int_match(ind2unmatch)		= [];
save([dir_mat 'layer_match_list.mat'], '-append', 'layer_match_list', 'depth_diff_match', 'int_match')
disp(i2u)
disp('DONE saving')

%%

if do_save
	id_match				= cell(1, num_file_pk);	
	for ii = 1:num_file_pk
		id_match{ii}		= NaN(num_layer(ii), 1);
	end
	% redo id_match based on bin totals
	for ii = 1:length(layer_bin)
    	for jj = 1:num_bin(ii)
        	id_match{layer_bin{ii}(jj, 1)}(layer_bin{ii}(jj, 2)) ...
							= ii;
    	end
	end
	save([dir_mat 'layer_bin_clean.mat'], '-v7.3', 'id_match', 'layer_bin', 'num_bin')
end

%% reference bin deletes

curr_bin = 9;
i2u = [183 78; 221 97; 240 104; 306 133; 397 232];

layer_bin{curr_bin}(ismember(layer_bin{curr_bin}, i2u, 'rows'), :) = [];

curr_bin = 3;
i2u = [347 75];

layer_bin{curr_bin}(ismember(layer_bin{curr_bin}, i2u, 'rows'), :) = [];

[num_bin, ~]				= cellfun(@size, layer_bin);


%% 2d plot radargrams for current bin and overlapping layers

close all

%%

curr_bin					= 1;

disp(['BIN ' num2str(curr_bin)])

ind_file_pk_unique			= unique(layer_bin{curr_bin}(:, 1));
num_file_pk_unique			= length(ind_file_pk_unique);
[amp_elev, dist_lin_radar, elev_radar] ...
							= deal(cell(1, num_file_pk_unique));
ind_decim					= NaN(1, num_file_pk_unique);
for ii = 1:num_file_pk_unique
	disp([num2str(ii) '/' num2str(num_file_pk_unique)])
	[amp_elev{ii}, dist_lin_radar{ii}, elev_radar{ii}, ind_decim(ii)] ...
							= radar_cat_to_elev(path_cat, ['Data_' file_pk{ind_file_pk_unique(ii)}], BM5, permitt_ice, speed_ice, ind_bed_cutoff, dist_smooth, num_trace_blank_min, num_color, num_std_color, dist_decim, sample_decim);
	% layer_curr				= layer_bin{curr_bin}((ind_file_pk_unique(ii) == layer_bin{curr_bin}(:, 1)), 2);	
	% ind_keep				= find(sum(~isnan(depth_pk{ind_file_pk_unique(ii)}(layer_curr, :))), 1):find(sum(~isnan(depth_pk{ind_file_pk_unique(ii)}(layer_curr, :))), 1, 'last');
	% [amp_elev{ii}, dist_lin_radar{ii}] ...
	% 						= deal(amp_elev{ii}(:, ind_keep), dist_lin_radar{ii}(ind_keep));
end

num_panel_max				= 6;
num_fig						= ceil(length(ind_file_pk_unique) / num_panel_max);
pl							= cell(1, num_fig);
for ii = 1:num_fig
	ind_file_curr			= ((ii - 1) * num_panel_max) + (1:num_panel_max);
	ind_file_curr(ind_file_curr > length(ind_file_pk_unique)) ...
							= [];
	pl{ii}					= gobjects(1, size(layer_bin{curr_bin}, 1));
	ll						= 1;
	figure('Color', 'w', 'Alphamap', [0 1])
	colormap(bone(num_color))
	tiledlayout(min([length(ind_file_curr) 3]), min([length(ind_file_curr) 2]), 'TileSpacing', 'tight')
	for jj = ind_file_curr
		nexttile
		hold on
		imagesc(dist_lin_radar{jj}, elev_radar{jj}, amp_elev{jj}, 'CDataMapping', 'direct', 'AlphaData', ~isnan(amp_elev{jj}))
		clim([1 num_color])
		axis xy
		[ind_decim_curr1, ind_decim_curr2, ~] ...
							= intersect(1:ind_decim(jj):num_trace(ind_file_pk_unique(jj)), ind_trace_layer{ind_file_pk_unique(jj)});
		if ~isempty(ind_decim_curr2)
			ind_layer_curr  = layer_bin{curr_bin}((layer_bin{curr_bin}(:, 1) == ind_file_pk_unique(jj)), 2);
			elev_pk_tmp		= elev_pk{ind_file_pk_unique(jj)}(ind_layer_curr, :);
			if (num_trace(ind_file_pk_unique(jj)) < num_trace(ind_file_pk_unique(jj)))
				elev_pk_curr= NaN(size(elev_pk_tmp, 1), num_trace(ind_file_pk_unique(jj)));
				elev_pk_curr(:, ind_trace_layer{ind_file_pk_unique(jj)}) ...
							= elev_pk_tmp;
			else
				elev_pk_curr= elev_pk_tmp;
			end
			elev_pk_curr	= elev_pk_curr(:, ind_decim_curr1);
			clear elev_pk_tmp
			[xlim_curr, ylim_curr] ...
							= deal(NaN(1, 2));
			for kk = 1:size(elev_pk_curr, 1)
				pl{ii}(ll)	= plot(dist_lin_radar{jj}(ind_decim_curr2), elev_pk_curr(kk, :), 'LineWidth', 2, 'Color', 'g', 'Tag', num2str(ind_layer_curr(kk)));
				pl{ii}(ll).DataTipTemplate.DataTipRows(end + 1) = dataTipTextRow(num2str(ind_layer_curr(kk)), 'Tag');
				xlim_curr	= [min([xlim_curr(1) dist_lin_radar{jj}(ind_decim_curr2(find(~isnan(elev_pk_curr(kk, :)), 1)))]) max([xlim_curr(2) dist_lin_radar{jj}(ind_decim_curr2(find(~isnan(elev_pk_curr(kk, :)), 1, 'last')))])];
				ylim_curr	= [min([ylim_curr(1) min(elev_pk_curr(kk, :))]) max([ylim_curr(2) max(elev_pk_curr(kk, :))])];
				ll			= ll + 1;
			end
			xlim(xlim_curr)
			ylim(ylim_curr + [-100 100])
		end
		set(gca, 'FontSize', 20, 'FontWeight', 'bold')
		title([num2str(ind_file_pk_unique(jj)) ': ' num2str(layer_bin{curr_bin}((layer_bin{curr_bin}(:, 1) == ind_file_pk_unique(jj)), 2)')])
		box on
	end
	pl{ii}					= pl{ii}(ishandle(pl{ii}));
	set(gcf, 'KeyPressFcn', {@keypress, [pl{ii}]})
end

figure(1)

beep

%%

function keypress(~, event, pl)
	if strcmp(event.Key, 'f')
		if strcmp(pl(1).Visible, 'on')
			[pl(:).Visible] = deal('off');
		else
			[pl(:).Visible] = deal('on');
		end
	end
end