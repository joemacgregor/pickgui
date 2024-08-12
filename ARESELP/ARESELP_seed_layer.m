function [layer, ind_layer, stop_flag] ...
							= ARESELP_seed_layer(pk_cwt, num_sample, num_trace, seed_pt, vert_sep_min, len_blk, vert_jump_max, angle_diff, vert_proxim_max, len_layer_min)
% ARESELP_SEED_LAYER Loop through layer seed points, trace individual layers and do some post-processing.
% 
%	Only use ARESELP_SEED_LAYER as called by ARESELP.
% 
%	ARESELP_SEED_LAYER is based on code archived at:
%	https://github.com/xiongsiting/ARESELP/
% 
% Joe MacGregor (NASA)
% Last updated: 12 August 2024

sz_data						= [num_sample num_trace]; % matrix size
layer						= NaN(sz_data); % matrix with layer numbers
layer_tmp					= false(sz_data); % temporary matrix for near-layer trimming

disp(['Evaluating ' num2str(size(seed_pt, 1)) ' seed points...']);

[ii, ind_layer]				= deal(1); % start with first seed pt for first layer
stop_flag					= {}; % record reasons for layers stopping

while true % keep going until no more non-NaN seed points
	
	[layer, ind_layer, stop_flag{end + 1}] ...
							= ARESELP_trace_layer(pk_cwt, num_sample, num_trace, seed_pt(ii, :), layer, ind_layer, vert_sep_min, len_blk, vert_jump_max, angle_diff, vert_proxim_max, len_layer_min); %#ok<AGROW> 
	
	% NaN out just-used seed point
	seed_pt(ii, :)			= NaN;
	
	% unsuccessful layer so move to next seed pt
	if isempty(find((layer == ind_layer), 1))
		ii					= find(~isnan(seed_pt(:, 1)), 1);
		if isempty(ii)
			ind_layer		= ind_layer - 1; % decrement ind_layer to account for failed last layer or end of loop
			break
		else
			continue
		end
	end
	
	% find seed pts that are close to current layer then get rid of them
	layer_tmp(layer == ind_layer) ...
							= true;
	[ind_close_i, ind_close_j] ...
							= find(bwdist(layer_tmp) < vert_sep_min);
	seed_pt(ismember(seed_pt, [ind_close_i ind_close_j], 'rows'), :) ...
							= NaN;
	layer_tmp(layer_tmp)	= false; % prep for next layer
	
	% done filtering, now move to next seed pt
	ii						= find(~isnan(seed_pt(:, 1)), 1);
	if isempty(ii)
		break
	end
	
	% increment layer number
	ind_layer				= max(layer, [], 'all') + 1;

	disp([num2str(length(find(~isnan(seed_pt(:, 1))))) ' left...'])
end