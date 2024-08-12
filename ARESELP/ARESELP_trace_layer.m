function [layer, ind_layer, stop_flag] ...
							= ARESELP_trace_layer(pk_cwt, num_sample, num_trace, seed_pt, layer, ind_layer, vert_sep_min, len_blk, vert_jump_max, angle_diff, vert_proxim_max, len_layer_min)
% ARESELP_TRACE_LAYER Trace an individual layer using Hough transforms.
% 
%	Only use ARESELP_TRACE_LAYER as called by ARESELP_SEED_LAYER.
% 
%	ARESELP_TRACE_LAYER is based on code archived at:
%	https://github.com/xiongsiting/ARESELP/
% 
% Joe MacGregor (NASA)
% Last updated: 12 August 2024

sz_data						= [num_sample num_trace]; % data matrix size

% houghlines gap filling parameter
gap_fill					= ceil(0.25 * len_blk);

% half-width test block size
half_blk					= floor(len_blk / 2);

% reference/starting Hough-transform-derived slopes
[hough_slope_old, hough_slope_first] ...
							= deal(0);

% start with seed pt
[ii, jj]					= deal(seed_pt(1), seed_pt(2));

% propagation direction from seed pt: right first then left
dir_prop					= 'right';

% default flag if everything goes well
stop_flag					= 'layer complete';

while true % first propagate layer to right, then left (see end of while loop)

	i_blk					= max([1 (ii - vert_sep_min)]):min([(ii + vert_sep_min) num_sample]); % vertical/i indices of block in data matrix
	j_blk					= max([1 (jj - half_blk)]):min([(jj + half_blk) num_trace]); % horizontal j indices of block in data matrix
	
	% identify most significant line from Hough transform
	[H, theta, rho]			= hough(pk_cwt(i_blk, j_blk));
	hough_lines_pk			= houghlines(pk_cwt(i_blk, j_blk), theta, rho, houghpeaks(H, 3), 'FillGap', gap_fill, 'MinLength', ceil(0.25 * length(j_blk)));
	if isempty(hough_lines_pk)
		stop_flag			= [dir_prop ': no lines'];
		break
	end
	
	if (length(hough_lines_pk) > 1) % evaluate slopes of up to top 3 lines
		
		[hough_slope, i_test] ...
							= deal(NaN(1, length(hough_lines_pk)));
		for kk = 1:length(hough_lines_pk)
			hough_slope(kk)	= (hough_lines_pk(kk).point2(2) - hough_lines_pk(kk).point1(2)) / (hough_lines_pk(kk).point2(1) - hough_lines_pk(kk).point1(1)); 
		end
		
		hough_slope(isinf(hough_slope)) ...
							= NaN; % NaN out bad slopes
		
		if isempty(find(~isnan(hough_slope), 1)) % stop if no good slopes
			stop_flag		= [dir_prop ': no good slopes'];
			break
		end
		
		for kk = find(~isnan(hough_slope)) % loop through good slopes to get i value at jj
			i_test(kk)		= i_blk(1) + round(interp1(((j_blk(1) - 1) + [hough_lines_pk(kk).point1(1) hough_lines_pk(kk).point2(1)]), [hough_lines_pk(kk).point1(2) hough_lines_pk(kk).point2(2)], jj)) - 1;
		end
		
		if all(abs(i_test(~isnan(i_test)) - ii) <= vert_jump_max) % all remaining slopes passed the test
			hough_slope		= mean(hough_slope, 'omitnan'); % average slope
		elseif isscalar(find(abs(i_test(~isnan(i_test)) - ii) <= vert_jump_max)) % only one good slope
			hough_slope		= hough_slope(abs(i_test(~isnan(i_test)) - ii) <= vert_jump_max);
		else % potentially 1-2 good slopes
			[diff_i_best, ind_slope_best] ...
							= min(abs(i_test - ii)); % pick slope that changes vertical offset at [ii,jj] the least
			if (diff_i_best > vert_jump_max) % stop if vertical jump at jj too high
				stop_flag	= [dir_prop ': slope with smallest jump still too big'];
				break
			end
			hough_slope		= hough_slope(ind_slope_best); % keep best slope
		end
		
	else % only one line, easy
		hough_slope			= (hough_lines_pk(1).point2(2) - hough_lines_pk(1).point1(2)) / (hough_lines_pk(1).point2(1) - hough_lines_pk(1).point1(1));
	end
	
	% stop if slope changed too much
	if (isnan(hough_slope) || ((abs(atand(hough_slope) - atand(hough_slope_old))) > angle_diff))
		stop_flag			= [dir_prop ': slope difference too big'];
		break
	end
	
	i_blk_line				= ii + round(hough_slope .* (j_blk - jj)); % vertical/i indices of line
	
	% stop and merge if crossing another layer
	ind_layer_old			= layer(sub2ind(sz_data, i_blk_line, j_blk)); % layer numbers in current line
	
	if ~isempty(find((~isnan(ind_layer_old) & (ind_layer_old ~= ind_layer)), 1))
		switch dir_prop
			case 'right'
				j_cross		= find((~isnan(ind_layer_old) & (ind_layer_old ~= ind_layer)), 1); % first out going right
				ind_layer_cross ...
							= ind_layer_old(j_cross); % layer 
				layer(sub2ind(sz_data, i_blk_line(1:j_cross), j_blk(1:j_cross)))...
							= ind_layer_cross; % assign NaNs in new line prior to crossing to old crossing layer
			case 'left'
				j_cross		= find((~isnan(ind_layer_old) & (ind_layer_old ~= ind_layer)), 1, 'last'); % first out going left
				ind_layer_cross ...
							= ind_layer_old(j_cross); % last out going left
				layer(sub2ind(sz_data, i_blk_line(j_cross:end), j_blk(j_cross:end)))...
							= ind_layer_cross; % assign NaNs in new line prior to crossing to old crossing layer
		end
		layer(layer == ind_layer) ...
							= ind_layer_cross; % change labeling of new layer to old crossing layer
		ind_layer			= ind_layer_cross; % change layer label to old layer
		stop_flag			= [dir_prop ': new layer crosses old layer so merged'];
		break
	end
	
	% traces of current block that have other layers in them
	[~, j_blk_test]			= find(~isnan(layer(:, j_blk)) & (layer(:, j_blk) ~= ind_layer));
	
	% stop if layer is getting too close to another in same trace
	if ~isempty(j_blk_test)
		j_blk_test			= unique(j_blk_test); % only need unique traces
		for kk = 1:length(j_blk_test)
			i_test			= find(~isnan(layer(:, j_blk(j_blk_test(kk)))) & (layer(:, j_blk(j_blk_test(kk))) ~= ind_layer));
			[min_val, min_ind] ...
							= min(abs(i_test - i_blk_line(j_blk_test(kk))));
			if (min_val < vert_proxim_max) % if layer gets too close, stop layer
				ind_layer_close ...
							= layer(i_test(min_ind), j_blk(j_blk_test(kk)));
				switch dir_prop
					case 'right'
						layer(sub2ind(sz_data, i_blk_line(1:j_blk_test(kk)), j_blk(1:j_blk_test(kk))))...
							= ind_layer_close; % assign NaNs in new line prior to crossing to old nearing layer
					case 'left'
						layer(sub2ind(sz_data, i_blk_line(j_blk_test(kk):end), j_blk(j_blk_test(kk):end)))...
							= ind_layer_close; % assign NaNs in new line prior to crossing to old nearing layer
				end
				layer(layer == ind_layer) ...
							= ind_layer_close;  % change labeling of new layer to old crossing layer
				ind_layer	= ind_layer_close;  % change layer label to old layer
				stop_flag	= [dir_prop ': new layer nearing old layer so merging'];
				break
			end
		end
	end
	
	% survived tests, so assign this slope/line to a layer
	layer(sub2ind(sz_data, i_blk_line, j_blk)) ...
							= ind_layer;
	
	% preserve previous slope
	hough_slope_old			= hough_slope;
	if ((jj == seed_pt(2)) && strcmp(dir_prop, 'right'))
		hough_slope_first	= hough_slope_old;
	end
	
	% conclude this instance of the loop based on direction
	switch dir_prop
		case 'right'
			if (j_blk(end) == num_trace) % stop if reached right end of radargram, go back to seed point and propagate left
				dir_prop	= 'left';
				[ii, jj]	= deal(seed_pt(1), seed_pt(2));
				hough_slope_old	= hough_slope_first;
			else
				[ii, jj]	= deal(i_blk_line(end), j_blk(end)); % move layer test block to end of slope
			end
		case 'left'
			if (j_blk(1) == 1) % stop if reached start (left end) of radargram
				break
			else
				% move layer test block to beginning of slope
				[ii, jj]	= deal(i_blk_line(1), j_blk(1));
			end
	end
end

% horizontal indices where current layer is
[~, ind_layer_j]			= find(layer == ind_layer);
ind_layer_j					= unique(ind_layer_j);

% delete layer if not long enough
if (length(ind_layer_j) < len_layer_min)
	layer(layer == ind_layer) ...
							= NaN;
	if ~strcmp(stop_flag, 'layer complete')
		stop_flag			= [stop_flag '; new layer too short'];
	end
	return
end

% resolve vertical repeats of current layer
for jj = ind_layer_j(sum(layer(:, ind_layer_j) == ind_layer) > 1)'
	layer((layer(:, jj) == ind_layer), jj) ...
							= Inf; % Inf out current layer in current column
	layer(round(mean(find(isinf(layer(:, jj))))), jj) ...
							= ind_layer; % average out layer position in current column
	layer(isinf(layer(:, jj)), jj) ...
							= NaN; % return to NaNs
end