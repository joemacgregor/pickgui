function [amp_elev, dist_lin, elev, ind_decim] ...
							= radar_cat_to_elev(path_cat, file_cat, bm_grd, permitt_ice, speed_ice, ind_bed_cutoff, dist_smooth, num_trace_blank_min, num_color, num_std_color, dist_decim, sample_decim)
% RADAR_CAT_TO_ELEV.M
% 
% Generate concatenated radargram elevation corrected for display.
% 
% Joe MacGregor (NASA)
% Last updated: 7 May 2024

data_cat					= load([path_cat file_cat], 'amp', 'dist', 'dist_lin', 'elev_air', 'elev_surf', 'num_sample', 'num_trace', 'twtt', 'twtt_bed', 'twtt_surf', 'x', 'y');

% check for unique locations and if not then trim
if (size(unique([data_cat.x' data_cat.y'], 'rows'), 1) < data_cat.num_trace)
	[~, ind_unique]			= unique([data_cat.x' data_cat.y'], 'rows');
	data_cat.amp(:, setdiff(1:data_cat.num_trace, ind_unique)) ...
							= NaN;
end

% along-track BMv5 surface elevation (geoid)
elev_surf_grd_curr			= interp2(bm_grd.x, bm_grd.y, bm_grd.elev_surf, data_cat.x, data_cat.y);

ind_surf_curr				= interp1(data_cat.twtt, 1:data_cat.num_sample, data_cat.twtt_surf, 'nearest', 'extrap'); % surface traveltime indices
if ~isempty(find(isnan(ind_surf_curr), 1))
    ind_surf_curr(isnan(ind_surf_curr)) ...
							= round(interp1(find(~isnan(ind_surf_curr)), ind_surf_curr(~isnan(ind_surf_curr)), find(isnan(ind_surf_curr)), 'linear', 'extrap'));
end
ind_surf_curr(ind_surf_curr < 1) ...
							= 1;
ind_surf_curr(ind_surf_curr > data_cat.num_sample) ...
							= data_cat.num_sample;

% shift data up to surface
amp_depth					= NaN(data_cat.num_sample, data_cat.num_trace, 'single');
ind_bottom_keep				= 1 + data_cat.num_sample - ind_surf_curr;
for jj = 1:data_cat.num_trace
	amp_depth(1:ind_bottom_keep(jj), jj) ...
							= data_cat.amp(ind_surf_curr(jj):data_cat.num_sample, jj);
end

dt							= median(diff(data_cat.twtt)); % traveltime sampling interval, s
depth						= (speed_ice / 2) .* (0:dt:((data_cat.num_sample - 1) * dt))'; % simple monotonically increasing depth vector, m

% do some basic amplitude corrections to amp_depth to get scaling ok, best to do before topographic correction
amp_depth					= amp_depth + (3e-2 .* depth(:, ones(1, data_cat.num_trace))); % correct for two-way attenuation (assume 15 dB/km one-way)
amp_depth					= amp_depth + repmat((2e1 .* log10(data_cat.elev_air - data_cat.elev_surf)), data_cat.num_sample, 1); % correct for geometric spreading to ice surface
amp_depth(2:end, :)			= amp_depth(2:end, :) + (2e1 .* log10(depth(2:end, ones(1, data_cat.num_trace)) ./ sqrt(permitt_ice))); % correct for geometric spreading beneath ice surface, skip depth=0
if ~isreal(amp_depth)
	amp_depth				= real(amp_depth);
end
amp_depth(isinf(amp_depth))	= NaN; % should no longer be an issue but just in case

% demean using a moving median mean
amp_depth					= amp_depth - smoothdata(amp_depth, 2, 'movmedian', round(dist_smooth / median(diff(data_cat.dist), 'omitnan')));

% blank out below bed
ind_bed						= interp1(data_cat.twtt, 1:data_cat.num_sample, data_cat.twtt_bed, 'nearest', 'extrap'); % bed traveltime indices
if ~isempty(find(isnan(ind_bed), 1)) % handle missing bed picks using BedMachine
    ind_bed(isnan(ind_bed)) = ind_surf_curr(isnan(ind_bed)) + round((2 / (speed_ice * dt)) .* interp2(bm_grd.x, bm_grd.y, bm_grd.thick, data_cat.x(isnan(ind_bed)), data_cat.y(isnan(ind_bed)), 'nearest', 0));
end
ind_bed_trim				= ind_bed - ind_surf_curr + ind_bed_cutoff;
ind_bed_trim(ind_bed_trim < 1) ...
							= 1;
ind_bed_trim(ind_bed_trim > data_cat.num_sample) ...
							= data_cat.num_sample;
for jj = 1:data_cat.num_trace
	amp_depth(min([data_cat.num_sample ind_bed_trim(jj)]):end, jj) ...
							= NaN;
end

% trim off NaN rows if present
if ~isempty(find((sum(~isnan(amp_depth), 2) < num_trace_blank_min), 1))
	depth					= depth(sum(~isnan(amp_depth), 2) >= num_trace_blank_min);
	amp_depth				= amp_depth((sum(~isnan(amp_depth), 2) >= num_trace_blank_min), :);
end

% normalize amplitudes for simpler display
range_curr					= linspace((mean(amp_depth, 'all', 'omitnan') - (num_std_color * std(amp_depth, 0, 'all', 'omitnan'))), (mean(amp_depth, 'all', 'omitnan') + (num_std_color * std(amp_depth, 0, 'all', 'omitnan'))), num_color);
amp_depth					= discretize(amp_depth, [-Inf range_curr(2:(end - 1)) Inf]);
elev						= flipud(max(elev_surf_grd_curr) - depth); % elevation vector, m

% topographically correct data and flip
amp_elev					= single(flipud(topocorr(amp_depth, depth, elev_surf_grd_curr)));

ind_decim					= round(dist_decim / median(diff(data_cat.dist)));
ind_decim_vec				= 1:ind_decim:data_cat.num_trace;
[amp_elev, dist_lin, elev]	= deal(amp_elev(1:sample_decim:end, ind_decim_vec), data_cat.dist_lin(ind_decim_vec), elev(1:sample_decim:end));