% GRIS_STRAT2_MOV Generate movie that shows all traced radargrams.
% 
% Joe MacGregor (NASA)
% Last updated: 13 November 2024

clear

plotting                    = false;

dir_mat						= '/Users/jamacgre/OneDrive - NASA/research/matlab/pickgui_v2/mat/';
dir_data					= '/Users/jamacgre/OneDrive - NASA/data/gris_strat2/arctic_cat_r1/';
path_mov_out				= '/Users/jamacgre/OneDrive - NASA/research/manuscripts/gris_strat_v2/fig/';
file_mov_out				= 'mov1';

load([dir_mat 'xyz_all.mat'], 'campaign', 'x', 'y')
pk_cat						= load([dir_mat 'pk_cat.mat'], 'depth', 'file_pk', 'ind_campaign', 'num_file_pk', 'num_layer', 'ind_trace_layer', 'num_trace', 'thick', 'x', 'y');
load([dir_mat 'core_int.mat'], 'name_core', 'int_core_cat', 'num_core', 'x_core', 'y_core')
load([dir_mat 'date_all.mat'], 'age')
load([dir_mat 'age_grd2_clean.mat'], 'x_grd', 'y_grd')
[xx_grd, yy_grd]			= meshgrid(x_grd, y_grd);

% BedMachine v5
BM5                         = struct;
BM5.x                       = double(ncread('/Users/jamacgre/OneDrive - NASA/research/data/greenland/BedMachine/BedMachineGreenland-v5.nc', 'x'))'; % projected x, m
BM5.y                       = double(flipud(ncread('/Users/jamacgre/OneDrive - NASA/research/data/greenland/BedMachine/BedMachineGreenland-v5.nc', 'y'))); % projected y, m
BM5.mask_gris               = double(rot90(ncread('/Users/jamacgre/OneDrive - NASA/research/data/greenland/BedMachine/BedMachineGreenland-v5.nc', 'mask'))); % ice mask
BM5.elev_surf				= double(rot90(ncread('/Users/jamacgre/OneDrive - NASA/research/data/greenland/BedMachine/BedMachineGreenland-v5.nc', 'surface'))); % surface elevation, m

% simplified masks for whole Greenland maps
BM5.mask_combo              = BM5.mask_gris;
BM5.mask_gris(BM5.mask_gris ~= 2) ...
                            = 0;
BM5.mask_gris(BM5.mask_gris == 2) ...
                            = 1;
BM5.mask_gris               = logical(BM5.mask_gris); % now mask_gris is for ice sheet only
BM5.mask_combo(BM5.mask_combo == 3) ...
                            = 0;
BM5.mask_combo(BM5.mask_combo == 4) ...
                            = 1; % mask_combo includes 0/ocean 1/land 2/ice
BM5.mask_combo              = uint8(BM5.mask_combo); % reduce size

% now BM5 mask 2=peripheral ice masses and 3=ice sheet (previously 2)
BM5.mask_combo_plot			= BM5.mask_combo;

BM5_mask_color				= cell(1, 3);
color_ocean					= [135 206 235] ./ 255;
color_land					= [159 89 39] ./ 255;
for ii = 1:3
	BM5_mask_color{ii}		= 0.9 .* ones(length(BM5.y(1:50:end)), length(BM5.x(1:50:end)));
	BM5_mask_color{ii}((BM5.mask_combo_plot(1:50:end, 1:50:end) == 0)) ...
							= color_ocean(ii);
	BM5_mask_color{ii}((BM5.mask_combo_plot(1:50:end, 1:50:end) == 1)) ...
							= color_land(ii);
end
BM5_mask_color				= cat(3, BM5_mask_color{1}, BM5_mask_color{2}, BM5_mask_color{3});

num_color_radar				= 2 ^ 8; % number of colors for display
param_anim					= [5 20 50];

speed_vacuum                = 299792458; % m/s
permitt_ice                 = 3.15; % dimensionless
speed_ice                   = speed_vacuum / sqrt(permitt_ice); % m/s

% radar system compilation (30 campaigns, ingoring ground)
radar_system				= {'ICORDS' 'ICORDS' 'ICORDS' 'ICORDS' 'ICORDSv2' 'ICORDSv2' 'ICORDSv2' 'ICORDSv2' 'ACORDS' 'ACORDS' 'MCRDS' 'MCRDS' '' 'MCRDS' 'MCRDS' 'MCoRDS' 'MCoRDS' 'MCoRDSv2' 'MCoRDSv2' 'MCoRDSv2' 'MCoRDSv3' 'MCoRDSv3' 'MCoRDSv5' '' 'MCoRDSv5' 'MCoRDSv5' '' 'MCoRDSv3' 'MCoRDSv3' 'MCoRDSv3'};

age_color					= 1e3 .* [-2 0 5 11.7 20 50 115 150];
color_age					= [1 1 1; parula(length(age_color) - 2)];

% simplify filenames and unsparsify depth (much simpler), and correct for firn
depth_pk					= cell(1, pk_cat.num_file_pk);
for ii = 1:pk_cat.num_file_pk
	depth_tmp				= full(pk_cat.depth{ii});
	depth_tmp(~depth_tmp)	= NaN;
	depth_pk{ii}			= NaN(pk_cat.num_layer(ii), pk_cat.num_trace(ii));
	depth_pk{ii}(:, pk_cat.ind_trace_layer{ii}) ...
							= depth_tmp;
end
pk_cat						= rmfield(pk_cat, 'depth');

%%
if plotting
%%

set(0, 'DefaultFigureWindowStyle', 'default') %#ok<UNRCH>
frad						= figure('Position', [100 100 960 540], 'Color', 'w');
colormap(bone(num_color_radar))
ax_radar					= subplot('Position', [0.08 0.11 0.89 0.82]);
hold on
p_rad						= imagesc([0 param_anim(3)], [0 1], NaN(2), [0 1]);
p_surf						= line(NaN, NaN, 'Marker', '.', 'Color', [0.5 0 0], 'MarkerSize', 6, 'LineStyle', 'none');
p_bed						= line(NaN, NaN, 'Marker', '.', 'Color', 'r', 'MarkerSize', 6, 'LineStyle', 'none');
p_pk						= gobjects(1, max(pk_cat.num_layer));
for ii = 1:max(pk_cat.num_layer)
	p_pk(ii)				= line(NaN, NaN, 'Marker', '.', 'Color', 'w', 'MarkerSize', 6, 'LineStyle', 'none');
end
p_lg						= gobjects(1, length(color_age));
for ii = 1:length(color_age)
	p_lg(ii)				= line(NaN, NaN, 'LineWidth', 3, 'Color', color_age(ii, :));		
end
axis([0 param_anim(3) 0 1])
axis xy
ax_radar.FontSize			= 16;
ax_radar.Layer				= 'top';
xlabel('Distance (km)')
ylabel('Elevation (m)')
ti							= title('', 'FontWeight', 'bold', 'Interpreter', 'none');
box on
annotation('textbox', [0.9677 0.3272 0.0906 0.0546], 'String', 'Age (ka)', 'FontSize', 16, 'Color', 'k', 'EdgeColor', 'none', 'Rotation', 270)
legend(p_lg, {'undated' '0-5' '5-11.7' '11.7-20' '20-50' '50-115' '>115'}, 'Position', [0.8656 0.1104 0.1042 0.3000], 'FontSize', 16, 'Color', [0.8 0.8 0.8], 'AutoUpdate', 'off')
ax_map						= axes('Position', [0.78 0.48 0.25 (0.25 * (diff(BM5.y([1 end])) / diff(BM5.x([1 end]))) * (ax_radar.Position(3) / ax_radar.Position(4)))], 'Color', 'w');
image(BM5.x(1:50:end), BM5.y(1:50:end), BM5_mask_color);
axis xy image
pl							= gobjects(1, pk_cat.num_file_pk);
for ii = 1:pk_cat.num_file_pk
	pl(ii)					= line(pk_cat.x{ii}(1:50:end), pk_cat.y{ii}(1:50:end), 'Color', 'k', 'LineWidth', 0.5);
end
plc							= line(NaN, NaN, 'Color', 'r', 'LineWidth', 4);
p_core_map					= line(x_core, y_core, 'Color', 'k', 'Marker', '^', 'MarkerSize', 8, 'MarkerFaceColor', 'm', 'LineStyle', 'none', 'LineWidth', 1);
[p_core, p_core_label]		= deal(gobjects(0));
axis xy image
[ax_map.XTick, ax_map.YTick]= deal([]);
box on

mp4							= VideoWriter([path_mov_out file_mov_out '.mp4'], 'MPEG-4');
mp4.FrameRate				= param_anim(1);
mp4.Quality					= 100;
open(mp4)
 
axes(ax_radar)
uistack(ax_map, 'top')
uistack(plc, 'top')

tic
for ii = 1:pk_cat.num_file_pk
	
	data_cat				= load([dir_data pk_cat.file_pk{ii}(1:(end - 7)) '.mat'], 'amp', 'dist_lin', 'elev_air', 'elev_surf', 'num_sample', 'num_trace', 'twtt', 'twtt_surf', 'x', 'y');
	data_cat.dist_lin		= 1e-3 .* data_cat.dist_lin; % m to km because we're using it so often below
	
    ind_surf				= NaN(1, data_cat.num_trace);
    ind_surf(~isnan(data_cat.twtt_surf) & ~isinf(data_cat.twtt_surf)) ...
							= interp1(data_cat.twtt, 1:data_cat.num_sample, data_cat.twtt_surf(~isnan(data_cat.twtt_surf) & ~isinf(data_cat.twtt_surf)), 'nearest', 'extrap');
    if ~isempty(find(isnan(ind_surf), 1))
        ind_surf(isnan(ind_surf)) ...
							= round(interp1(find(~isnan(ind_surf)), ind_surf(~isnan(ind_surf)), find(isnan(ind_surf)), 'linear', 'extrap'));
    end
    ind_surf(ind_surf <= 0) = 1;
	
	% get to depth
	amp_curr				= NaN(data_cat.num_sample, data_cat.num_trace, 'single');
	for jj = find(~isnan(ind_surf))
    	amp_curr(1:(data_cat.num_sample - ind_surf(jj) + 1), jj) ...
							= data_cat.amp(ind_surf(jj):data_cat.num_sample, jj); % shift data up to surface
	end

	dt						= median(diff(data_cat.twtt)); % traveltime sampling interval, s
	depth_vec				= (speed_ice / 2) .* (0:dt:((data_cat.num_sample - 1) * dt))'; % simple monotonically increasing depth vector, m
	
	% % do some amplitude corrections
	% amp_curr				= amp_curr + (3e-2 .* depth_vec(:, ones(1, data_cat.num_trace))); % correct for two-way attenuation (assume 15 dB/km one-way)
	% amp_curr				= amp_curr + repmat((2e1 .* log10(data_cat.elev_air - data_cat.elev_surf)), data_cat.num_sample, 1); % correct for geometric spreading to ice surface
	% amp_curr(2:end, :)		= amp_curr(2:end, :) + (2e1 .* log10(depth_vec(2:end, ones(1, data_cat.num_trace)) ./ sqrt(permitt_ice))); % correct for geometric spreading beneath ice surface, skip depth=0
	% if ~isreal(amp_curr)
	% 	amp_curr			= real(amp_curr);
	% end
	% amp_curr(isinf(amp_curr)) ...
	% 						= NaN; % should no longer be an issue but just in case
	
	% topographically correct data and flip
	elev_surf_grd_curr		= interp2(BM5.x, BM5.y, BM5.elev_surf, data_cat.x, data_cat.y);
	elev_bed_curr			= elev_surf_grd_curr - pk_cat.thick{ii};
	elev_bed_curr(isinf(elev_bed_curr)) ...
							= NaN;
	elev_vec				= flipud(max(elev_surf_grd_curr) - depth_vec); % elevation vector, m
	
	% handle cutoff
	if (min(elev_vec) > min(elev_bed_curr))
		amp_curr			= single(flipud(topocorr(amp_curr, depth_vec, elev_surf_grd_curr, false)));
		depth_vec			= (speed_ice / 2) .* (0:dt:((size(amp_curr, 1) - 1) * dt))';
		elev_vec			= flipud(max(elev_surf_grd_curr) - depth_vec);
		ind_good			= find(elev_vec >= min(elev_bed_curr));
		[amp_curr, elev_vec]= deal(amp_curr(ind_good, :), elev_vec(ind_good));
	else
		amp_curr			= single(flipud(topocorr(amp_curr, depth_vec, elev_surf_grd_curr, true)));
	end
	
	% plot radargram
	[p_rad.XData, p_rad.YData, p_rad.CData] ...
							= deal(data_cat.dist_lin, elev_vec, amp_curr);
	[p_surf.XData, p_surf.YData] ...
							= deal(data_cat.dist_lin, elev_surf_grd_curr);
	[p_bed.XData, p_bed.YData] ...
							= deal(data_cat.dist_lin, elev_bed_curr);
	
	age_curr				= age{ii};
	age_curr(isnan(age_curr)) ...
							= -1;
	age_color_curr			= discretize(age_curr, age_color);
	
	for jj = 1:pk_cat.num_layer(ii)
		ind_good			= find(~isnan(depth_pk{ii}(jj, :)));
		[p_pk(jj).XData, p_pk(jj).YData, p_pk(jj).Color] ...
							= deal(data_cat.dist_lin(ind_good), (elev_surf_grd_curr(ind_good) - depth_pk{ii}(jj, ind_good)), color_age(age_color_curr(jj), :));
	end
	if ((ii > 1) && (pk_cat.num_layer(ii - 1) > pk_cat.num_layer(ii))) % clear out previous layers
		for jj = (pk_cat.num_layer(ii) + 1):(pk_cat.num_layer(ii - 1))
			[p_pk(jj).XData, p_pk(jj).YData] ....
							= deal(NaN);
		end
	end
	
	ti.String				= [pk_cat.file_pk{ii}(6:(end - 7)) ' : ' campaign{pk_cat.ind_campaign(ii)} ' : ' radar_system{pk_cat.ind_campaign(ii)} ' : ' num2str(ii) '/' num2str(pk_cat.num_file_pk) '      '];
	
	ax_radar.XLim			= [0 min([param_anim(3) max(data_cat.dist_lin)])];
	ind_dist_curr			= interp1(data_cat.dist_lin, 1:data_cat.num_trace, ax_radar.XLim, 'nearest', 'extrap');	
	if isempty(find(~isnan(elev_bed_curr(ind_dist_curr(1):ind_dist_curr(2))), 1))
		ax_radar.YLim		= [elev_vec(1) max(elev_surf_grd_curr(ind_dist_curr(1):ind_dist_curr(2)))];		
	else
		ax_radar.YLim		= [min(elev_bed_curr(ind_dist_curr(1):ind_dist_curr(2)), [], 'omitnan') max(elev_surf_grd_curr(ind_dist_curr(1):ind_dist_curr(2)))];
	end
	ax_radar.CLim			= mean(amp_curr(1:10:end, ind_dist_curr(1):10:ind_dist_curr(2)), 'all', 'omitnan') + (std(amp_curr(1:10:end, ind_dist_curr(1):10:ind_dist_curr(2)), 0, 'all', 'omitnan') .* [-2 2]);
	
	if (ii > 1)
		pl(ii - 1).Color	= 'k';
		pl(ii - 1).LineWidth= 0.5;
	end
	pl(ii).Color			= 'b';
	pl(ii).LineWidth		= 3;
	uistack(pl(ii), 'top')
	
	[plc.XData, plc.YData]	= deal(pk_cat.x{ii}(ind_dist_curr), pk_cat.y{ii}(ind_dist_curr));
	uistack(plc, 'top')
	
	delete(p_core(isgraphics(p_core)))
	delete(p_core_label(isgraphics(p_core_label)))
	
	% add core labels
	if ~isempty(find((int_core_cat(:, 1) == ii), 1))
		[p_core, p_core_label] ...
							= deal([]);
		for jj = find(int_core_cat(:, 1) == ii)'
			p_core			= [p_core line(data_cat.dist_lin(int_core_cat([jj jj], 3)), [elev_bed_curr(int_core_cat(jj, 3)) elev_surf_grd_curr(int_core_cat(jj, 3))], 'LineWidth', 3, 'LineStyle', '--', 'Color', 'm')];
			p_core_label	= [p_core_label text((data_cat.dist_lin(int_core_cat(jj, 3)) + 1), (elev_surf_grd_curr(int_core_cat(jj, 3)) - 100), name_core(int_core_cat(jj, 2)), 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'm', 'Rotation', 270)];
		end
	end
	
	rect_curr				= frad.Position(3:4);
	writeVideo(mp4, getframe(frad, [0 0 rect_curr]));
	
	jj						= 0;

	% loop along segment until reach end
	while true
		pause(0.1)
    	dist_lim_curr       = ax_radar.XLim + param_anim(2);
		if (dist_lim_curr(end) >= data_cat.dist_lin(end))
        	dist_lim_curr   = [(data_cat.dist_lin(end) - param_anim(3)) data_cat.dist_lin(end)];
			if (dist_lim_curr(1) < 0)
				dist_lim_curr(1) ...
							= 0;
			end
		end
		ind_dist_curr		= interp1(data_cat.dist_lin, 1:data_cat.num_trace, dist_lim_curr, 'nearest', 'extrap');
    	ax_radar.XLim		= dist_lim_curr;
		if isnan(min(elev_bed_curr(ind_dist_curr(1):ind_dist_curr(2)), [], 'omitnan'))
			ax_radar.YLim	= [elev_vec(1) max(elev_surf_grd_curr(ind_dist_curr(1):ind_dist_curr(2)))];
		else
			ax_radar.YLim	= [min(elev_bed_curr(ind_dist_curr(1):ind_dist_curr(2)), [], 'omitnan') max(elev_surf_grd_curr(ind_dist_curr(1):ind_dist_curr(2)))];
		end
		[plc.XData, plc.YData] ...
							= deal(pk_cat.x{ii}(ind_dist_curr), pk_cat.y{ii}(ind_dist_curr));
		ax_radar.CLim		= mean(amp_curr(1:10:end, ind_dist_curr(1):10:ind_dist_curr(2)), 'all', 'omitnan') + (std(amp_curr(1:10:end, ind_dist_curr(1):10:ind_dist_curr(2)), 0, 'all', 'omitnan') .* [-2 2]);
		disp(jj)
		if (jj == 1)
			if strcmp(p_pk(1).Visible, 'on')
				[p_pk(1:pk_cat.num_layer(ii)).Visible] = deal('off');
			else
				[p_pk(1:pk_cat.num_layer(ii)).Visible] = deal('on');
			end
			jj				= 0;
		end
    	writeVideo(mp4, getframe(frad, [0 0 rect_curr]))
		if (dist_lim_curr(end) == data_cat.dist_lin(end)) % reached end of segment
			writeVideo(mp4, getframe(frad, [0 0 rect_curr]))
			break
		end
		jj					= jj + 1;
	end
	
	uistack(p_core_map, 'top')
end
close(mp4)
toc
%%	
end