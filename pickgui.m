 function pickgui(varargin)
% PICKGUI Interactive radar-layer picker.
%   
%   PICKGUI loads a GUI for tracing layers in single or concatenated CReSIS
%   radargrams. This second major version of the GUI includes
%   semi-automatic tracing of layers and layer prediction using independent
%   CWT/Hough transform layer prediction (ARESELP). Layers from
%   intersecting radargrams can later be matched between segments using
%   FENCEGUI.
%   
%   Refer to pickgui_hotkeys.pdf for shortcuts.
%	
%   If loading individual CReSIS data frames, PICKGUI requires the Mapping
%   Toolbox be licensed and available.
%	
%	If the Parallel Computing Toolbox is licensed and available, then
%	several calculations related to data flattening will be parallelized.
%   
%   Including any argument in the PICKGUI call (e.g., "pickgui(1)") speeds
%   up intialization by not starting the parallel pool, even if the
%   Parallel Computing Toolbox is available. However, this turns off
%   parallelization for those loops that can use it.
%   
% Joe MacGregor (NASA)
% Last updated: 12 August 2025

%% Intialize variables

pk							= struct;
pk.layer                    = struct;
[pk.num_layer, pk.ind_trim_start] ...
                            = deal(0, 1);

% two-way traveltime defaults
[twtt_min_ref, twtt_max_ref]= deal(0, 1);
[twtt_min, twtt_max]        = deal(twtt_min_ref, twtt_max_ref);

% distance defaults
[dist_min_ref, dist_max_ref]= deal(0, 1e3);
[dist_min, dist_max]        = deal(dist_min_ref, dist_max_ref);

% dB default
[db_min_ref, db_max_ref]    = deal(-130, 0);
[db_min, db_max]            = deal(-100, -20);

% default values for several parameters
speed_vacuum                = 299792458; % speed of light in the vacuum, m/s
permitt_ice                 = 3.15; % real part of relative permittivity of ice at radio frequencies, dimensionless
speed_ice                   = speed_vacuum / sqrt(permitt_ice); % speed of light in ice, m/s
length_chunk                = 100; % display chunk length, km
length_pk_max				= 25; % maximum length of semi-automatic tracing in either direction, km
pk.num_win                  = 1; % +/- number of vertical indices in window within which to search for peak in adjacent traces
ord_poly                    = 2; % order of polynomial fit for flattening
length_smooth				= 3; % smoothing length for flattening polynomials, km

num_ind_layer_min			= 5; % minimum number of non-NaN values for any given layer
num_win_ref					= pk.num_win;
disp_type                   = 'twtt';
cmaps                       = {bone(2 ^ 8); jet(2 ^ 8)};

pk_color_def				= [0    0       0.75;
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

if (license('checkout', 'distrib_computing_toolbox') && ~nargin)
    pool_check				= gcp('nocreate');
    if isempty(pool_check)
        try
            parpool(4);
        catch
            parpool;
        end
    end
    parallel_check          = true;
else
    parallel_check          = false;
end

% pre-allocate a bunch of variables
[bed_avail, cbfix_check1, clutter_avail, cross_pass, depth_avail, distfix_check, flat_done, grid_check, load_done, load_flat, norm_done, pk_check, pk_done, surf_avail, surfbed_check, trim_done, twttfix_check] ...
                            = deal(false);
[cbfix_check2, cross_check, int_core_check]	...
							= deal(true);
[amp, amp_depth, amp_flat, amp_norm, data_cat, button, curr_chunk, curr_layer, dist_chunk, ii, ind_bed, ind_bed_flat, ind_bed_norm, ind_bed_smooth, ind_num_trace, ind_surf, ind_surf_flat, ind_surf_norm, ind_surf_smooth, ...
 ind_thick_smooth, ind_x_pk, ind_z_curr, ind_z_flat, ind_z_mat, ind_z_pk, jj, kk, num_chunk, num_sample_trim, pk_color, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8] ...
							= deal(NaN);
[p_bed, p_data, p_surf, pk_fig] ...
                            = deal(gobjects(1));
[p_core, p_int, p_pk]		= deal(gobjects(0));
[pk_ind_z, pk_ind_z_norm, pk_ind_z_flat] ...
							= deal([]);
[file_data, file_pk, path_data, path_pk, path_int] ...
                            = deal('');
% path_data					= '/Users/jamacgre/OneDrive - NASA/data/AntArchitecture/ant_cat/';
% path_pk						= '/Users/jamacgre/OneDrive - NASA/data/AntArchitecture/ant_cat_pk/';
% path_int					= '/Users/jamacgre/OneDrive - NASA/research/matlab/pickgui_v2/mat/';
cmap_curr					= 1;
narrow_ax					= [dist_min dist_max twtt_min twtt_max];

if ~isempty(path_int)
	load([path_int 'int_cat_core_ant.mat'], 'file_cat', 'int_all', 'int_core')
end

if ispc % windows switch GUI element size
    size_font               = 14;
    width_slide_z           = 0.01;
else
    size_font               = 18;
    width_slide_z           = 0.02;
end

%% Draw the GUI

set(0, 'DefaultFigureWindowStyle', 'docked')
pk_gui                      = figure('Toolbar', 'figure', 'Name', 'PICKGUI_v2', 'MenuBar', 'none', 'KeyPressFcn', @keypress, 'WindowButtonDownFcn', @mouse_click);
ax_radar                    = subplot('Position', [0.065 0.06 0.86 0.81]);
disableDefaultInteractivity(ax_radar)
hold on
colormap(bone)
clim([db_min db_max])
axis ij tight
[ax_radar.FontSize, ax_radar.Layer, ax_radar.Toolbar.Visible, ax_radar.LineWidth, ax_radar.GridColor, ax_radar.GridAlpha, ax_radar.GridLineStyle] = deal(size_font, 'top', 'off', 1, 'w', 0.75, '--');
xlabel('Distance (km)')
ylabel('Traveltime ({\mu}s), 1 {\mu}s ~ 85 m in ice')
colorbar('FontSize', size_font)
box on
% pan/zoom callbacks
h_pan                       = pan;
h_zoom                      = zoom;
[h_pan.ActionPostCallback, h_zoom.ActionPostCallback] = deal(@pan_zoom);

% sliders
twtt_min_slide              = uicontrol(pk_gui, 'Style', 'slider', 'Units', 'normalized', 'Position', [0.005 0.50 width_slide_z 0.32], 'CallBack', @slide_twtt_min, 'Min', 0, 'Max', 1, 'Value', twtt_max_ref, 'sliderstep', [0.01 0.1]);
twtt_max_slide              = uicontrol(pk_gui, 'Style', 'slider', 'Units', 'normalized', 'Position', [0.005 0.07 width_slide_z 0.32], 'CallBack', @slide_twtt_max, 'Min', 0, 'Max', 1, 'Value', twtt_min_ref, 'sliderstep', [0.01 0.1]);
cb_min_slide                = uicontrol(pk_gui, 'Style', 'slider', 'Units', 'normalized', 'Position', [0.96 0.07 width_slide_z 0.32], 'CallBack', @slide_db_min, 'Min', -150, 'Max', 0, 'Value', db_min_ref, 'sliderstep', [0.01 0.1]);
cb_max_slide                = uicontrol(pk_gui, 'Style', 'slider', 'Units', 'normalized', 'Position', [0.96 0.50 width_slide_z 0.32], 'CallBack', @slide_db_max, 'Min', -150, 'Max', 0, 'Value', db_max_ref, 'sliderstep', [0.01 0.1]);
dist_min_slide              = uicontrol(pk_gui, 'Style', 'slider', 'Units', 'normalized', 'Position', [0.12 0.005 0.27 0.02], 'CallBack', @slide_dist_min, 'Min', 0, 'Max', 1, 'Value', (1e-3 * dist_min_ref), 'sliderstep', [0.01 0.1]);
dist_max_slide              = uicontrol(pk_gui, 'Style', 'slider', 'Units', 'normalized', 'Position', [0.64 0.005 0.27 0.02], 'CallBack', @slide_dist_max, 'Min', 0, 'Max', 1, 'Value', (1e-3 * dist_max_ref), 'sliderstep', [0.01 0.1]);

% slider values
twtt_min_edit               = annotation('textbox', [0.005 0.82 0.04 0.03], 'String', num2str(twtt_min_ref), 'FontSize', size_font, 'Color', 'k', 'EdgeColor', 'none');
twtt_max_edit               = annotation('textbox', [0.005 0.39 0.04 0.03], 'String', num2str(twtt_max_ref), 'FontSize', size_font, 'Color', 'k', 'EdgeColor', 'none');
cb_min_edit                 = annotation('textbox', [0.965 0.39 0.04 0.03], 'String', num2str(db_min_ref), 'FontSize', size_font, 'Color', 'k', 'EdgeColor', 'none');
cb_max_edit                 = annotation('textbox', [0.9665 0.82 0.04 0.03], 'String', num2str(db_max_ref), 'FontSize', size_font, 'Color', 'k', 'EdgeColor', 'none');
dist_min_edit               = annotation('textbox', [0.08 0.005 0.04 0.03], 'String', num2str(1e-3 * dist_min_ref), 'FontSize', size_font, 'Color', 'k', 'EdgeColor', 'none');
dist_max_edit               = annotation('textbox', [0.59 0.005 0.04 0.03], 'String', num2str(1e-3 * dist_max_ref), 'FontSize', size_font, 'Color', 'k', 'EdgeColor', 'none');

% push buttons, comment out unneeded GUI elements for faster performance (quicker refresh)
push_param					= {	'Load data (l)'				[0.005 0.965 0.06 0.03]		@load_data		'b';
								'Trim z'					[0.07 0.965 0.04 0.03]		@trim_z			'b';
								'Load picks (i)'			[0.30 0.965 0.06 0.03]		@pk_load		'b';
								'Pop fig. (q)'				[0.92 0.925 0.04 0.03]		@pop_fig		'g';
								'Flatten (f)'				[0.42 0.925 0.04 0.03]		@flatten		'm';
								'Semi-automatically (p/3)'	[0.46 0.925 0.09 0.03]		@pk_auto		'm';
								'Delete (d)'				[0.80 0.885 0.05 0.03]		@pk_del			'r';
								'Adjust (a)'				[0.7625 0.925 0.0425 0.03]	@pk_adj			'm';
								'Merge (m)'					[0.805 0.925 0.04 0.03]		@pk_merge		'm';
								'Focus (c)'					[0.715 0.925 0.045 0.03]	@pk_focus		'm';
								'Surf./bed (b/4)'			[0.72 0.885 0.05 0.03]		@pk_surfbed		'm';
								'Next (n)'					[0.6775 0.925 0.035 0.03]	@pk_next		'm';
								'Last (b)'					[0.64 0.925 0.035 0.03]		@pk_last		'm';
								'Split (x)'					[0.8475 0.925 0.035 0.03]	@pk_split		'm';
								'Shift z (h)'				[0.885 0.925 0.035 0.03]	@pk_shift		'm';
								'Test (t)'					[0.855 0.885 0.035 0.03]	@misctest		'r';
								'Cross (v)'					[0.96 0.965 0.035 0.03]		@pk_cross		'm';
								'Save (s)'					[0.965 0.925 0.03 0.03]		@pk_save		'g';
								'Reset x/y (e)'				[0.945 0.885 0.05 0.03]		@reset_xz		'r';
								'Reset dB'					[0.955 0.03 0.04 0.03]		@reset_db		'r';
								'Zoom on (z)'				[0.21 0.885 0.045 0.03]		@zoom_on		'b';
								'Zoom off'					[0.26 0.885 0.04 0.03]		@zoom_off		'b';};

for ii = 1:size(push_param, 1)
	uicontrol(pk_gui, 'Style', 'pushbutton', 'String', push_param{ii, 1}, 'Units', 'normalized', 'Position', push_param{ii, 2}, 'CallBack', push_param{ii, 3}, 'FontSize', size_font, 'ForegroundColor', push_param{ii, 4})
end
areselp_push                = uicontrol(pk_gui, 'Style', 'pushbutton', 'String', 'ARESELP', 'Units', 'normalized', 'Position', [0.30 0.925 0.05 0.03], 'CallBack', @use_areselp, 'FontSize', size_font, 'ForegroundColor', 'm', 'Visible', 'off');

% fixed text annotations
text_param					= { [0.21 0.965 0.03 0.03]		'Chunk'			'b';
								[0.21 0.925 0.03 0.03]		'L_{chunk}'		'b';
								[0.965 0.42 0.03 0.03]		'dB_{min}'		'k';
								[0.965 0.85 0.03 0.03]		'dB_{max}'		'k';
								[0.56 0.925 0.03 0.03]		'N_{win}'		'm';
								[0.30 0.925 0.03 0.03]		'Grid (g)'		'k';
								[0.005 0.85 0.03 0.03]		't_{min}'		'k';
								[0.005 0.42 0.03 0.03]		't_{max}'		'k';
								[0.04 0.005 0.03 0.03]		'dist_{min}'	'k';
								[0.55 0.005 0.03 0.03]		'dist_{max}'	'k';
								[0.45 0.965 0.03 0.03]		'Pick'			'm';
								[0.64 0.88 0.04 0.03]		'Layer'			'm';
								[0.56 0.885 0.04 0.03]		'L_{pk,max}'	'm'};
for ii = 1:size(text_param, 1)
	if ~ispc
		annotation('textbox', text_param{ii, 1}, 'String', text_param{ii, 2}, 'FontSize', size_font, 'Color', text_param{ii, 3}, 'EdgeColor', 'none', 'FontWeight', 'bold')
	else
		annotation('textbox', text_param{ii, 1}, 'String', text_param{ii, 2}, 'FontSize', size_font, 'Color', text_param{ii, 3}, 'EdgeColor', 'none')
	end
end

% variable text annotations
file_box                    = annotation('textbox', [0.005 0.925 0.20 0.03], 'String', '', 'Color', 'k', 'FontSize', size_font, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'interpreter', 'none', 'LineWidth', 1);
status_box                  = annotation('textbox', [0.64 0.965 0.315 0.03], 'String', '', 'Color', 'k', 'FontSize', size_font, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'interpreter', 'none', 'LineWidth', 1);
cbl                         = annotation('textbox', [0.93 0.03 0.03 0.03], 'String', '', 'FontSize', size_font, 'Color', 'k', 'EdgeColor', 'none', 'FontWeight', 'bold');

% value boxes
length_chunk_edit           = uicontrol(pk_gui, 'Style', 'edit', 'String', num2str(length_chunk), 'Units', 'normalized', 'Position', [0.25 0.925 0.03 0.03], 'FontSize', size_font, 'ForegroundColor', 'k', 'BackgroundColor', 'w', ...
												'CallBack', @adj_length_chunk);
length_pk_max_edit          = uicontrol(pk_gui, 'Style', 'edit', 'String', num2str(length_pk_max), 'Units', 'normalized', 'Position', [0.595 0.885 0.03 0.03], 'FontSize', size_font, 'ForegroundColor', 'k', 'BackgroundColor', 'w', ...
												'CallBack', @adj_length_pk_max);
num_win_edit                = uicontrol(pk_gui, 'Style', 'edit', 'String', num2str(pk.num_win), 'Units', 'normalized', 'Position', [0.595 0.925 0.03 0.03], 'FontSize', size_font, 'ForegroundColor', 'k', 'BackgroundColor', 'w', ...
												'CallBack', @adj_num_win);

% menus
chunk_list                  = uicontrol(pk_gui, 'Style', 'popupmenu', 'String', 'N/A', 'Value', 1, 'Units', 'normalized', 'Position', [0.25 0.955 0.045 0.04], 'FontSize', size_font, 'ForegroundColor', 'k', 'CallBack', @plot_chunk);
layer_list                  = uicontrol(pk_gui, 'Style', 'popupmenu', 'String', 'N/A', 'Value', 1, 'Units', 'normalized', 'Position', [0.665 0.875 0.05 0.04], 'FontSize', size_font, 'ForegroundColor', 'k', 'CallBack', @pk_select);

% check boxes
length_pk_max_check         = uicontrol(pk_gui, 'Style', 'checkbox', 'Units', 'normalized', 'Position', [0.625 0.885 0.01 0.02], 'FontSize', size_font, 'Value', 1, 'BackgroundColor', pk_gui.Color);

% display buttons
disp_group                  = uibuttongroup('Position', [0.005 0.885 0.20 0.03], 'SelectionChangeFcn', @disp_radio);
uicontrol(pk_gui, 'Style', 'text', 'Parent', disp_group, 'Units', 'normalized', 'Position', [0 0.6 0.9 0.3], 'FontSize', size_font)
disp_check_param			= {'twtt'		0.01;
							   'depth'		0.15;
							   'norm'		0.33;
							   'clutter'	0.55;
							   'flat'		0.77};
disp_check					= gobjects(1, size(disp_check_param, 1));
for ii = 1:size(disp_check_param, 1) %#ok<*FXUP>
	disp_check(ii)          = uicontrol(pk_gui, 'Style', 'radio', 'String', disp_check_param{ii, 1}, 'Units', 'normalized', 'Position', [disp_check_param{ii, 2} 0.1 0.2 0.8], 'Parent', disp_group, 'FontSize', size_font, ...
		'HandleVisibility', 'off');
end
disp_group.SelectedObject = disp_check(1);
[disp_check(2:end).Visible] = deal('off');

%% Clear plots

    function clear_plots(src, event)
        delete(p_pk(isgraphics(p_pk)))
		surfbed_check = true;
        pk_check = false;
        disp_group.SelectedObject = disp_check(1);
        [disp_check(2:end).Visible] = deal('off');
        layer_list.String = 'N/A'; layer_list.Value = 1;
        status_box.EdgeColor = 'k'; status_box.LineWidth = 1;
		pause(0.1)
    end

%% Clear data and picks

    function clear_data(src, event)
        pk					= struct;
        pk.layer            = struct;
        [pk.num_layer, amp, amp_depth, amp_norm, button, curr_chunk, curr_layer, dist_chunk, ii, ind_bed, ind_bed_flat, ind_bed_norm, ind_bed_smooth, ind_num_trace, ind_surf, ind_surf_flat, ind_surf_norm, ...
		 ind_surf_smooth, ind_thick_smooth, ind_x_pk, ind_z_curr, ind_z_flat, ind_z_mat, ind_z_pk, jj, kk, num_chunk, num_sample_trim, pk_color, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8] ...
                            = deal(0);
		pk.ind_trim_start	= 1;
        [bed_avail, clutter_avail, cross_pass, depth_avail, flat_done, load_done, load_flat, norm_done, pk_done, surf_avail, trim_done] ...
                            = deal(false);
		[cross_check, int_core_check] ...
							= deal(true);
        [pk_ind_z, pk_ind_z_norm, pk_ind_z_flat] ...
							= deal([]);
		amp_flat		    = NaN;
		[p_bed, p_data, p_surf, pk_fig] ...
                            = deal(gobjects(1));
		[p_core, p_int, p_pk] ...
							= deal(line([], []));
        file_pk				= '';
        pk.num_win		    = num_win_ref;
        num_win_edit.String = num2str(pk.num_win);
    end

%% Load radar data

    function load_data(src, event) %#ok<*INUSD>
        
        [pk_gui.KeyPressFcn, pk_gui.WindowButtonDownFcn] = deal('');
        
        tmp1                = file_data;
        
        % see if user just wants to move on to next segment
        if ~isempty(file_data)
			tmp2			= dir([path_data '*.mat']);
			tmp2			= {tmp2.name};
			tmp3			= find(contains(tmp2, file_data));
			if (tmp3 < length(tmp2))
				status_box.String = 'Load next segment? Y: yes; otherwise: no...';
				waitforbuttonpress
				if strcmpi(pk_gui.CurrentCharacter, 'Y')
					file_data ...
							= tmp2{tmp3 + 1};
				end
			end
        end
		
        % dialog box to choose radar data file to load
        if (isempty(file_data) || strcmp(tmp1, file_data))
            [tmp1, tmp2]    = deal(file_data, path_data);
            if ~isempty(path_data)
                [file_data, path_data] = uigetfile('*.mat', 'Load radar data:', path_data);
            elseif ~isempty(path_pk)
                [file_data, path_data] = uigetfile('*.mat', 'Load radar data:', path_pk);
            else
                [file_data, path_data] = uigetfile('*.mat', 'Load radar data:');
            end
            if ~ischar(path_data)
                [file_data, path_data] ...
                            = deal('', tmp2);
            end
        end
        
		% loading failed early so bail
        if isempty(file_data)
            file_data       = tmp1;
            status_box.String = 'No data loaded.';
            pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
            return
        end
        
        % attempt to load data
        status_box.String	= 'Loading data...';
        pause(0.1)
		
		tmp1				= load([path_data file_data]);
		
		if isfield(tmp1, 'amp')
			
			data_cat		= tmp1;
			tmp1			= 0;
			file_box.String = file_data(1:(end - 4));
            
		elseif isfield(tmp1, 'Data')
			
			status_box.String = 'Loading possible CReSIS data frame...';
			pause(0.1)
			
            data_cat		= struct;
            
			try
                [data_cat.amp, data_cat.lat, data_cat.lon, data_cat.num_trace, data_cat.twtt, data_cat.elev_air, data_cat.num_sample, data_cat.file_in] ...
							= deal(single(tmp1.Data), tmp1.Latitude, tmp1.Longitude, length(tmp1.Latitude), tmp1.Time, tmp1.Elevation, length(tmp1.Time), file_data(1:(end - 4)));
				file_box.String ...
							= file_data(1:(end - 4));
            catch
                status_box.String = 'Selected file does not contain expected variables.';
                pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
                return
			end
			
			try
                data_cat.twtt_surf ...
							= tmp1.Surface;
            catch
                data_cat.twtt_surf ...
							= NaN(1, data_cat.num_trace);
                status_box.String = 'Selected file does not contain the Surface variable. Setting surface pick to NaN.';
			end
			
			try
                data_cat.twtt_bed ...
							= tmp1.Bottom;
            catch
                data_cat.twtt_bed ...
							= NaN(1, data_cat.num_trace);
                status_box.String = 'Selected file does not contain the Bottom variable. Setting bed pick to NaN.';
			end
			
			% add cluttergram if available			
			if isfield(tmp1, 'Clutter')
				data_cat.clutter ...
							= tmp1.Clutter;
				clutter_avail ...
							= true;
			end
			
        	% convert to dB
        	data_cat.amp((data_cat.amp == 0) | isinf(data_cat.amp)) ...
                            = NaN;
			data_cat.amp	= 10 .* log10(abs(data_cat.amp));
			
			if isrow(data_cat.twtt) % want column-vector fast-time/traveltime
                data_cat.twtt ...
							= data_cat.twtt';
			end
            data_cat.dist	= cumsum([0 distance([data_cat.lat(1:(end - 1))' data_cat.lon(1:(end - 1))'], [data_cat.lat(2:end)' data_cat.lon(2:end)'], wgs84Ellipsoid)']); % great-circle distance, m
            data_cat.dist_lin ...
							= interp1([1 data_cat.num_trace], data_cat.dist([1 end]), 1:data_cat.num_trace);
			[data_cat.ind_trim_surf, data_cat.ind_trim_bed] ...
							= deal(NaN);
            tmp1			= 0;
            
        else
            status_box.String = 'Selected file does not contain expected variables.';
            pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
            return
		end
		
        if load_done
			if isgraphics(p_data)
				delete(p_data)
			end
            if isgraphics(p_bed)
                delete(p_bed)
            end
            if isgraphics(p_surf)
                delete(p_surf)
            end
            [disp_check(2:end).Visible, areselp_push.Visible] = deal('off');
            clear_plots
			if any(isgraphics(p_core))
				delete(p_core(isgraphics(p_core)))
			end
			if any(isgraphics(p_int))
				delete(p_int(isgraphics(p_int)))
			end
            clear_data
            load_done       = false;
			disp_type       = 'twtt';
            twttfix_check	= false;
			disp_group.SelectedObject = disp_check(1);
			if isfield(data_cat, 'Clutter')
				clutter_avail ...
							= true;
			end			
            pause(0.1)
        end
		
        if isfield(data_cat, 'layer_ARESELP') % no ARESELP
            areselp_push.Visible = 'on';
        end
		
		amp					= data_cat.amp;
		data_cat			= rmfield(data_cat, 'amp');
        num_sample_trim     = data_cat.num_sample;
		ind_num_trace		= 1:data_cat.num_trace;
		
        % make chunks
        adj_length_chunk
        
        % assign traveltime and distance reference values/sliders based on data
        [twtt_min_ref, twtt_max_ref, twtt_min, twtt_max, db_min_ref, db_max_ref, db_min, db_max, dist_min_ref, dist_max_ref, dist_min, dist_max] ...
                            = deal(data_cat.twtt(1), data_cat.twtt(end), data_cat.twtt(1), data_cat.twtt(end), min(amp, [], 'all'), max(amp, [], 'all'), min(amp, [], 'all'), max(amp, [], 'all'), data_cat.dist_lin(1), ...
								   data_cat.dist_lin(end), data_cat.dist_lin(1), data_cat.dist_lin(end));
        [twtt_min_slide.Min, twtt_max_slide.Min] = deal(1e6 * twtt_min_ref); [twtt_min_slide.Max, twtt_max_slide.Max] = deal(1e6 * twtt_max_ref); twtt_min_slide.Value = (1e6 * twtt_max_ref); twtt_max_slide.Value = (1e6 * twtt_min_ref);
        twtt_min_edit.String = sprintf('%3.1f', (1e6 * twtt_min_ref)); twtt_max_edit.String = sprintf('%3.1f', (1e6 * twtt_max_ref));		
        [dist_min_slide.Min, dist_max_slide.Min] = deal(1e-3 * dist_min_ref); [dist_min_slide.Max, dist_max_slide.Max] = deal(1e-3 * dist_max_ref); dist_min_slide.Value = 1e-3 * dist_min_ref; dist_max_slide.Value = 1e-3 * dist_max_ref;
        dist_min_edit.String = sprintf('%3.1f', (1e-3 * dist_min_ref)); dist_max_edit.String = sprintf('%3.1f', (1e-3 * dist_max_ref));
        [cb_min_slide.Min, cb_max_slide.Max] = deal(db_min_ref); [cb_min_slide.Max, cb_max_slide.Max] = deal(db_max_ref); cb_min_slide.Value = db_min; cb_max_slide.Value = db_max;
		
		update_dist_range
        update_twtt_range
        
        [ind_surf, ind_bed] = deal(NaN(1, data_cat.num_trace));
        if isfield(data_cat, 'twtt_surf')
            if ~isempty(find((~isnan(data_cat.twtt_surf) & ~isinf(data_cat.twtt_surf)), 1))
                ind_surf(~isnan(data_cat.twtt_surf) & ~isinf(data_cat.twtt_surf)) ...
                            = interp1(data_cat.twtt, 1:num_sample_trim, data_cat.twtt_surf(~isnan(data_cat.twtt_surf) & ~isinf(data_cat.twtt_surf)), 'nearest', 'extrap');
                if ~isempty(find(isnan(ind_surf), 1))
                    ind_surf(isnan(ind_surf)) ...
                            = round(interp1(find(~isnan(ind_surf)), ind_surf(~isnan(ind_surf)), find(isnan(ind_surf)), 'linear', 'extrap'));
                end
                ind_surf(ind_surf <= 0) ...
                            = 1;
                [surf_avail, surfbed_check] = deal(true);
                p_surf      = line((1e-3 .* data_cat.dist_lin), (1e6 .* data_cat.twtt_surf), 'LineStyle', 'none', 'Marker', '.', 'Color', 'm', 'MarkerSize', 8, 'Visible', 'off');
            	amp_depth   = NaN(num_sample_trim, data_cat.num_trace, 'single');
            	for ii = find(~isnan(ind_surf))
                	amp_depth(1:(num_sample_trim - ind_surf(ii) + 1), ii) ...
                        	= amp(ind_surf(ii):num_sample_trim, ii); % shift data up to surface
            	end
            	disp_check(2).Visible = 'on';
            	depth_avail = true;
            else
                p_surf      = line(NaN, NaN, 'LineStyle', 'none', 'Marker', '.', 'Color', 'm', 'MarkerSize', 8, 'Visible', 'off');
            end
        else
            p_surf			= line(NaN, NaN, 'LineStyle', 'none', 'Marker', '.', 'Color', 'm', 'MarkerSize', 8, 'Visible', 'off');
        end
        
        if isfield(data_cat, 'twtt_bed')
            if ~isempty(find((~isnan(data_cat.twtt_bed) & ~isinf(data_cat.twtt_bed)), 1))
                ind_bed(~isnan(data_cat.twtt_bed) & ~isinf(data_cat.twtt_bed)) ...
                            = interp1(data_cat.twtt, 1:num_sample_trim, data_cat.twtt_bed(~isnan(data_cat.twtt_bed) & ~isinf(data_cat.twtt_bed)), 'nearest', 'extrap');
                ind_bed(ind_bed > num_sample_trim) ...
                            = num_sample_trim;
                bed_avail   = true;
                surfbed_check = true;
                p_bed       = line((1e-3 .* data_cat.dist_lin), (1e6 .* data_cat.twtt_bed), 'LineStyle', 'none', 'Marker', '.', 'Color', 'm', 'MarkerSize', 8, 'Visible', 'off');
            else
                p_bed       = line(NaN, NaN, 'LineStyle', 'none', 'Marker', '.', 'Color', 'm', 'MarkerSize', 8, 'Visible', 'off');
            end
        else
            p_bed           = line(NaN, NaN, 'LineStyle', 'none', 'Marker', '.', 'Color', 'm', 'MarkerSize', 8, 'Visible', 'off');
        end
		
		if (surf_avail && bed_avail) % thickness-normalization can be done
			disp_check(3).Visible = 'on';
			do_norm
		end
		
		if (~isnan(data_cat.ind_trim_surf) && ~isnan(data_cat.ind_trim_bed))
			trim_done		= true;
		end
		
		if ~isempty(path_int)
			tmp1			= find(contains(file_cat, file_data));
			if ~isempty(tmp1)
				tmp2		= [int_all((int_all(:, 1) == tmp1), 2); int_all((int_all(:, 3) == tmp1), 4)];
				if ~isempty(tmp2)
					p_int	= gobjects(1, length(tmp2));
					for ii = 1:length(tmp2)
						p_int(ii) ...
							= line((1e-3 .* data_cat.dist_lin(tmp2([ii ii]))), (1e6 .* [twtt_min_ref twtt_max_ref]), 'Color', 'w', 'LineWidth', 2, 'LineStyle', '--', 'Visible', 'on');
					end
				end
				tmp2		= int_core((int_core(:, 1) == tmp1), 2);
				if ~isempty(tmp2)
					p_core	= gobjects(1, length(tmp2));
					for ii = 1:length(tmp2)
						p_core(ii) ...
							= line((1e-3 .* data_cat.dist_lin(tmp2([ii ii]))), (1e6 .* [twtt_min_ref twtt_max_ref]), 'Color', 'm', 'LineWidth', 2, 'LineStyle', '--', 'Visible', 'on');
					end
				end
			end
		end
		
		if clutter_avail
            disp_check(4).Visible = 'on';
            data_cat.clutter(isinf(data_cat.clutter)) ...
                            = NaN;
			data_cat.clutter= interp1([min(data_cat.clutter, [], 'all') max(data_cat.clutter, [], 'all')], [db_min db_max], data_cat.clutter);
			if (num_sample_trim > size(data_cat.clutter, 1))
                data_cat.clutter ...
							= [data_cat.clutter; NaN((num_sample_trim - size(data_cat.clutter, 1)), size(data_cat.clutter, 2))];
			end
		end
		
        layer_list.String = {'surface' 'bed'}; layer_list.Value = 1;
        pk_select
        
        % plot data
        disp_group.SelectedObject = disp_check(1);
        load_done           = true;
        plot_twtt
        status_box.String = 'Data loaded successfully.';
		pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
    end

%% Trim excess data within GUI (e.g., before surface reflection and after deepest bed reflection)

    function trim_z(src, event)
        
        if ~load_done
            status_box.String = 'No data to trim.';
            return
        end
        if (trim_done && flat_done)
            status_box.String = 'Data that has already been flattened should not be trimmed twice.';
            return
        end
        if ~strcmp(disp_type, 'twtt')
            disp_group.SelectedObject = disp_check(1);
            disp_type       = 'twtt';
            plot_twtt
        end
        
        status_box.String = 'Trimming data...';
        pause(0.1)
        
        % trim column vector
        tmp1                = (interp1(data_cat.twtt, 1:num_sample_trim, twtt_min, 'nearest', 'extrap'):interp1(data_cat.twtt, 1:num_sample_trim, twtt_max, 'nearest', 'extrap'))';
        if trim_done
            tmp2            = pk.ind_trim_start + tmp1(1) - 1;
        end
        pk.ind_trim_start   = tmp1(1);
        
        % adjust and trim surface indices
        if surf_avail
            ind_surf        = ind_surf - pk.ind_trim_start + 1;
            ind_surf(ind_surf < 1) ...
                            = NaN;
        end
        
        % trim data
        amp					= amp(tmp1, :);
        if depth_avail
            tmp3            = amp_depth;
            amp_depth       = NaN(size(amp), 'single');
            for ii = find(~isnan(ind_surf))
                amp_depth(1:(length(tmp1) - ind_surf(ii) + 1), ii) ...
                            = tmp3(1:(length(tmp1) - ind_surf(ii) + 1), ii);
            end
        end
        
        num_sample_trim     = length(tmp1);
        
        % adjust and trim surface and bed indices
        if surf_avail
            ind_surf(ind_surf > num_sample_trim) ...
                            = NaN;
        end
        if bed_avail
            ind_bed         = ind_bed - pk.ind_trim_start + 1;
            ind_bed((ind_bed < 1) | (ind_bed > num_sample_trim)) ...
                            = NaN;
        end
        
        % if layers already exist, then trim them too
		if pk_done
			pk_ind_z		= pk_ind_z - pk.ind_trim_start + 1;
		end
		
		if trim_done
            pk.ind_trim_start ...
                            = tmp2;
        else
            trim_done       = true;
		end
		
        % save traveltime for the end
        data_cat.twtt       = data_cat.twtt(tmp1);
		
		if norm_done
			do_norm
			if pk_done
				tmp1		= (1:num_sample_trim)'; % normalization vector
				tmp2		= (tmp1(:, ones(1, data_cat.num_trace)) - ind_surf_smooth(ones(num_sample_trim, 1), :)) ./ ind_thick_smooth(ones(num_sample_trim, 1), :); % normalization matrix
				for ii = 1:pk.num_layer
					for jj = find(~isnan(pk_ind_z(ii, :)) & ~isnan(tmp2(1, :)))
						pk_ind_z_norm(ii, jj) ...
							= interp1(tmp1, tmp2(:, jj), pk_ind_z(ii, jj));
					end
				end
				pk_ind_z_norm ...
							= round(pk_ind_z_norm .* num_sample_trim);
				pk_ind_z_norm((pk_ind_z_norm < 1) | (pk_ind_z_norm > num_sample_trim)) ...
							= NaN;
			end
		end
		
		[p_data.YData, p_data.CData] = deal((1e6 .* data_cat.twtt), amp);
        [twtt_min_ref, twtt_max_ref] ...
                            = deal(data_cat.twtt(1), data_cat.twtt(end));
        [twtt_min_slide.Min, twtt_max_slide.Min] = deal(1e6 * twtt_min_ref);  [twtt_min_slide.Max, twtt_max_slide.Max] ...
							= deal(1e6 * twtt_max_ref); twtt_min_slide.Value = (1e6 * twtt_max_ref); twtt_max_slide.Value = (1e6 * twtt_min_ref);
		
		if any(isgraphics(p_int))
			for ii = find(isgraphics(p_int))
				p_int(ii).YData = 1e6 .* [twtt_min_ref twtt_max_ref];
			end
		end
		if any(isgraphics(p_core))
			for ii = find(isgraphics(p_core))
				p_core(ii).YData = 1e6 .* [twtt_min_ref twtt_max_ref];
			end
		end

        update_twtt_range
		update_pk_plot
        status_box.String = ['Data trimmed off before ' num2str((1e6 * data_cat.twtt(1)), '%2.1f') ' us and after ' num2str((1e6 * data_cat.twtt(end)), '%2.1f') ' us.'];
    end

%% Use ARESELP predicted layers to get started

	function use_areselp(src, event)
		
		if ~load_done
            status_box.String = 'Load data.';
            return
		end
		if pk_done
			status_box.String = 'Cannot use ARESELP until all layers are deleted.';
			return
		end
		if ~strcmp(disp_type, 'twtt')
            disp_group.SelectedObject = disp_check(1);
            disp_type       = 'twtt';
            plot_twtt
		end
		
        [pk_gui.KeyPressFcn, pk_gui.WindowButtonDownFcn] = deal('');
		
		pk.num_layer		= data_cat.num_layer_ARESELP;
		pk_ind_z			= reshape([data_cat.layer_ARESELP(:).ind_z], data_cat.num_trace, pk.num_layer)';
		
		update_pk_color
		p_pk				= deal(gobjects(1, pk.num_layer));
        [pk_ind_z_flat, pk_ind_z_norm] ...
							= deal(NaN(pk.num_layer, data_cat.num_trace));
        tmp1                = [];
        for ii = 1:pk.num_layer
			if ~isempty(find(~isnan(pk_ind_z(ii, :)), 1))
                try
					p_pk(ii)= line((1e-3 .* data_cat.dist_lin(~isnan(pk_ind_z(ii, :)))), (1e6 .* data_cat.twtt(pk_ind_z(ii, ~isnan(pk_ind_z(ii, :))))'), ...
								   'LineStyle', 'none', 'Marker', '.', 'Color', pk_color(ii, :), 'MarkerSize', 8, 'Visible', 'off');
                catch
                    tmp1    = [tmp1 ii]; %#ok<AGROW>
                    continue
                end
            else
                tmp1        = [tmp1 ii]; %#ok<AGROW>
                continue
			end
        end
        
        % get rid of empty layers
        if ~isempty(tmp1)
            tmp2            = setdiff(1:pk.num_layer, tmp1);
            for ii = tmp1
                if isgraphics(p_pk(ii))
                    delete(p_pk(ii))
                end
            end
            if (length(tmp2) > 1)
                curr_layer  = interp1(1:length(tmp2), tmp2, curr_layer, 'nearest', 'extrap');
            else
                curr_layer  = 1;
            end
            [pk.num_layer, pk_ind_z, p_pk, pk_color] ...
                            = deal(length(tmp2), pk_ind_z(tmp2, :), p_pk(tmp2), pk_color(tmp2, :));
            layer_list.String = [num2cell(1:pk.num_layer) 'surface' 'bed']; layer_list.Value = curr_layer;
        end
		
		if flat_done
			flat_done		= false;
		end
		[pk_check, pk_done] = deal(true);
        show_pk
		pk_sort
		pk_select
        pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
		status_box.String = 'ARESELP layers loaded as new layers.';
	 end

%% Load existing layer picks for this dataset

    function pk_load(src, event)
        
        if ~load_done
            status_box.String = 'Load data before picks.';
            return
        end
        if ~strcmp(disp_type, 'twtt')
            disp_group.SelectedObject = disp_check(1);
            disp_type       = 'twtt';
            plot_twtt
        end
        
        [pk_gui.KeyPressFcn, pk_gui.WindowButtonDownFcn] = deal('');
        
        if exist([path_data file_data(1:(end - 4)) '_pk.mat'], 'file')
			path_pk			= path_data;
            file_pk         = [file_data(1:(end - 4)) '_pk.mat'];
		elseif (~isempty(path_pk) && exist([path_pk file_data(1:(end - 4)) '_pk.mat'], 'file'))
			file_pk         = [file_data(1:(end - 4)) '_pk.mat'];
        else % Dialog box to choose picks file to load
            if ~isempty(path_pk)
                [file_pk, path_pk] = uigetfile('*.mat', 'Load picks:', path_pk);
            elseif ~isempty(path_data)
                [file_pk, path_pk] = uigetfile('*.mat', 'Load picks:', path_data);
            else
                [file_pk, path_pk] = uigetfile('*.mat', 'Load picks:');
            end
            if ~ischar(file_pk)
                [file_pk, path_pk] = deal('', path_data);
            end
        end
        
        if isempty(file_pk)
            status_box.String = 'No picks loaded.';
            pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
            return
        end
        
        % check against data file
        try %#ok<TRYNC>
            if ~strcmp(file_data(1:(end - 4)), file_pk(1:(end - 7)))
                status_box.String = ['Selected picks file (' file_pk(1:(end - 4)) ') might not match data file. Continue loading? Y: yes; otherwise: no...'];
                waitforbuttonpress
                if ~strcmpi(pk_gui.CurrentCharacter, 'Y')
                    pk_gui.KeyPressFcn =@keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
                    status_box.String = 'Loading of picks file canceled.';
                    return
                end
            end
        end
        
        % load picks file
        tmp1                = load([path_pk file_pk]);
        try
            pk              = tmp1.pk;
            tmp1            = 0;
        catch
            status_box.String = 'Selected file does not contained a pk structure.';
            pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
            return
        end
        
        if ~pk.num_layer
            status_box.String = 'Picks file has no picks.';
            pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
            return
        end
        
        clear_plots
        if (surf_avail || bed_avail)
            surfbed_check = true;
        end
        
		% replace surface and bed with pk values if different from data
        if (isfield(data_cat, 'twtt_surf') && isfield(pk, 'twtt_surf'))
            if (~isequal(data_cat.twtt_surf, pk.twtt_surf) && ~isempty(find(~isnan(pk.twtt_surf), 1)))
                data_cat.twtt_surf ...
                            = pk.twtt_surf;
				if isgraphics(p_surf)
					[p_surf.XData, p_surf.YData] = deal((1e-3 .* data_cat.dist_lin), (1e6 .* data_cat.twtt_surf));
				end
            end
        end
        if (isfield(data_cat, 'twtt_bed') && isfield(pk, 'twtt_bed'))
            if (~isequal(data_cat.twtt_bed, pk.twtt_bed) && ~isempty(find(~isnan(pk.twtt_bed), 1)))
                data_cat.twtt_bed...
                            = pk.twtt_bed;
				if isgraphics(p_bed)
					[p_bed.XData, p_bed.YData] = deal((1e-3 .* data_cat.dist_lin), (1e6 .* data_cat.twtt_bed));
				end
            end
        end
        
        [flat_done, load_flat, pk_done] ...
                            = deal(false);
        
		if isfield(pk, 'poly_flat')
            load_flat       = true;
            ord_poly        = size(pk.poly_flat, 1) - 1;
		end
		
		% check for empty layers and if so remove them
		tmp1				= [];
		for ii = 1:length(pk.layer)
			if isempty(pk.layer(ii).ind_z)
				tmp1		= [tmp1 ii]; %#ok<AGROW> 
			end
		end
		if ~isempty(tmp1)
			pk.layer		= pk.layer(setdiff(1:length(pk.layer), tmp1));
			pk.num_layer	= length(pk.layer);
		end
		
		pk_ind_z			= reshape([pk.layer(:).ind_z], data_cat.num_trace, pk.num_layer)';
		pk.layer			= struct;
		
        num_win_edit.String = num2str(pk.num_win);
        layer_list.String = [num2cell(1:pk.num_layer) 'surface' 'bed']; layer_list.Value = 1;
        
        if ~isfield(pk, 'twtt_min_ref')
            pk.twtt_min_ref = data_cat.twtt(1);
        end
        if ~isfield(pk, 'twtt_max_ref')
            pk.twtt_max_ref = data_cat.twtt(end);
        end
        
        % trim data as it was before
        [twtt_min, twtt_min_ref] ...
                            = deal(pk.twtt_min_ref);
        if surf_avail
            if (twtt_min_ref > min(data_cat.twtt_surf(~isinf(data_cat.twtt_surf))))
				[twtt_min, twtt_min_ref] ...
                            = deal(min(data_cat.twtt_surf(~isinf(data_cat.twtt_surf))));
            end
        end
        [twtt_max, twtt_max_ref] ...
                            = deal(pk.twtt_max_ref);
        if bed_avail
            if (twtt_max_ref < max(data_cat.twtt_bed(~isinf(data_cat.twtt_bed))))
				[twtt_max, twtt_max_ref] ...
                            = deal(max(data_cat.twtt_bed(~isinf(data_cat.twtt_bed))));
            end
        end
        twtt_min_edit.String = sprintf('%3.1f', (1e6 * twtt_min_ref));
        twtt_max_edit.String = sprintf('%3.1f', (1e6 * twtt_max_ref));
		ax_radar.YLim		= 1e6 .* [twtt_min_ref twtt_max_ref];
		
		if any(isgraphics(p_int))
			for ii = find(isgraphics(p_int))
				p_int(ii).YData = 1e6 .* [twtt_min_ref twtt_max_ref];
			end
		end
		if any(isgraphics(p_core))
			for ii = find(isgraphics(p_core))
				p_core(ii).YData = 1e6 .* [twtt_min_ref twtt_max_ref];
			end
		end
		
        pk_done             = true;
        status_box.String = 'Continuing pick loading...';
		pause(0.1)
		
        % plot layers in twtt
        update_pk_color
		p_pk				= deal(gobjects(1, pk.num_layer));
        tmp1                = [];
		for ii = 1:pk.num_layer
			if (length(find(~isnan(pk_ind_z(ii, :)))) >= num_ind_layer_min) 
                try
					p_pk(ii)= line((1e-3 .* data_cat.dist_lin(~isnan(pk_ind_z(ii, :)))), (1e6 .* data_cat.twtt(pk_ind_z(ii, ~isnan(pk_ind_z(ii, :))))'), ...
								   'LineStyle', 'none', 'Marker', '.', 'Color', pk_color(ii, :), 'MarkerSize', 8, 'Visible', 'off');
                catch
                    tmp1    = [tmp1 ii]; %#ok<AGROW>
                    continue
                end
            else
                tmp1        = [tmp1 ii]; %#ok<AGROW>
                continue
			end
		end
		
        % get rid of nearly empty layers
        if ~isempty(tmp1)
            tmp2            = setdiff(1:pk.num_layer, tmp1);
            for ii = tmp1
                if isgraphics(p_pk(ii))
                    delete(p_pk(ii))
                end
            end
            if (length(tmp2) > 1)
                curr_layer  = interp1(1:length(tmp2), tmp2, curr_layer, 'nearest', 'extrap');
            else
                curr_layer  = 1;
            end
            [pk.num_layer, pk_ind_z, p_pk, pk_color] ...
                            = deal(length(tmp2), pk_ind_z(tmp2, :), p_pk(tmp2), pk_color(tmp2, :));
            layer_list.String = [num2cell(1:pk.num_layer) 'surface' 'bed']; layer_list.Value = curr_layer;
        end
        
		if (pk.ind_trim_start > 1)
			trim_done		= true;
			pk.ind_trim_start ...
							= 1;
			trim_z
		elseif ~pk.ind_trim_start
			pk.ind_trim_start ...
							= 1;
		end
		
        pk_check			= true;
		if depth_avail
            disp_check(2).Visible = 'on';
		end
		
		[pk_ind_z_norm, pk_ind_z_flat] ...
							= deal(NaN(pk.num_layer, data_cat.num_trace));
		
        % if load_flat % reflatten if possible
        %     pause(0.1)
        %     status_box.String = 'Re-flattening data...';
        %     flatten
        %     [pk_gui.KeyPressFcn, pk_gui.WindowButtonDownFcn] = deal('');
        % end
		
		if norm_done
			disp_check(3).Visible = 'on';
			tmp1			= (1:num_sample_trim)'; % normalization vector
			tmp2			= (tmp1(:, ones(1, data_cat.num_trace)) - ind_surf_smooth(ones(num_sample_trim, 1), :)) ./ ind_thick_smooth(ones(num_sample_trim, 1), :); % normalization matrix
			for ii = 1:pk.num_layer
				for jj = find(~isnan(pk_ind_z(ii, :)) & ~isnan(tmp2(1, :)))
					pk_ind_z_norm(ii, jj) ...
							= interp1(tmp1, tmp2(:, jj), pk_ind_z(ii, jj));
				end
			end
			pk_ind_z_norm = round(pk_ind_z_norm .* num_sample_trim);
			pk_ind_z_norm((pk_ind_z_norm < 1) | (pk_ind_z_norm > num_sample_trim)) ...
							= NaN;
		end
        pk_cross
        show_pk
		pk_sort
		pk_select
        pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
        if load_flat
            status_box.String = ['Picks loaded from ' file_pk(1:(end - 4)) '.'];
        else
            status_box.String = ['Picks loaded from ' file_pk(1:(end - 4)) ' (no flattening).'];
        end
    end

%% Flatten data

    function flatten(src, event)
        
    	if (pk.num_layer < ord_poly)
        	status_box.String = 'Not enough layers to attempt flattening.';
        	return
    	end
        
		[pk_gui.KeyPressFcn, pk_gui.WindowButtonDownFcn] = deal('');
		
        flat_done           = false;
        
        status_box.String = 'Starting polynomial fitting...';
        pause(0.1)
        
        pk.poly_flat        = NaN((ord_poly + 1), data_cat.num_trace, 'single');
        % 2nd-order polynomial fits between kept layers at pk.ind_x_start (x) and kept layers at each trace (y)
        
		if parallel_check
            pctRunOnAll warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale')
            pctRunOnAll warning('off', 'MATLAB:polyfit:PolyNotUnique')
        else
            warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale')
            warning('off', 'MATLAB:polyfit:PolyNotUnique')
		end
        
        ind_z_curr			= pk_ind_z;
        if surf_avail
            ind_z_curr		= [smoothdata(ind_surf, 'rlowess', ceil(length_smooth / median(diff(1e-3 .* data_cat.dist_lin), 'omitnan'))); ind_z_curr];
        end
		if bed_avail
			ind_z_curr		= [ind_z_curr; smoothdata(ind_bed, 'rlowess', ceil(length_smooth / median(diff(1e-3 .* data_cat.dist_lin), 'omitnan')))];
		end
        
        ind_x_pk			= sum(~isnan(ind_z_curr), 2); % length of each layer
		ind_x_pk			= sum(~isnan(ind_z_curr) .* ind_x_pk(:, ones(1, data_cat.num_trace))); % layer locations multiplied with length of each layer
        ind_x_pk			= find((ind_x_pk == max(ind_x_pk)), 1); % first trace in the record that has the most+longest layers
        ind_z_pk			= ind_z_curr(:, ind_x_pk); % depth of all layers at reference trace (including NaN)
        tmp1				= ind_z_curr(~isnan(ind_z_pk), :); % depth of ALL layers that are not NaN at the reference trace
        tmp2				= tmp1(:, ind_x_pk); % depths of non-NaN layers at reference trace
		
        % polyfit to layers
        tmp3                = find(sum(~isnan(tmp1)) > ord_poly); % traces where it will be worthwhile to do the polynomial
        if (parallel_check && (length(tmp3) >= 1e3))
            tmp1            = tmp1(:, tmp3);
            tmp4            = pk.poly_flat(:, tmp3);
			[tmp6, tmp7]	= deal(cell(1, length(tmp4)));
			for ii = 1:length(tmp4)
				tmp6{ii}	= tmp2(~isnan(tmp1(:, ii)));
				tmp7{ii}	= tmp1(~isnan(tmp1(:, ii)), ii);
			end
            parfor ii = 1:length(tmp3)
                tmp4(:, ii) = polyfit(tmp6{ii}, tmp7{ii}, ord_poly)';
            end
            pk.poly_flat(:, tmp3) ...
                            = tmp4;
            [tmp4, tmp6, tmp7] ...
							= deal(0);
        else
			for ii = tmp3
                pk.poly_flat(:, ii) ...
                            = polyfit(tmp2(~isnan(tmp1(:, ii))), tmp1(~isnan(tmp1(:, ii)), ii), ord_poly)';
			end
        end
        
		pk.poly_flat		= smoothdata(pk.poly_flat, 2, 'rlowess', ceil(length_smooth / median(diff(1e-3 .* data_cat.dist_lin), 'omitnan')), 'includenan'); % smooth polynomials
		pk.poly_flat		= fillmissing(pk.poly_flat, 'spline', 2, 'EndValues', 'extrap', 'MaxGap', ceil(length_smooth / median(diff(1e-3 .* data_cat.dist_lin), 'omitnan'))); % fill small gaps (if any)
		
        % iterate using other layers if any are available
        if any(isnan(ind_z_pk))
            
            status_box.String = ['Iterating flattening for ' num2str(length(find(isnan(ind_z_pk)))) ' layers not fitted initially...'];
            pause(0.1)
            
            [kk, tmp8]		= deal(0, []);
			
			while true
                
                kk                  = kk + 1;
				
				% determine which layers overlap with original polyfitted layers and order further addition based on the maximum overlap
        		tmp1				= zeros(1, length(ind_z_pk));
				tmp1(isnan(ind_z_pk)) ...
									= sum(~isnan(ind_z_curr(isnan(ind_z_pk), sum(~isnan(ind_z_curr(~isnan(ind_z_pk), :))) > ord_poly)), 2)';
				if isempty(setdiff(find(tmp1 > 0), tmp8))
					break
				end
        		[tmp1, tmp2]		= sort(tmp1, 'descend');
				ii					= tmp2(logical(tmp1));
				if ~isempty(tmp8) % skip tricky layers
					tmp2			= setdiff(tmp2, tmp8, 'stable');
					if isempty(tmp2)
						break
					end
				end
				ii					= tmp2(1); % for each iteration, pick layer with maximum overlap to add
				
                status_box.String = ['Adding layer #' num2str(ii) ' (#' num2str(kk) ' so far) to flattening...'];
                pause(0.1)
                
                % calculate flattening matrix based on current polynomials
                tmp1                = find(~isnan(ind_z_curr(ii, :)));
                ind_z_mat           = single((1:num_sample_trim)');
                ind_z_mat           = ind_z_mat(:, ones(1, length(tmp1))); % matrix of y indices
                switch ord_poly
                    case 2
                        ind_z_flat  = ((ind_z_mat .^ 2) .* pk.poly_flat(ones(num_sample_trim, 1), tmp1)) + (ind_z_mat .* (pk.poly_flat((2 .* ones(num_sample_trim, 1)), tmp1))) + pk.poly_flat((3 .* ones(num_sample_trim, 1)), tmp1);
                    case 3
                        ind_z_flat  = ((ind_z_mat .^ 3) .* pk.poly_flat(ones(num_sample_trim, 1), tmp1)) + ((ind_z_mat .^ 2) .* pk.poly_flat((2 .* ones(num_sample_trim, 1)), tmp1)) + ...
                                      (ind_z_mat .* (pk.poly_flat((3 .* ones(num_sample_trim, 1)), tmp1))) + pk.poly_flat((4 .* ones(num_sample_trim, 1)), tmp1);
                end
                ind_z_mat           = 0;
                ind_z_flat(ind_z_flat < 1) ...
                                    = 1; % limit too-low indices
                ind_z_flat(ind_z_flat > num_sample_trim) ...
                                    = num_sample_trim; % limit too-high indices
                tmp3                = find(sum(~isnan(ind_z_flat)));
                tmp3                = tmp3(~isnan(ind_z_curr(ii, tmp1(tmp3)))); % indices where both flattening and new layer exist
                
                if isempty(tmp3)
                    status_box.String = ['Layer #' num2str(ii) ' has no overlap...'];
                    tmp8			= [tmp8 ii]; %#ok<AGROW> 
					continue
                end
                
                % index of current layer at reference trace for all overlapping traces
                tmp4                = NaN(1, length(tmp1));
				for jj = tmp3
                    [~, tmp2]       = unique(ind_z_flat(:, jj), 'last');
                    tmp2            = intersect((1 + find(diff(ind_z_flat(:, jj)) > 0)), tmp2);
                    if (length(tmp2) > 1)
                        tmp4(jj)    = interp1(ind_z_flat(tmp2, jj), tmp2, ind_z_curr(ii, tmp1(jj)), 'linear', NaN);
                    end
				end
				
				if isempty(find(~isnan(tmp4), 1))
                    status_box.String = ['Layer #' num2str(ii) ' has no good overlap...'];
                    tmp8			= [tmp8 ii]; %#ok<AGROW> 
					continue
				end
				
                ind_z_pk(ii)        = mean(tmp4, 'omitnan'); % best guess y index at reference trace 
                
                % extract best layers again, now including the new layer
                tmp4                = tmp1(tmp3);
                tmp1                = ind_z_curr(~isnan(ind_z_pk), tmp4);
                tmp2                = ind_z_pk(~isnan(ind_z_pk));
                tmp3                = find(sum(~isnan(tmp1)) > ord_poly);
                tmp4                = tmp4(tmp3);
                
                % new polynomials using additional "depth" for this layer
                if (parallel_check && (length(tmp4) >= 1e3))
                    tmp1            = tmp1(:, tmp3);
                    tmp3            = pk.poly_flat(:, tmp3);
					[tmp6, tmp7]	= deal(cell(1, length(tmp4)));
					for jj = 1:length(tmp4)
						tmp6{jj}	= tmp2(~isnan(tmp1(:, jj)));
						tmp7{jj}	= tmp1(~isnan(tmp1(:, jj)), jj);
					end	
					parfor jj = 1:length(tmp4)
                        tmp3(:, jj) = polyfit(tmp6{jj}, tmp7{jj}, ord_poly)';
					end
                    pk.poly_flat(:, tmp4) ...
                                    = tmp3;
					[tmp6, tmp7]	= deal(0);
                else
                    for jj = 1:length(tmp3)
                        pk.poly_flat(:, tmp4(jj)) ...
                                    = polyfit(tmp2(~isnan(tmp1(:, tmp3(jj)))), tmp1(~isnan(tmp1(:, tmp3(jj))), tmp3(jj)), ord_poly)';
                    end
                end
                
                % smooth polynomials again but not as much each time
				if ((ceil(mean(length_smooth) / mean(diff(1e-3 .* data_cat.dist_lin), 'omitnan')) - kk) > 2)
					pk.poly_flat	= smoothdata(pk.poly_flat, 2, 'rlowess', (ceil(mean(length_smooth) / mean(diff(1e-3 .* data_cat.dist_lin), 'omitnan')) - kk), 'includenan'); % decrease amount of smoothing with each iteration
				end
				pk.poly_flat		= fillmissing(pk.poly_flat, 'spline', 2, 'EndValues', 'extrap', 'MaxGap', ceil(length_smooth / median(diff(1e-3 .* data_cat.dist_lin), 'omitnan'))); % fill small gaps (if any)				
				pk.poly_flat(:, (sum(isnan(pk.poly_flat)) >= 1)) ...
									= NaN; % deal with uneven gap-filling
			end
        end
        
        if parallel_check
            pctRunOnAll warning('on', 'MATLAB:polyfit:RepeatedPointsOrRescale')
            pctRunOnAll warning('on', 'MATLAB:polyfit:PolyNotUnique')
        else
            warning('on', 'MATLAB:polyfit:RepeatedPointsOrRescale')
            warning('on', 'MATLAB:polyfit:PolyNotUnique')
        end
        
        % flattened y indices based on fits to kept layers
        ind_z_mat           = single((1:num_sample_trim)');
        ind_z_mat           = ind_z_mat(:, ones(1, data_cat.num_trace)); % matrix of y indices
        switch ord_poly
            case 2
                ind_z_flat  = ((ind_z_mat .^ 2) .* pk.poly_flat(ones(num_sample_trim, 1), :)) + (ind_z_mat .* (pk.poly_flat((2 .* ones(num_sample_trim, 1)), :))) + pk.poly_flat((3 .* ones(num_sample_trim, 1)), :);
            case 3
                ind_z_flat  = ((ind_z_mat .^ 3) .* pk.poly_flat(ones(num_sample_trim, 1), :)) + ((ind_z_mat .^ 2) .* pk.poly_flat((2 .* ones(num_sample_trim, 1)), :)) + ...
					(ind_z_mat .* (pk.poly_flat((3 .* ones(num_sample_trim, 1)), :))) + pk.poly_flat((4 .* ones(num_sample_trim, 1)), :);
        end
        ind_z_mat           = 0;
        
        % prevent out of bounds
        ind_z_flat((ind_z_flat < 1) | (ind_z_flat > num_sample_trim)) ...
                            = NaN;
        
        % prevent non-unique flattening
		tmp3				= (1:num_sample_trim)';
		tmp4				= diff(ind_z_flat) <= 0;
        for ii = find(sum(tmp4))
            ind_z_flat((1 + find(diff(ind_z_flat(:, ii)) <= 0)), ii) ...
                            = NaN;
            [tmp1, tmp2]    = unique(ind_z_flat(:, ii));
            ind_z_flat(setdiff(tmp3, tmp2(~isnan(tmp1)), 'stable'), ii) ...
                            = deal(NaN); % reduce to unique values
        end
        
        status_box.String = 'Done polynomial fitting. Now flattening radargram...';
        pause(0.1)
        
        % flattened radargram based on layer fits
        amp_flat            = NaN(num_sample_trim, data_cat.num_trace, 'single');
        tmp2                = find(sum(~isnan(ind_z_flat)) > ord_poly);
        if isempty(tmp2)
            status_box.String = 'Flattening canceled because of insufficient constraints.';
            pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
            return
        end
        
        if parallel_check
            tmp1            = amp(:, tmp2);
            tmp3            = ind_z_flat(:, tmp2);
            tmp4            = NaN(num_sample_trim, length(tmp2), 'single');
            parfor ii = 1:length(tmp2)
                tmp4(:, ii) = interp1(tmp1(:, ii), tmp3(:, ii));
            end
            amp_flat(:, tmp2) ...
                            = tmp4;
            [tmp1, tmp3, tmp4] ...
                            = deal(0);
        else
			for ii = tmp2
                amp_flat(:, ii) ...
                            = interp1(amp(:, ii), ind_z_flat(:, ii));
			end
        end
        
        status_box.String = 'Flattening layers...';
        pause(0.1)
        
        % flatten surface pick
        ind_surf_flat       = NaN(1, data_cat.num_trace);
        if surf_avail
            for ii = tmp2
                ind_surf_flat(ii) ...
                            = interp1(ind_z_flat(~isnan(ind_z_flat(:, ii)), ii), find(~isnan(ind_z_flat(:, ii))), ind_surf(ii), 'nearest', 'extrap');
            end
            ind_surf_flat((ind_surf_flat < 1) | (ind_surf_flat > num_sample_trim)) ...
                            = NaN;
        end
        
        % flatten bed pick
        ind_bed_flat        = NaN(1, data_cat.num_trace);
        if bed_avail
            for ii = tmp2
                ind_bed_flat(ii) ...
                            = interp1(ind_z_flat(~isnan(ind_z_flat(:, ii)), ii), find(~isnan(ind_z_flat(:, ii))), ind_bed(ii), 'nearest', 'extrap');
            end
            ind_bed_flat((ind_bed_flat < 1) | (ind_bed_flat > num_sample_trim)) ...
                            = NaN;
        end
        
        flat_done           = true;
        
        % re-flatten layers if any exist
        if (pk_done && pk.num_layer)
            pk_ind_z_flat = NaN(pk.num_layer, data_cat.num_trace);
			for ii = find(sum(~isnan(ind_z_flat)) >= 2)
                pk_ind_z_flat(~isnan(pk_ind_z(:, ii)), ii) ...
                            = interp1(ind_z_flat(~isnan(ind_z_flat(:, ii)), ii), find(~isnan(ind_z_flat(:, ii))), pk_ind_z(~isnan(pk_ind_z(:, ii)), ii), 'nearest');
			end
        end
        
        disp_check(5).Visible = 'on';
		disp_group.SelectedObject = disp_check(5);
        pk_select
        disp_type           = 'flat';
        plot_flat
        status_box.String = 'Flattened radargram and layers.';
        pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
    end

%% Semi-automatic layer picking

    function pk_auto(src, event)
        
        if ~any(strcmp(disp_type, {'twtt' 'depth' 'norm' 'flat'}))
            status_box.String = 'Layers can only be traced in twtt, depth, norm or flat.';
            return
        end
        
        if ~load_done
            status_box.String = 'Data not yet loaded.';
            return
        end
        
        if ~pk.num_layer % initialize for no layers
            [pk_ind_z_flat, pk_ind_z_norm] ...
							= deal([]);
        end
        
        % switch to layer picking callbacks
        tmp3                = pk.num_layer + 1; % used to shorten real/flatten conversion loop after this while picking loop
        [pk_gui.KeyPressFcn, pk_gui.WindowButtonDownFcn] = deal('');
		
        while true
            
            status_box.String = 'Left-click: Pick; D: delete; U: undo; L: cut left; R: cut right; C: cut chunk; M: merge; Q: done...';
            
            % get pick and convert to indices
            [ind_x_pk, ind_z_pk, button] ...
                            = ginput(1);
			ind_x_pk		= interp1((1e-3 .* data_cat.dist_lin), ind_num_trace, ind_x_pk, 'nearest', 'extrap');
			switch disp_type
                case {'twtt' 'norm' 'flat'}
                    ind_z_pk= interp1(data_cat.twtt, 1:num_sample_trim, (1e-6 * ind_z_pk), 'nearest', 'extrap');
                case 'depth'
                    ind_z_pk= interp1(data_cat.twtt, 1:num_sample_trim, (1e-6 * (ind_z_pk + (1e6 * (data_cat.twtt_surf(ind_x_pk) - data_cat.twtt(1))))), 'nearest', 'extrap');
			end
			switch disp_type
				case 'norm'
					if (strcmp(disp_type, 'flat') && isnan(amp_norm(ind_z_pk, ind_x_pk)))
						status_box.String = 'Cannot propagate from empty normalized space.';
						pause(0.5)
						continue
					end
				case 'flat'
					if (strcmp(disp_type, 'flat') && isnan(amp_flat(ind_z_pk, ind_x_pk)))
						status_box.String = 'Cannot propagate from empty flattened space.';
						pause(0.5)
						continue
					end
			end
			
            if (button == 1) % trace layer
                
                pk.num_layer= pk.num_layer + 1;
                curr_layer  = pk.num_layer;
                [pk_ind_z, pk_ind_z_flat, pk_ind_z_norm] ...
                            = deal([pk_ind_z; NaN(1, data_cat.num_trace)], [pk_ind_z_flat; NaN(1, data_cat.num_trace)], [pk_ind_z_norm; NaN(1, data_cat.num_trace)]);
				
                pk_prop
				
				if ((strcmp(disp_type, 'twtt') && (length(find(~isnan(pk_ind_z(curr_layer, :)))) <= num_ind_layer_min)) || (strcmp(disp_type, 'norm') && (length(find(~isnan(pk_ind_z_norm(curr_layer, :)))) <= num_ind_layer_min)) || ...
					(strcmp(disp_type, 'flat') && (length(find(~isnan(pk_ind_z_flat(curr_layer, :)))) <= num_ind_layer_min)))
					delete(p_pk(curr_layer))
					[pk_ind_z, pk.num_layer, curr_layer, p_pk, pk_ind_z_flat, pk_ind_z_norm] ...
							= deal(pk_ind_z(1:(end - 1), :), (pk.num_layer - 1), (curr_layer - 1), p_pk(1:(end - 1)), pk_ind_z_flat(1:(end - 1), :), pk_ind_z_norm(1:(end - 1), :));
					status_box.String = 'Layer too small to preserve.';
					pause(0.5)
					continue
				end
				
				if (pk.num_layer > 1) % check if just-picked layer overlaps with an existing one, if so stop it there
					switch disp_type
						case {'twtt' 'depth'}
							tmp1 = (pk_ind_z((curr_layer .* ones((pk.num_layer + 1), 1)), :) - [ind_surf; pk_ind_z(1:(end - 1), :); ind_bed]) == 0; % find overlap with existing layers
							tmp2 = false;
							[~, tmp4] ...
									= find(tmp1(:, 1:(ind_x_pk - 1))); % any overlap to left of pick
							if ~isempty(tmp4)
								tmp4= max(tmp4); % rightmost overlap of left-hand section
								pk_ind_z(curr_layer, 1:tmp4(1)) ...
									= NaN;
								tmp2= true;
							end
							[~, tmp4] ...
									= find(tmp1(:, (ind_x_pk + 1):end)); % now to right of pick
							if ~isempty(tmp4)
								tmp4= min(tmp4); % leftmost overlap of right-hand section
								pk_ind_z(curr_layer, (ind_x_pk + tmp4(1)):end) ...
									= NaN;
								tmp2= true;
							end
							if tmp2
								switch disp_type
									case 'twtt'
										[p_pk(curr_layer).XData, p_pk(curr_layer).YData] ...
											= deal((1e-3 .* data_cat.dist_lin(~isnan(pk_ind_z(curr_layer, :)))), (1e6 .* data_cat.twtt(pk_ind_z(curr_layer, ~isnan(pk_ind_z(curr_layer, :))))'));
									case 'depth'
                                		tmp1= find(~isnan(pk_ind_z(curr_layer, :)) & ~isnan(ind_surf));
                                		tmp2= pk_ind_z(curr_layer, tmp1) - ind_surf(tmp1) + 1;
                                		tmp2((tmp2 < 1) | (tmp2 > num_sample_trim)) ...
                                    		= NaN;
										if ~isempty(find(~isnan(tmp2), 1))
											[p_pk(curr_layer).XData, p_pk(curr_layer).YData] = deal((1e-3 .* data_cat.dist_lin(tmp1(~isnan(tmp2)))), (1e6 .* data_cat.twtt(tmp2(~isnan(tmp2)))'));
										else
											[p_pk(curr_layer).XData, p_pk(curr_layer).YData] = deal(NaN);
										end
								end
							end
						case 'norm'
							tmp1 = (pk_ind_z_norm((pk.num_layer .* ones((pk.num_layer + 1), 1)), :) - [ind_surf_norm; pk_ind_z_norm(1:(pk.num_layer - 1), :); ind_bed_norm]) == 0;
							tmp2 = false;
							[~, tmp4] ...
								= find(tmp1(:, 1:(ind_x_pk - 1)));
							if ~isempty(tmp4)
								tmp4= max(tmp4);
								pk_ind_z_norm(curr_layer, 1:tmp4(1)) ...
									= NaN;
								tmp2= true;
							end
							[~, tmp4] ...
								= find(tmp1(:, (ind_x_pk + 1):end));
							if ~isempty(tmp4)
								tmp4= min(tmp4); % leftmost overlap of right-hand section
								pk_ind_z_norm(curr_layer, (ind_x_pk + tmp4(1)):end) ...
									= NaN;
								tmp2= true;
							end
							if tmp2
								[p_pk(curr_layer).XData, p_pk(curr_layer).YData] = deal((1e-3 .* data_cat.dist_lin(~isnan(pk_ind_z_norm(curr_layer, :)))), ...
																						(1e6 .* data_cat.twtt(pk_ind_z_norm(curr_layer, ~isnan(pk_ind_z_norm(curr_layer, :))))'));
							end
						case 'flat'
							tmp1 = (pk_ind_z_flat((pk.num_layer .* ones((pk.num_layer + 1), 1)), :) - [ind_surf_flat; pk_ind_z_flat(1:(pk.num_layer - 1), :); ind_bed_flat]) == 0;
							tmp2 = false;
							[~, tmp4] ...
								= find(tmp1(:, 1:(ind_x_pk - 1)));
							if ~isempty(tmp4)
								tmp4= max(tmp4);
								pk_ind_z_flat(curr_layer, 1:tmp4(1)) ...
									= NaN;
								tmp2= true;
							end
							[~, tmp4] ...
								= find(tmp1(:, (ind_x_pk + 1):end));
							if ~isempty(tmp4)
								tmp4= min(tmp4); % leftmost overlap of right-hand section
								pk_ind_z_flat(curr_layer, (ind_x_pk + tmp4(1)):end) ...
									= NaN;
								tmp2= true;
							end
							if tmp2
								[p_pk(curr_layer).XData, p_pk(curr_layer).YData] ...
									= deal(1e-3 .* data_cat.dist_lin(~isnan(pk_ind_z_flat(curr_layer, :))), (1e6 .* data_cat.twtt(pk_ind_z_flat(curr_layer, ~isnan(pk_ind_z_flat(curr_layer, :))))'));
							end
					end
				end
                status_box.String = ['Layer #' num2str(curr_layer) ' picked.'];
                
            elseif any(strcmpi(char(button), {'D' 'backspace'})) % delete layer
                
				if ~pk.num_layer
                    continue
				end
				switch disp_type
                    case 'twtt'
						tmp1= NaN((pk.num_layer - tmp3 + 1), 1);
						for ii = 1:(pk.num_layer - tmp3 + 1)
							tmp1(ii) = pk_ind_z((tmp3 + ii - 1), ind_x_pk); % y index at x index pick for each layer
						end
                    case 'depth'
						tmp1= NaN((pk.num_layer - tmp3 + 1), 1);
						for ii = 1:(pk.num_layer - tmp3 + 1)
							tmp1(ii) = pk_ind_z((tmp3 + ii - 1), ind_x_pk) - ind_surf(ind_x_pk) + 1;
						end
					case 'norm'
                        tmp1= pk_ind_z_norm(tmp3:end, ind_x_pk);
                    case 'flat'
                        tmp1= pk_ind_z_flat(tmp3:end, ind_x_pk);
				end
                [tmp2, tmp4]= unique(tmp1);
				if (length(find(~isnan(tmp2))) > 1)
                    tmp1	= interp1(tmp2(~isnan(tmp2)), tmp4(~isnan(tmp2)), ind_z_pk, 'nearest', 'extrap') + tmp3 - 1;
                elseif isscalar(find(~isnan(tmp2)))
                    tmp1	= find(~isnan(tmp2)) + tmp3 - 1;
                else
                    status_box.String = 'Cannot determine which layer to delete. Pick a more distinct x index.';
                    pause(0.1)
                    continue
				end
				delete(p_pk(tmp1))
                tmp2        = setdiff(1:pk.num_layer, tmp1);
				p_pk		= p_pk(tmp2);
                [pk_ind_z, pk.num_layer, pk_ind_z_flat, pk_ind_z_norm] ...
                            = deal(pk_ind_z(tmp2, :), (pk.num_layer - 1), pk_ind_z_flat(tmp2, :), pk_ind_z_norm(tmp2, :));
                status_box.String = ['Deleted layer #' num2str(tmp1) '.'];
                % pause(0.1)
                
            elseif strcmpi(char(button), 'U') % undo
                
				if (pk.num_layer < tmp3)
                    continue
				end
                delete(p_pk(end))
                [pk_ind_z, p_pk, pk.num_layer, pk_ind_z_flat, pk_ind_z_norm] ...
                            = deal(pk_ind_z(1:(end - 1), :), p_pk(1:(end - 1)), (pk.num_layer - 1), pk_ind_z_flat(1:(end - 1), :), pk_ind_z_norm(1:(end - 1), :));
                status_box.String = 'Undid last layer.';
                
            elseif (strcmpi(char(button), 'L') || strcmpi(char(button), 'R') || strcmpi(char(button), 'C')) % delete portion of layer
                
				% no new layers yet
				if (pk.num_layer < tmp3)
                    continue
				end
				
				try
					switch disp_type
						case {'twtt' 'depth'}
							tmp1 = pk_ind_z(tmp3:end, ind_x_pk);
						case 'norm'
							tmp1 = pk_ind_z_norm(tmp3:end, ind_x_pk);
						case 'flat'
							tmp1 = pk_ind_z_flat(tmp3:end, ind_x_pk);
					end
				catch
					status_box.String = 'Something went wrong compiling new layer z-indices. Try again.';
					pause(0.5)
					continue
				end
				tmp2        = find(~isnan(tmp1));
				
				if ((length(tmp2) > 1) && (length(unique(tmp1(tmp2))) == length(tmp2)))
                    tmp1    = interp1(tmp1(tmp2), tmp2, ind_z_pk, 'nearest', 'extrap') + tmp3 - 1;
                elseif isscalar(tmp2)
                    tmp1    = tmp3 - 1 + tmp2;
                else
                    status_box.String = 'Cannot determine which new layer to edit. Pick a more distinct x index.';
                    pause(0.5)
                    continue
				end
				
				if strcmpi(char(button), 'C')
                    status_box.String = 'Now choose right end of cut...';
					pause(0.1)
                    [ind_x_pk(2), ~] ...
                            = ginput(1);
					ind_x_pk(2)	...
							= interp1((1e-3 .* data_cat.dist_lin), ind_num_trace, ind_x_pk(2), 'nearest', 'extrap');
				end
				
				switch char(button)
                    case {'l' 'L'}
                        tmp2= 1:ind_x_pk;
                    case {'r' 'R'}
                        tmp2 = ind_x_pk:data_cat.num_trace;
                    case {'c' 'C'}
                        tmp2= ind_x_pk(1):ind_x_pk(2);
				end
				
				switch disp_type
                    case {'twtt' 'depth'}
                        pk_ind_z(tmp1, tmp2) ...
                            = NaN;
					case 'norm'
                        pk_ind_z_norm(tmp1, tmp2) ...
                            = NaN;
                    case 'flat'
                        pk_ind_z_flat(tmp1, tmp2) ...
                            = NaN;
				end
                switch disp_type
                    case {'twtt' 'depth'}
						if isempty(find(~isnan(pk_ind_z(tmp1, :)), 1))
                            tmp2 = true;
                        else
                            tmp2 = false;
						end
					case 'norm'
						if isempty(find(~isnan(pk_ind_z_norm(tmp1, :)), 1))
                            tmp2 = true;
                        else
                            tmp2 = false;
						end
                    case 'flat'
                        if isempty(find(~isnan(pk_ind_z_flat(tmp1, :)), 1))
                            tmp2 = true;
                        else
                            tmp2 = false;
                        end
                end
                if tmp2
                    tmp4    = setdiff(1:pk.num_layer, tmp1);
                    [pk_ind_z, p_pk, pk.num_layer, pk_ind_z_norm, pk_ind_z_flat] ...
                            = deal(pk_ind_z(tmp4, :), p_pk(tmp4), (pk.num_layer - 1), pk_ind_z_norm(tmp4, :), pk_ind_z_flat(tmp4, :));
                    status_box.String = 'Deleted edited layer because it is now empty.';
                    pause(0.1)
                else
                    switch disp_type
                        case 'twtt'
                            [p_pk(tmp1).XData, p_pk(tmp1).YData] = deal((1e-3 .* data_cat.dist_lin(~isnan(pk_ind_z(tmp1, :)))), (1e6 .* data_cat.twtt(pk_ind_z(tmp1, ~isnan(pk_ind_z(tmp1, :))))'));
                        case 'depth'
                            tmp2 = find(~isnan(pk_ind_z(tmp1, :)) & ~isnan(ind_surf));
                            tmp4 = pk_ind_z(tmp1, tmp2) - ind_surf(tmp2) + 1;
                            tmp4((tmp4 < 1) | (tmp4 > num_sample_trim)) ...
                                 = NaN;
                            [p_pk(tmp1).XData, p_pk(tmp1).YData] = deal((1e-3 .* data_cat.dist_lin(tmp2(~isnan(tmp4)))), (1e6 .* data_cat.twtt(tmp4(~isnan(tmp4)))'));
						case 'norm'
                            [p_pk(tmp1).XData, p_pk(tmp1).YData] = deal((1e-3 .* data_cat.dist_lin(~isnan(pk_ind_z_norm(tmp1, :)))), (1e6 .* data_cat.twtt(pk_ind_z_norm(tmp1, ~isnan(pk_ind_z_norm(tmp1, :))))'));
                        case 'flat'
                            [p_pk(tmp1).XData, p_pk(tmp1).YData] = deal((1e-3 .* data_cat.dist_lin(~isnan(pk_ind_z_flat(tmp1, :)))), (1e6 .* data_cat.twtt(pk_ind_z_flat(tmp1, ~isnan(pk_ind_z_flat(tmp1, :))))'));
                    end
                end
                
			elseif strcmpi(char(button), 'M') % merge two new layers (only new layers that have not yet been finalized)
                
				if (pk.num_layer < (tmp3 + 1))
                    status_box.String = 'Not enough new layers to merge.';
                    continue
				end
				
				switch disp_type
					case {'twtt' 'depth'}
						tmp1= NaN((pk.num_layer - tmp3 + 1), 1);
						for ii = 1:(pk.num_layer - tmp3 + 1)
							tmp1(ii) = pk_ind_z((tmp3 + ii - 1), ind_x_pk);
						end
					case 'norm'
						tmp1= pk_ind_z_norm(tmp3:end, ind_x_pk);
					case 'flat'
						tmp1= pk_ind_z_flat(tmp3:end, ind_x_pk);
				end
				
                tmp2        = find(~isnan(tmp1));
                if ((length(tmp2) > 1) && (length(unique(tmp1(tmp2))) == length(tmp2)))
                    tmp1    = interp1(tmp1(tmp2), tmp2, ind_z_pk, 'nearest', 'extrap') + tmp3 - 1;
                elseif isscalar(tmp2)
                    tmp1    = tmp2 + tmp3 - 1;
                else
                    status_box.String = 'Cannot determine which new layer to edit. Pick a more distinct x index.';
                    pause(0.5)
                    continue
                end
                
				p_pk(tmp1).Color = 'y';
                
                status_box.String = 'Now choose layer to merge with (Left-click to select or Q: quit ONLY)...';
                
				[ind_x_pk, ind_z_pk, button] ...
                            = ginput(1);
				ind_x_pk	= interp1((1e-3 .* data_cat.dist_lin), ind_num_trace, ind_x_pk, 'nearest', 'extrap');
				if strcmpi(char(button), 'Q')
					p_pk(tmp1).Color = [1 0.7 0.7];
					continue
				elseif (button ~= 1)
					p_pk(tmp1).Color = [1 0.7 0.7];
					continue
				end
				
                switch disp_type
                    case {'twtt' 'flat' 'norm'}
                        ind_z_pk ...
                            = interp1(data_cat.twtt, 1:num_sample_trim, (1e-6 * ind_z_pk), 'nearest', 'extrap');
                    case 'depth'
                        ind_z_pk ...
                            = interp1(data_cat.twtt, 1:num_sample_trim, (1e-6 * (ind_z_pk + (1e6 * (data_cat.twtt_surf(ind_x_pk) - data_cat.twtt(1))))), 'nearest', 'extrap');
                end
				switch disp_type
					case {'twtt' 'depth'}
						tmp2= NaN((pk.num_layer - tmp3 + 1), 1);
						for ii = 1:(pk.num_layer - tmp3 + 1)
							tmp2(ii) ...
								= pk_ind_z((tmp3 + ii - 1), ind_x_pk);
						end
					case 'norm'
						tmp2= pk_ind_z_norm(tmp3:end, ind_x_pk);
					case 'flat'
						tmp2= pk_ind_z_flat(tmp3:end, ind_x_pk);
				end
                tmp4        = find(~isnan(tmp2));
                if ((length(tmp4) > 1) && (length(unique(tmp2(tmp4))) == length(tmp4)))
                    tmp2    = interp1(tmp2(tmp4), tmp4, ind_z_pk, 'nearest', 'extrap') + tmp3 - 1;
                elseif isscalar(tmp4)
                    tmp2    = tmp4 + tmp3 - 1;
                else
                    status_box.String = 'Cannot determine which layer to merge with. Pick a more distinct x index.';
					p_pk(tmp1).Color = [1 0.7 0.7];
                    pause(0.5)
                    continue
                end
				
				p_pk(tmp2).Color = 'y';
                pause(0.1)
                
                tmp4        = setdiff(1:pk.num_layer, tmp2');
				
				delete(p_pk(tmp2))
                switch disp_type
                    case {'twtt' 'depth'}
                        pk_ind_z(tmp1, isnan(pk_ind_z(tmp1, :))) ...
                            = pk_ind_z(tmp2, isnan(pk_ind_z(tmp1, :)));
						switch disp_type
							case 'twtt'
								[p_pk(tmp1).XData, p_pk(tmp1).YData] = deal((1e-3 .* data_cat.dist_lin(~isnan(pk_ind_z(tmp1, :)))), (1e6 .* data_cat.twtt(pk_ind_z(tmp1, ~isnan(pk_ind_z(tmp1, :))))'));
							case 'depth'
								tmp2	= find(~isnan(pk_ind_z(tmp1, :)) & ~isnan(ind_surf));
								tmp5	= pk_ind_z(tmp1, tmp2) - ind_surf(tmp2) + 1;
								tmp5((tmp5 < 1) | (tmp5 > num_sample_trim)) ...
										= NaN;
								[p_pk(tmp1).XData, p_pk(tmp1).YData] = deal((1e-3 .* data_cat.dist_lin(tmp2(~isnan(tmp5)))), (1e6 .* data_cat.twtt(tmp5(~isnan(tmp5))')));
						end
					case 'norm'
                        pk_ind_z_norm(tmp1, isnan(pk_ind_z_norm(tmp1, :))) ...
                            = pk_ind_z_norm(tmp2, isnan(pk_ind_z_norm(tmp1, :)));
                        [p_pk(tmp1).XData, p_pk(tmp1).YData] = deal((1e-3 .* data_cat.dist_lin(~isnan(pk_ind_z_norm(tmp1, :)))), (1e6 .* data_cat.twtt(pk_ind_z_norm(tmp1, ~isnan(pk_ind_z_norm(tmp1, :))))'));
                    case 'flat'
                        pk_ind_z_flat(tmp1, isnan(pk_ind_z_flat(tmp1, :))) ...
                            = pk_ind_z_flat(tmp2, isnan(pk_ind_z_flat(tmp1, :)));
                        [p_pk(tmp1).XData, p_pk(tmp1).YData] = deal((1e-3 .* data_cat.dist_lin(~isnan(pk_ind_z_flat(tmp1, :)))), (1e6 .* data_cat.twtt(pk_ind_z_flat(tmp1, ~isnan(pk_ind_z_flat(tmp1, :))))'));
                end
                
                [pk_ind_z, p_pk, pk.num_layer, pk_ind_z_flat, pk_ind_z_norm] ...
                            = deal(pk_ind_z(tmp4, :), p_pk(tmp4), (pk.num_layer - 1), pk_ind_z_flat(tmp4, :), pk_ind_z_norm(tmp4, :));
				
				p_pk(tmp4 == tmp1).Color = [1 0.7 0.7];
				
                status_box.String = 'Layers merged.';
                
            elseif strcmpi(char(button), 'E')
                
                reset_xz
                
            elseif strcmpi(char(button), 'W')
                
                pk.num_win  = pk.num_win + 1;
                num_win_edit.String = num2str(pk.num_win);
                status_box.String = ['Vertical search window widened to +/- ' num2str(pk.num_win) ' sample(s).'];
                pause(0.1)
            
            elseif strcmpi(char(button), 'S')
                
                if (pk.num_win > 1)
                    pk.num_win ...
                            = pk.num_win - 1;
                    num_win_edit.String = num2str(pk.num_win);
                    status_box.String = ['Vertical search window narrowed to +/- ' num2str(pk.num_win) ' sample(s).'];
                    pause(0.1)
                end
                
            elseif (button == 28)
                pan_left
            elseif (button == 29)
                pan_right
            elseif (button == 30)
                pan_up
            elseif (button == 31)
                pan_down
            elseif strcmpi(char(button), 'I')
                zoom_in
            elseif strcmpi(char(button), 'O')
                zoom_out
				
			elseif strcmpi(char(button), '3')
				
				if ~pk.num_layer
					continue
				end
				if strcmp(p_pk(1).Visible, 'on')
					[p_pk(:).Visible] = deal('off');
				else
					[p_pk(:).Visible] = deal('on');
				end
                
            elseif strcmpi(char(button), 'Q') % done picking lines
                
                status_box.String = 'Done picking layers, exiting picking mode...';
				pause(0.1)
                break
            end
        end
        
        if ~pk.num_layer
            pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
            status_box.String = 'No layers picked/left.';
            return
        end
        
		switch disp_type
			
			case {'twtt' 'depth'}
				
				if norm_done
					tmp1	= (1:num_sample_trim)'; % normalization vector
					tmp2	= (tmp1(:, ones(1, data_cat.num_trace)) - ind_surf_smooth(ones(num_sample_trim, 1), :)) ./ ind_thick_smooth(ones(num_sample_trim, 1), :); % normalization matrix
					for ii = tmp3:pk.num_layer
						tmp4= find(~isnan(pk_ind_z(ii, :)) & ~isnan(tmp2(1, :)));
						for jj = tmp4
							pk_ind_z_norm(ii, jj) ...
								= interp1(tmp1, tmp2(:, jj), pk_ind_z(ii, jj));
						end
					end
					pk_ind_z_norm(tmp3:pk.num_layer, :) ...
							= round(pk_ind_z_norm(tmp3:pk.num_layer, :) .* num_sample_trim);
					pk_ind_z_norm((pk_ind_z_norm < 1) | (pk_ind_z_norm > num_sample_trim)) ...
							= NaN;
				end
				if flat_done
					for ii = find(sum(~isnan(ind_z_flat)) > ord_poly)
						pk_ind_z_flat(~isnan(pk_ind_z(:, ii)), ii) ...
							= interp1(ind_z_flat(~isnan(ind_z_flat(:, ii)), ii), find(~isnan(ind_z_flat(:, ii))), pk_ind_z(~isnan(pk_ind_z(:, ii)), ii), 'nearest', 'extrap');
					end
					pk_ind_z_flat((pk_ind_z_flat < 1) | (pk_ind_z_flat > num_sample_trim)) ...
							= NaN;
				end
				
			case 'norm'
				
				tmp1		= (1:num_sample_trim)'; % normalization vector
				for ii = tmp3:pk.num_layer
					tmp2	= find(~isnan(pk_ind_z_norm(ii, :)));
					tmp4	= (tmp1(:, ones(1, length(tmp2))) - ind_surf_smooth(ones(num_sample_trim, 1), tmp2)) ./ ind_thick_smooth(ones(num_sample_trim, 1), tmp2); % normalization matrix
					tmp5	= pk_ind_z_norm(ii, tmp2) ./ num_sample_trim;
					for jj = 1:length(tmp2)
						pk_ind_z(ii, tmp2(jj)) ...
							= interp1(tmp4(:, jj), tmp1, tmp5(jj), 'nearest');
					end
				end
				pk_ind_z((pk_ind_z < 1) | (pk_ind_z > num_sample_trim)) ...
							= NaN;
				if flat_done
					for ii = find(sum(~isnan(ind_z_flat)) > ord_poly)
						pk_ind_z_flat(~isnan(pk_ind_z(:, ii)), ii) ...
							= interp1(ind_z_flat(~isnan(ind_z_flat(:, ii)), ii), find(~isnan(ind_z_flat(:, ii))), pk_ind_z(~isnan(pk_ind_z(:, ii)), ii), 'nearest', 'extrap');
					end
				end
				
			case 'flat'
				
                % loop through each and reproject
				for ii = tmp3:pk.num_layer
					pk_ind_z(ii, ~isnan(pk_ind_z_flat(ii, :))) ...
                            = round(ind_z_flat(sub2ind([num_sample_trim data_cat.num_trace], pk_ind_z_flat(ii, ~isnan(pk_ind_z_flat(ii, :))), find(~isnan(pk_ind_z_flat(ii, :))))));
				end
                pk_ind_z((pk_ind_z < 1) | (pk_ind_z > num_sample_trim)) ...
                            = NaN;
				if norm_done
					tmp1	= (1:num_sample_trim)'; % normalization vector
					tmp2	= (tmp1(:, ones(1, data_cat.num_trace)) - ind_surf_smooth(ones(num_sample_trim, 1), :)) ./ ind_thick_smooth(ones(num_sample_trim, 1), :); % normalization matrix
					for ii = tmp3:pk.num_layer
						for jj = find(~isnan(pk_ind_z(ii, :)) & ~isnan(tmp2(1, :)))
							pk_ind_z_norm(ii, jj) ...
								= interp1(tmp1, tmp2(:, jj), pk_ind_z(ii, jj));
						end
					end
					pk_ind_z_norm(tmp3:pk.num_layer, :) ...
							= round(pk_ind_z_norm(tmp3:pk.num_layer, :) .* num_sample_trim);
					pk_ind_z_norm((pk_ind_z_norm < 1) | (pk_ind_z_norm > num_sample_trim)) ...
							= NaN;
				end
		end
		
        pk_check			= true;
        update_pk_color
		
        % make colors bolder
		for ii = tmp3:pk.num_layer
			p_pk(ii).Color = pk_color(ii, :);
		end
        curr_layer      = pk.num_layer; % just starting out
        if (isempty(curr_layer) || any(curr_layer < 1) || any(curr_layer > pk.num_layer))
            curr_layer      = 1;
        end
        
        [layer_list.String, layer_list.Value] = deal([num2cell(1:pk.num_layer) 'surface' 'bed'], curr_layer);
        
		pk_done				= true;
        pk_select
        pk_cross
        pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
    end

%% Propagate layer from pick

    function pk_prop(src, event)
		
        switch disp_type
            case 'twtt'
                [~, pk_ind_z(curr_layer, ind_x_pk)] ...
                            = max(amp((ind_z_pk - pk.num_win):(ind_z_pk + pk.num_win), ind_x_pk)); % y index of nearest max
                pk_ind_z(curr_layer, ind_x_pk) ...
                            = ind_z_pk - ((pk.num_win + 1) - pk_ind_z(curr_layer, ind_x_pk)); % correct y index of max because search was done in a narrow window
            case 'depth'
                tmp1        = NaN(1, data_cat.num_trace);
                tmp2        = [ind_x_pk (ind_z_pk - ind_surf(ind_x_pk) + 1)];
                [~, tmp1(tmp2(1))] ...
                            = max(amp_depth((tmp2(2) - pk.num_win):(tmp2(2) + pk.num_win), tmp2(1)));
                tmp1(tmp2(1)) ...
                            = tmp2(2) - ((pk.num_win + 1) - tmp1(tmp2(1)));
			case 'norm'
                [~, pk_ind_z_norm(curr_layer, ind_x_pk)] ...
							= max(amp_norm((ind_z_pk - pk.num_win):(ind_z_pk + pk.num_win), ind_x_pk));
                pk_ind_z_norm(curr_layer, ind_x_pk) ...
							= ind_z_pk - ((pk.num_win + 1) - pk_ind_z_norm(curr_layer, ind_x_pk));
            case 'flat'
                [~, pk_ind_z_flat(curr_layer, ind_x_pk)] ...
							= max(amp_flat((ind_z_pk - pk.num_win):(ind_z_pk + pk.num_win), ind_x_pk));
                pk_ind_z_flat(curr_layer, ind_x_pk) ...
							= ind_z_pk - ((pk.num_win + 1) - pk_ind_z_flat(curr_layer, ind_x_pk));
        end
        
		% decide how far to go to left
		if length_pk_max_check.Value
			switch disp_type
				case 'twtt'
					tmp4	= (ind_x_pk - 1):-1:interp1((1e-3 .* data_cat.dist_lin), ind_num_trace, ((1e-3 * data_cat.dist_lin(ind_x_pk)) - length_pk_max), 'nearest', 'extrap');
				case 'depth'
					tmp4	= (tmp2(1) - 1):-1:interp1((1e-3 .* data_cat.dist_lin), ind_num_trace, ((1e-3 * data_cat.dist_lin(ind_x_pk)) - length_pk_max), 'nearest', 'extrap');
				case 'norm'
					tmp4	= (ind_x_pk(1) - 1):-1:interp1((1e-3 .* data_cat.dist_lin), ind_num_trace, ((1e-3 * data_cat.dist_lin(ind_x_pk)) - length_pk_max), 'nearest', 'extrap');
					if ~isempty(find(sum(isnan(amp_norm(:, tmp4))), 1)) % stop at normalization gaps
						tmp4 ...
							= tmp4(1:(find(sum(isnan(amp_norm(:, tmp4))), 1) - 1));
					end
				case 'flat'
					tmp4	= (ind_x_pk - 1):-1:interp1((1e-3 .* data_cat.dist_lin), ind_num_trace, ((1e-3 * data_cat.dist_lin(ind_x_pk)) - length_pk_max), 'nearest', 'extrap');
					if ~isempty(find(isnan(pk.poly_flat(1, tmp4)), 1)) % stop at flattening gaps
						tmp4	= tmp4(1:(find(isnan(pk.poly_flat(1, tmp4)), 1) - 1));
					end
			end
		else
			switch disp_type
				case {'twtt' 'norm' 'flat'}
					tmp4	= (ind_x_pk - 1):-1:1;
				case 'depth'
					tmp4	= (tmp2(1) - 1):-1:interp1((1e-3 .* data_cat.dist_lin), ind_num_trace, (1e-3 * data_cat.dist_lin(ind_x_pk) - length_pk_max), 'nearest', 'extrap');
			end
		end
		
        % loop for left of ind_x_pk
		switch disp_type
            case {'twtt' 'norm' 'flat'}
				try
                    switch disp_type
						case 'twtt'
							for ii = tmp4
                    			[~, pk_ind_z(curr_layer, ii)] ...
                        			= max(amp((pk_ind_z(curr_layer, (ii + 1)) - pk.num_win):(pk_ind_z(curr_layer, (ii + 1)) + pk.num_win), ii));
								pk_ind_z(curr_layer, ii) ...
									= pk_ind_z(curr_layer, (ii + 1)) - ((pk.num_win + 1) - pk_ind_z(curr_layer, ii));
							end
						case 'norm'
							for ii = tmp4
								[~, pk_ind_z_norm(curr_layer, ii)] ...
									= max(amp_norm((pk_ind_z_norm(curr_layer, (ii + 1)) - pk.num_win):(pk_ind_z_norm(curr_layer, (ii + 1)) + pk.num_win), ii));
								pk_ind_z_norm(curr_layer, ii) ...
									= pk_ind_z_norm(curr_layer, (ii + 1)) - ((pk.num_win + 1) - pk_ind_z_norm(curr_layer, ii));
							end
                        case 'flat'
							for ii = tmp4
								[~, pk_ind_z_flat(curr_layer, ii)] ...
									= max(amp_flat((pk_ind_z_flat(curr_layer, (ii + 1)) - pk.num_win):(pk_ind_z_flat(curr_layer, (ii + 1)) + pk.num_win), ii));
								pk_ind_z_flat(curr_layer, ii) ...
									= pk_ind_z_flat(curr_layer, (ii + 1)) - ((pk.num_win + 1) - pk_ind_z_flat(curr_layer, ii));
							end
                    end
				catch
					status_box.String = 'Propagation left stopped.';
					pause(0.1)
				end
            case 'depth'
				try
					for ii = tmp4
                        [~, tmp1(ii)] ...
                            = max(amp_depth((tmp1(ii + 1) - pk.num_win):(tmp1(ii + 1) + pk.num_win), ii));
                        tmp1(ii) ...
                            = tmp1(ii + 1) - ((pk.num_win + 1) - tmp1(ii));
                        
					end
				catch
					status_box.String = 'Propagation left stopped.';
					pause(0.1)
				end
		end
		
		% decide how far to go to right
		if length_pk_max_check.Value
			switch disp_type
				case 'twtt'
					tmp4	= (ind_x_pk + 1):interp1((1e-3 .* data_cat.dist_lin), ind_num_trace, (1e-3 * data_cat.dist_lin(ind_x_pk) + length_pk_max), 'nearest', 'extrap');
				case 'depth'
					tmp4	= (tmp2(1) + 1):interp1((1e-3 .* data_cat.dist_lin), ind_num_trace, (1e-3 * data_cat.dist_lin(ind_x_pk) + length_pk_max), 'nearest', 'extrap');
				case 'norm'
					tmp4	= (ind_x_pk(1) + 1):interp1((1e-3 .* data_cat.dist_lin), ind_num_trace, (1e-3 * data_cat.dist_lin(ind_x_pk) + length_pk_max), 'nearest', 'extrap');
					if ~isempty(find(sum(isnan(amp_norm(:, tmp4))), 1)) % stop at normalization gaps
						tmp4= tmp4(1:(find(sum(isnan(amp_norm(:, tmp4))), 1) - 1));
					end
				case 'flat'
					tmp4	= (ind_x_pk + 1):interp1((1e-3 .* data_cat.dist_lin), ind_num_trace, (1e-3 * data_cat.dist_lin(ind_x_pk) + length_pk_max), 'nearest', 'extrap');
					if ~isempty(find(isnan(pk.poly_flat(1, tmp4)), 1))
						tmp4= tmp4(1:(find(isnan(pk.poly_flat(1, tmp4)), 1) - 1));
					end
			end
		else
			switch disp_type
				case {'twtt' 'norm' 'flat'}
					tmp4	= (ind_x_pk + 1):data_cat.num_trace;
				case 'depth'
					tmp4	= (tmp2(1) + 1):data_cat.num_trace;
			end
		end
		
        % loop for right of ind_x_pk
        switch disp_type
            case {'twtt' 'norm' 'flat'}
				try
					switch disp_type
						case 'twtt'
							for ii = tmp4
								[~, pk_ind_z(curr_layer, ii)] ...
									= max(amp((pk_ind_z(curr_layer, (ii - 1)) - pk.num_win):(pk_ind_z(curr_layer, (ii - 1)) + pk.num_win), ii));
								pk_ind_z(curr_layer, ii) ...
									= pk_ind_z(curr_layer, (ii - 1)) - ((pk.num_win + 1) - pk_ind_z(curr_layer, ii));
							end
						case 'norm'
							for ii = tmp4
								[~, pk_ind_z_norm(curr_layer, ii)] ...
									= max(amp_norm((pk_ind_z_norm(curr_layer, (ii - 1)) - pk.num_win):(pk_ind_z_norm(curr_layer, (ii - 1)) + pk.num_win), ii));
								pk_ind_z_norm(curr_layer, ii) ...
									= pk_ind_z_norm(curr_layer, (ii - 1)) - ((pk.num_win + 1) - pk_ind_z_norm(curr_layer, ii));
							end
						case 'flat'
							for ii = tmp4
								[~, pk_ind_z_flat(curr_layer, ii)] ...
									= max(amp_flat((pk_ind_z_flat(curr_layer, (ii - 1)) - pk.num_win):(pk_ind_z_flat(curr_layer, (ii - 1)) + pk.num_win), ii));
								pk_ind_z_flat(curr_layer, ii) ...
									= pk_ind_z_flat(curr_layer, (ii - 1)) - ((pk.num_win + 1) - pk_ind_z_flat(curr_layer, ii));
							end
					end
				catch
					status_box.String = 'Propagation right stopped.';
					pause(0.5)
				end
            case 'depth'
				try
					for ii = tmp4
                        [~, tmp1(ii)] ...
                            = max(amp_depth((tmp1(ii - 1) - pk.num_win):(tmp1(ii - 1) + pk.num_win), ii));
                        tmp1(ii) ...
                            = tmp1(ii - 1) - ((pk.num_win + 1) - tmp1(ii));
					end
				catch
					status_box.String = 'Propagation right stopped.';
					pause(0.5)
				end
        end
		
		if strcmp(disp_type, 'depth')
			pk_ind_z(curr_layer, :) ...
							= tmp1 + ind_surf - 1;
		end
		
        % plot new layer
		switch disp_type
            case 'twtt'
                p_pk(curr_layer) = line(1e-3 .* data_cat.dist_lin(~isnan(pk_ind_z(curr_layer, :))), (1e6 .* data_cat.twtt(pk_ind_z(curr_layer, ~isnan(pk_ind_z(curr_layer, :))))'), ...
									    'LineStyle', 'none', 'Marker', '.', 'Color', [1 0.7 0.7], 'MarkerSize', 8);
            case 'depth'
                tmp1        = find(~isnan(pk_ind_z(curr_layer, :)) & ~isnan(ind_surf));
                tmp2        = pk_ind_z(curr_layer, tmp1) - ind_surf(tmp1) + 1;
                tmp2((tmp2 < 1) | (tmp2 > num_sample_trim)) ...
                            = NaN;
                p_pk(curr_layer) = line(1e-3 .* data_cat.dist_lin(tmp1(~isnan(tmp2))), (1e6 .* data_cat.twtt(tmp2(~isnan(tmp2)))'), 'LineStyle', 'none', 'Marker', '.', 'Color', [1 0.7 0.7], 'MarkerSize', 8);
			case 'norm'
                p_pk(curr_layer) = line(1e-3 .* data_cat.dist_lin(~isnan(pk_ind_z_norm(curr_layer, :))), (1e6 .* data_cat.twtt(pk_ind_z_norm(curr_layer, ~isnan(pk_ind_z_norm(curr_layer, :))))'), ...
										'LineStyle', 'none', 'Marker', '.', 'Color', [1 0.7 0.7], 'MarkerSize', 8);
            case 'flat'
                p_pk(curr_layer) = line(1e-3 .* data_cat.dist_lin(~isnan(pk_ind_z_flat(curr_layer, :))), (1e6 .* data_cat.twtt(pk_ind_z_flat(curr_layer, ~isnan(pk_ind_z_flat(curr_layer, :))))'), ...
										'LineStyle', 'none', 'Marker', '.', 'Color', [1 0.7 0.7], 'MarkerSize', 8);
		end
    end

%% Check for crossing picked layers (that's bad)

    function pk_cross(src, event)
		if ((pk.num_layer < 2) || ~cross_check)
        	return
		end
		tmp2				= find(sum(~isnan(pk_ind_z)) > 1); % traces with >1 layer
		tmp1				= pk_ind_z(:, tmp2);
    	for ii = 1:(pk.num_layer - 1)
        	tmp3			= find(sum(tmp1((ii + 1):end, ~isnan(tmp1(ii, :))), 2, 'omitnan')); % skip this layer if subsequent layers have no overlap where this layer exists
        	if isempty(tmp3)
            	continue
        	end
        	for jj = (ii + tmp3')
            	tmp4		= intersect(find(~isnan(tmp1(ii, :))), find(~isnan(tmp1(jj, :))));
            	if (~isempty(find(diff(sign(tmp1(ii, tmp4) - tmp1(jj, tmp4))), 1)) || all(diff(tmp1([ii jj], tmp4), 1, 1) == 0))
                	[status_box.EdgeColor, status_box.LineWidth, status_box.String] = ...
							deal('r', 3, ['Layer #' num2str(ii) ' crosses layer #' num2str(jj) ' at ~' sprintf('%3.1f', (1e-3 .* data_cat.dist_lin(tmp2(tmp4(find(diff(sign(tmp1(ii, tmp4) - tmp1(jj, tmp4))), 1)))))) ' km.']);
                	cross_pass ...
							= false;
					[layer_list.Value, curr_layer] ...
							= deal(ii);
					return
            	end
        	end
    	end
    	cross_pass       = true;
		[status_box.EdgeColor, status_box.LineWidth] = deal('g', 3); % no crossed layers (that's good)
    end

%% Sort layers from top to bottom based on their mean vertical index

    function pk_sort(src, event)
        [~, tmp1]           = sort(mean(pk_ind_z, 2, 'omitnan'));
        [pk_ind_z, p_pk, pk_ind_z_norm, pk_ind_z_flat] ...
                            = deal(pk_ind_z(tmp1, :), p_pk(tmp1), pk_ind_z_norm(tmp1, :), pk_ind_z_flat(tmp1, :));
		for ii = 1:pk.num_layer
			p_pk(ii).Color = pk_color(ii, :);
		end
    end

%% Choose/highlight the current layer

    function pk_select(src, event)
        
        if ~(pk.num_layer || surf_avail || bed_avail)
            status_box.String = 'No picked layers to select.';
            return
        end
        
        curr_layer          = layer_list.Value;
        
        if isgraphics(p_bed)
            p_bed.MarkerSize = 8;
        end
		if isgraphics(p_surf)
            p_surf.MarkerSize = 8;
		end
		if pk_done
			[p_pk(isgraphics(p_pk)).MarkerSize] = deal(8);
		end
		
		if (curr_layer == (pk.num_layer + 1))
        
            if isgraphics(p_surf)
                p_surf.MarkerSize = 24;
            end
            status_box.String = 'Surface selected.';
            
        elseif (curr_layer == (pk.num_layer + 2))
            
            if isgraphics(p_bed)
                p_bed.MarkerSize = 24;
            end
            status_box.String = 'Bed selected.';
            
        else
            
            if isgraphics(p_pk(curr_layer))
                p_pk(curr_layer).MarkerSize = 24;
            end
            status_box.String = ['Layer #' num2str(curr_layer) ' selected.'];
            
		end
		
    end

%% Choose/select a layer interactively

    function pk_select_gui(src, event)
		
		ind_x_pk    = interp1((1e-3 .* data_cat.dist_lin), ind_num_trace, ind_x_pk, 'nearest', 'extrap');
		
		switch disp_type
			case {'twtt' 'norm' 'clutter' 'flat'}
                ind_z_pk   = interp1(data_cat.twtt, 1:num_sample_trim, (1e-6 * ind_z_pk), 'nearest', 'extrap');
            case 'depth'
                ind_z_pk   = interp1(data_cat.twtt, 1:num_sample_trim, (1e-6 * (ind_z_pk + (1e6 * (data_cat.twtt_surf(ind_x_pk) - data_cat.twtt(1))))), 'nearest', 'extrap');
		end
		
		if pk.num_layer
			switch disp_type
				case {'twtt' 'depth' 'clutter'}
					tmp1	= pk_ind_z(:, ind_x_pk);
				case 'norm'
					tmp1	= pk_ind_z_norm(:, ind_x_pk);
				case 'flat'
					tmp1	= pk_ind_z_flat(:, ind_x_pk);
			end
		else
			tmp1			= [];
		end
		
		switch disp_type
			case {'twtt' 'clutter'}
                if (surf_avail && bed_avail)
                    tmp1    = [tmp1; ind_surf(ind_x_pk); ind_bed(ind_x_pk)];
				elseif surf_avail
                    tmp1    = [tmp1; ind_surf(ind_x_pk)];
				elseif bed_avail
					tmp1	= [tmp1; NaN; ind_bed(ind_x_pk)];
                end
            case 'depth'
				if bed_avail
                    tmp1    = [tmp1; NaN; (ind_bed(ind_x_pk) - ind_surf(ind_x_pk) + 1)];
				end
			case 'norm'
				tmp1		= [tmp1; ind_surf_norm(ind_x_pk); ind_bed_norm(ind_x_pk)];
            case 'flat'
                if (surf_avail && bed_avail)
                    tmp1    = [tmp1; ind_surf_flat(ind_x_pk); ind_bed_flat(ind_x_pk)];
				elseif surf_avail
                    tmp1    = [tmp1; ind_surf_flat(ind_x_pk)];
				elseif bed_avail
                    tmp1    = [tmp1; NaN; ind_bed_flat(ind_x_pk)];
                end
		end
		
        tmp2                = find(~isnan(tmp1));
        tmp1                = tmp1(tmp2);
        [tmp1, tmp3]        = unique(tmp1);
        if (length(tmp1) > 1)
            curr_layer      = interp1(tmp1, tmp2(tmp3), ind_z_pk, 'nearest', 'extrap');
        else
            curr_layer      = tmp2;
        end
        if (isempty(curr_layer) || any(curr_layer < 1) || (length(curr_layer) > 1))
            curr_layer      = 1;
        end
        layer_list.Value = curr_layer;
        pk_select
    end

%% Focus on a layer

    function pk_focus(src, event)
        if ~pk_done
            status_box.String = 'No picked layers to focus on.';
            return
        end
		switch disp_type
			case {'twtt' 'clutter'}
                if (curr_layer == (pk.num_layer + 1))
					reset_dist_min
					reset_dist_max
                    ax_radar.YLim = (1e6 .* [min(data_cat.twtt_surf) max(data_cat.twtt_surf)]);
                elseif (curr_layer == (pk.num_layer + 2))
					reset_dist_min
					reset_dist_max
                    ax_radar.YLim = (1e6 .* [min(data_cat.twtt_bed) max(data_cat.twtt_bed)]);
				else
					ax_radar.XLim = 1e-3 .* data_cat.dist_lin([find(~isnan(pk_ind_z(curr_layer, :)), 1) find(~isnan(pk_ind_z(curr_layer, :)), 1, 'last')]);
                    ax_radar.YLim = (1e6 .* data_cat.twtt([min(pk_ind_z(curr_layer, :)) max(pk_ind_z(curr_layer, :))]));
                end
            case 'depth'
				if (curr_layer == (pk.num_layer + 2))
					reset_dist_min
					reset_dist_max
                    ax_radar.YLim = (1e6 .* [min(data_cat.twtt_bed - data_cat.twtt(ind_surf + 1)) max(data_cat.twtt_bed - data_cat.twtt(ind_surf + 1))]);
                elseif (curr_layer <= pk.num_layer)
					ax_radar.XLim = 1e-3 .* data_cat.dist_lin([find(~isnan(pk_ind_z(curr_layer, :)), 1) find(~isnan(pk_ind_z(curr_layer, :)), 1, 'last')]);
                    ax_radar.YLim = (1e6 .* data_cat.twtt([min(pk_ind_z(curr_layer, :) - ind_surf + 1) max(pk_ind_z(curr_layer, :) - ind_surf + 1)]));
				end
			case 'norm'
				if (curr_layer == (pk.num_layer + 1))
					reset_dist_min
					reset_dist_max
                    ax_radar.YLim = (1e6 .* [min(data_cat.twtt(ind_surf_norm(~isnan(ind_surf_norm)))) max(data_cat.twtt(ind_surf_norm(~isnan(ind_surf_norm))))]);
                elseif (curr_layer == (pk.num_layer + 2))
                    ax_radar.YLim = (1e6 .* [min(data_cat.twtt(ind_bed_norm(~isnan(ind_bed_norm)))) max(data_cat.twtt(ind_bed_norm(~isnan(ind_bed_norm))))]);
				else
					ax_radar.XLim = 1e-3 .* data_cat.dist_lin([find(~isnan(pk_ind_z(curr_layer, :)), 1) find(~isnan(pk_ind_z(curr_layer, :)), 1, 'last')]);
                    ax_radar.YLim = (1e6 .* data_cat.twtt([min(pk_ind_z_norm(curr_layer, :)) max(pk_ind_z_norm(curr_layer, :))]));
				end
            case 'flat'
                if (curr_layer == (pk.num_layer + 1))
					reset_dist_min
					reset_dist_max
                    ax_radar.YLim = (1e6 .* [min(data_cat.twtt(ind_surf_flat(~isnan(ind_surf_flat)))) max(data_cat.twtt(~isnan(ind_surf_flat)))]);
                elseif (curr_layer == (pk.num_layer + 2))
					reset_dist_min
					reset_dist_max
                    ax_radar.YLim = (1e6 .* [min(data_cat.twtt(ind_bed_flat(~isnan(ind_bed_flat)))) max(data_cat.twtt(ind_bed_flat(~isnan(ind_bed_flat))))]);
				else
					ax_radar.XLim = 1e-3 .* data_cat.dist_lin([find(~isnan(pk_ind_z_flat(curr_layer, :)), 1) find(~isnan(pk_ind_z_flat(curr_layer, :)), 1, 'last')]);
                    ax_radar.YLim = (1e6 .* data_cat.twtt([min(pk_ind_z_flat(curr_layer, :)) max(pk_ind_z_flat(curr_layer, :))]));
                end
            otherwise
                return
		end
        tmp1                = ax_radar.XLim;
        [tmp1(1), tmp1(2)]  = deal((tmp1(1) - diff(tmp1)), (tmp1(2) + diff(tmp1)));
        if (tmp1(1) < (1e-3 * dist_min_ref))
            tmp1(1)         = 1e-3 * dist_min_ref;
        end
        if (tmp1(2) > (1e-3 * dist_max_ref))
            tmp1(2)         = 1e-3 * dist_max_ref;
        end
        ax_radar.XLim = tmp1;
        [dist_min, dist_max]= deal((1e3 * tmp1(1)), (1e3 * tmp1(2)));
		dist_min			= max([dist_min dist_min_ref]);
		dist_min_slide.Value= max([(1e-3 * dist_min) dist_min_slide.Min]);
		dist_max			= min([dist_max dist_max_ref]);
		dist_max_slide.Value = min([(1e-3 * dist_max) dist_max_slide.Max]);
        dist_min_edit.String = sprintf('%3.1f', (1e-3 * dist_min));
        dist_max_edit.String = sprintf('%3.1f', (1e-3 * dist_max));
        tmp1                = ax_radar.YLim;
        [tmp1(1), tmp1(2)]  = deal((tmp1(1) - diff(tmp1)), (tmp1(2) + diff(tmp1)));
        if (tmp1(1) < (1e6 * twtt_min_ref))
            tmp1(1)         = 1e6 * twtt_min_ref;
        end
        if (tmp1(2) > (1e6 * twtt_max_ref))
            tmp1(2)         = 1e6 * twtt_max_ref;
        end
        ax_radar.YLim = tmp1;
        [twtt_min, twtt_max]= deal((1e-6 * tmp1(1)), (1e-6 * tmp1(2)));
        if ((1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref))) < twtt_min_slide.Min)
            twtt_min_slide.Value = twtt_min_slide.Min;
        elseif ((1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref))) > twtt_min_slide.Max)
            twtt_max_slide.Value = twtt_min_slide.Max;
        else
            twtt_min_slide.Value = (1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref)));
        end
        if ((1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref))) < twtt_max_slide.Min)
            twtt_max_slide.Value = twtt_max_slide.Min;
        elseif ((1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref))) > twtt_max_slide.Max)
            twtt_max_slide.Value = twtt_max_slide.Max;
        else
            twtt_max_slide.Value = (1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref)));
        end
        twtt_min_edit.String = sprintf('%3.1f', (1e6 * twtt_min));
        twtt_max_edit.String = sprintf('%3.1f', (1e6 * twtt_max));
        if (curr_layer == (pk.num_layer + 1))
            status_box.String = 'Focused on surface.';
        elseif (curr_layer == (pk.num_layer + 2))
            status_box.String = 'Focused on bed.';
        else
            status_box.String = ['Focused on layer #' num2str(curr_layer) '.'];
        end
        narrow_cb
    end

%% Switch to previous layer

    function pk_last(src, event)
        if (curr_layer > 1)
            curr_layer      = curr_layer - 1;
            layer_list.Value = curr_layer;
            pk_select
        end
    end

%% Switch to next layer

    function pk_next(src, event)
        if (curr_layer < pk.num_layer)
            curr_layer      = curr_layer + 1;
            layer_list.Value = curr_layer;
            pk_select
        end
    end

%% Delete layer

    function pk_del(src, event)
        if (~pk_done && (pk.num_layer == 0) && ~surf_avail && ~bed_avail)
            status_box.String = 'No picked layers or surface/bed to delete yet.';
            return
        end
        [pk_gui.KeyPressFcn, pk_gui.WindowButtonDownFcn] = deal('');
        status_box.String = 'Delete current layer? Y: yes; A: ALL traced layers (not surface or bed); otherwise: no.';
        waitforbuttonpress
		if strcmpi(pk_gui.CurrentCharacter, 'A')
			status_box.String = 'Confirm deleting ALL layers? Y: yes; otherwise: no.';
			waitforbuttonpress
			if ~strcmpi(pk_gui.CurrentCharacter, 'Y')
				status_box.String = 'Layer deletion canceled.';
				pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
				return
			end
			if pk.num_layer	
				if any(isgraphics(p_pk))
					delete(p_pk(isgraphics(p_pk)))
				end
				[curr_layer, pk_done, pk_color, pk.num_layer, pk_ind_z, p_pk, pk_ind_z_norm, pk_ind_z_flat] ...
                            = deal(NaN, false, NaN, 0, [], gobjects(0), [], []);
				status_box.String = 'All layers deleted.';
        		pause(0.1)
                layer_list.String = {'surface' 'bed'}; layer_list.Value = 1;
        		pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
				return
			end
		elseif ~strcmpi(pk_gui.CurrentCharacter, 'Y')
            status_box.String = 'Layer deletion canceled.';
            pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
            return
		end
        pk_del_breakout
        if (curr_layer == (pk.num_layer + 1))
            status_box.String = 'Surface deleted.';
        elseif (curr_layer == (pk.num_layer + 2))
            status_box.String = 'Bed deleted.';
        else
            status_box.String = ['Layer #' num2str(curr_layer) ' deleted.'];
        end
        pause(0.1)
        pk_select
        pk_cross
		update_pk_color
        pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
    end

    function pk_del_breakout(src, event)
        if (curr_layer == (pk.num_layer + 1))
            if isgraphics(p_surf)
                [p_surf.XData, p_surf.YData] = deal(NaN);
            end
            [ind_surf, ind_surf_flat, ind_surf_smooth, data_cat.twtt_surf] ...
                            = deal(NaN(1, data_cat.num_trace));
            [norm_done, surf_avail] ...
							= deal(false);
        elseif (curr_layer == (pk.num_layer + 2))
            if isgraphics(p_bed)
                [p_bed.XData, p_bed.YData] = deal(NaN);
            end
            [ind_bed, ind_bed_flat, ind_bed_smooth, data_cat.twtt_bed] ...
                            = deal(NaN(1, data_cat.num_trace));
            [bed_avail, norm_done] ...
							= deal(false);
        else
            tmp1            = setdiff(1:pk.num_layer, curr_layer);
            if isgraphics(p_pk(curr_layer))
                [p_pk(curr_layer).XData, p_pk(curr_layer).YData] = deal(NaN);
            end
            [pk.num_layer, pk_ind_z, p_pk, pk_ind_z_norm, pk_ind_z_flat] ...
                            = deal((pk.num_layer - 1), pk_ind_z(tmp1, :), p_pk(tmp1), pk_ind_z_norm(tmp1, :), pk_ind_z_flat(tmp1, :));
            curr_layer      = max([1 (curr_layer - 1)]);
            if pk.num_layer
                layer_list.String = [num2cell(1:pk.num_layer) 'surface' 'bed']; layer_list.Value = curr_layer;
            else
                layer_list.String = {'surface' 'bed'}; layer_list.Value = 1;
            end
        end
    end

%% Adjust/edit current layer

    function pk_adj(src, event)
        
        if ~(pk_done || surf_avail || bed_avail)
            status_box.String = 'No layers to adjust yet.';
            return
        end
        
        tmp5                = 0;
		
        [pk_gui.KeyPressFcn, pk_gui.WindowButtonDownFcn] = deal('');
		
        while true
            
			if isempty(find(strcmp(button, {'L1' 'R1' 'LR2'}), 1))
            	status_box.String ...
							= 'L: delete left; R: delete right; C: cut (left start); U: undo; Q: quit.';
            	% get pick and convert to indices
            	[ind_x_pk, ~, button] ...
                            = ginput(1);
				button		= char(button);
			end
			ind_x_pk		= interp1((1e-3 .* data_cat.dist_lin), ind_num_trace, ind_x_pk, 'nearest', 'extrap');
			
            if ~isempty(find(strcmpi(button, {'L' 'L1' 'R' 'R1' 'C'}), 1))
                
                switch button
                    case {'l' 'L' 'L1'}
                        tmp1= 1:ind_x_pk;
                    case {'r' 'R' 'R1'}
                        tmp1= ind_x_pk:data_cat.num_trace;
                    case {'c' 'C'}
                        status_box.String = 'Now choose right end of cut...';
                        [tmp2, ~] ...
							= ginput(1);
						tmp2= interp1((1e-3 .* data_cat.dist_lin), ind_num_trace, tmp2, 'nearest', 'extrap');
                        tmp4= tmp2; % in case of undo
                        tmp1= ind_x_pk:tmp2;
                end
                
                [tmp2, tmp3]= deal(button, cell(1, 4));
                
				if (curr_layer == (pk.num_layer + 1))
                    [tmp3{1}, tmp3{3}] ...
                            = deal(ind_surf(tmp1), data_cat.twtt_surf(tmp1));
                    [data_cat.twtt_surf(tmp1), ind_surf(tmp1)] ...
                            = deal(NaN);
					if norm_done
                        tmp3{4} ...
                            = ind_surf_norm(tmp1);
                        ind_surf_norm(tmp1) ...
							= NaN;
					end
                    if flat_done
                        tmp3{2} ...
                            = ind_surf_flat(tmp1);
                        ind_surf_flat(tmp1) ...
                            = NaN;
                    end
                elseif (curr_layer == (pk.num_layer + 2))
                    [tmp3{1}, tmp3{3}] ...
                            = deal(ind_bed(tmp1), data_cat.twtt_bed(tmp1));
                    [data_cat.twtt_bed(tmp1), ind_bed(tmp1)] ...
                            = deal(NaN);
					if norm_done
						tmp3{4} ...
							= ind_bed_norm(tmp1);
						ind_bed_norm(tmp1) ...
							= NaN;
					end
                    if flat_done
                        tmp3{2} ...
                            = ind_bed_flat(tmp1);
                        ind_bed_flat(tmp1) ...
                            = NaN;
                    end
                else
                    tmp3{1} = pk_ind_z(curr_layer, tmp1);
					tmp3{2} = pk_ind_z_flat(curr_layer, tmp1);
					tmp3{4} = pk_ind_z_norm(curr_layer, tmp1);
                    [pk_ind_z(curr_layer, tmp1), pk_ind_z_norm(curr_layer, tmp1), pk_ind_z_flat(curr_layer, tmp1)] ...
                            = deal(NaN);
				end
				
                pk_fix
                
                tmp5        = ind_x_pk;
                switch button
                    case {'l' 'L' 'L1'}
                        if (curr_layer == (pk.num_layer + 1))
                            status_box.String = ['Surface cut left at ' num2str((1e-3 * data_cat.dist_lin(ind_x_pk)), '%3.1f') ' km.'];
                        elseif (curr_layer == (pk.num_layer + 2))
                            status_box.String = ['Bed cut left at ' num2str((1e-3 * data_cat.dist_lin(ind_x_pk)), '%3.1f') ' km.'];
                        else
                            status_box.String = ['Layer #' num2str(curr_layer) ' cut left at ' num2str((1e-3 * data_cat.dist_lin(ind_x_pk)), '%3.1f') ' km.'];
                        end
                    case {'r' 'R' 'R1'}
                        if (curr_layer == (pk.num_layer + 1))
                            status_box.String = ['Surface cut right at ' num2str((1e-3 * data_cat.dist_lin(ind_x_pk)), '%3.1f') ' km.'];
                        elseif (curr_layer == (pk.num_layer + 2))
                            status_box.String = ['Bed cut right at ' num2str((1e-3 * data_cat.dist_lin(ind_x_pk)), '%3.1f') ' km.'];
                        else
                            status_box.String = ['Layer #' num2str(curr_layer) ' cut right at ' num2str((1e-3 * data_cat.dist_lin(ind_x_pk)), '%3.1f') ' km.'];
                        end
                    case {'c' 'C'}
                        if (curr_layer == (pk.num_layer + 1))
                            status_box.String = ['Surface cut starting at ' num2str((1e-3 * data_cat.dist_lin(ind_x_pk)), '%3.1f') ' km.'];
                        elseif (curr_layer == (pk.num_layer + 2))
                            status_box.String = ['Bed cut starting at ' num2str((1e-3 * data_cat.dist_lin(ind_x_pk)), '%3.1f') ' km.'];
                        else
                            status_box.String = ['Layer #' num2str(curr_layer) ' cut starting at ' num2str((1e-3 * data_cat.dist_lin(ind_x_pk)), '%3.1f') ' km.'];
                        end
                end
                
				if ~isempty(find(strcmp(button, {'L1' 'R1'}), 1))
					button		= 'LR2';
				end
				
            elseif strcmpi(button, 'U') % undo adjustment done in current set
                
                if ~tmp5
                    continue
                end
                
                if ~isempty(find(strcmpi(tmp2, {'L' 'R' 'C'}), 1))
                    
                    tmp1    = cell(1, 3);
                    switch tmp2
                        case {'l' 'L'}
                            tmp1= 1:tmp5;
                        case {'r' 'R'}
                            tmp1= tmp5:data_cat.num_trace;
                        case {'c' 'C'}
                            tmp1= tmp5:tmp4;
                    end
                    if (curr_layer == (pk.num_layer + 1))
                        [ind_surf(tmp1), data_cat.twtt_surf(tmp1)] ...
                            = deal(tmp3, tmp3{3});
						if norm_done
                            ind_surf_norm(tmp1) ...
                                = tmp3{4};
						end
                        if flat_done
                            ind_surf_flat(tmp1) ...
                                = tmp3{2};
                        end
                    elseif (curr_layer == (pk.num_layer + 2))
                        [ind_bed(tmp1), data_cat.twtt_bed(tmp1)] ...
                            = deal(tmp3{1}, tmp3{3});
						if norm_done
							ind_bed_norm(tmp1) ...
								= tmp3{4};
						end
                        if flat_done
                            ind_bed_flat(tmp1) ...
                                = tmp3{2};
                        end
                    else
                        pk_ind_z(curr_layer, tmp1) ...
                            = tmp3{1};
						pk_ind_z_norm(curr_layer, tmp1) ...
							= tmp3{4};
                        pk_ind_z_flat(curr_layer, tmp1) ...
                            = tmp3{2};
                    end
                end
                
                tmp3        = 0;
                pk_fix
                status_box.String = 'Undid previous adjustment.';
                
			elseif ~isempty(find(strcmpi(button, {'Q' 'LR2'}), 1))
                
                if (curr_layer == (pk.num_layer + 1))
                    status_box.String = 'Done adjusting surface.';
                elseif (curr_layer == (pk.num_layer + 2))
                    status_box.String = 'Done adjusting bed.';
                else
					tmp1	= curr_layer;
                    if ((tmp1 ~= curr_layer) && ~isempty(tmp1))
                        status_box.String = ['Done adjusting layer #' num2str(curr_layer) ' (now layer #' num2str(tmp1) ').'];
                        curr_layer ...
                            = tmp1;
                    else
                        status_box.String = ['Done adjusting layer #' num2str(curr_layer) '.'];
                    end
                    if (isempty(curr_layer) || any(curr_layer < 1))
                        curr_layer ...
                            = 1;
                    end
                    layer_list.Value = curr_layer;
                    pk_select
                end
                pk_cross
                pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
                return
            end
        end
    end

%% Fix layer based on interactive adjustments

    function pk_fix(src, event)
        
		if (curr_layer == (pk.num_layer + 1))
            tmp1            = find(~isnan(ind_surf));
        elseif (curr_layer == (pk.num_layer + 2))
            tmp1            = find(~isnan(ind_bed));
        else
            tmp1            = find(~isnan(pk_ind_z(curr_layer, :)));
		end
		
		if (length(find(~isnan(tmp1))) < num_ind_layer_min)
            pk_del
            status_box.String = 'Layer now too small so it was deleted.';
            return
		end
        
        if (curr_layer == (pk.num_layer + 1))
            
            amp_depth       = NaN(num_sample_trim, data_cat.num_trace, 'single');
			for ii = find(~isnan(ind_surf))
                amp_depth(1:(num_sample_trim - ind_surf(ii) + 1), ii) ...
                            = amp(ind_surf(ii):num_sample_trim, ii); % shift data up to surface
			end
			norm_done		= false;
			switch disp_type
				case 'depth'
					plot_depth
				case 'norm'
					plot_norm
				otherwise
					show_surfbed
					update_pk_plot
			end
            
        elseif (curr_layer == (pk.num_layer + 2))

			norm_done		= false;
			if strcmp(disp_type, 'norm')
				plot_norm
			else
				show_surfbed
			end
            
		else
			if isgraphics(p_pk(curr_layer))
				update_pk_plot
			end
        end
        
        pk_select
        tmp2                = button;
	end

%% Cut all layers (deals with weird spots: 2010-3)

    function pk_trimx(src, event)
        
        if ~pk_done
            status_box.String = 'No layers to trim x values for yet.';
            return
        end
        
        [pk_gui.KeyPressFcn, pk_gui.WindowButtonDownFcn] = deal('');
		
        while true
            
            status_box.String = 'Pick cut (left start); Q: quit.';
            
            % get pick and convert to indices
            [ind_x_pk, ~, button] ...
                            = ginput(1);
			ind_x_pk		= interp1((1e-3 .* data_cat.dist_lin), ind_num_trace, ind_x_pk, 'nearest', 'extrap');
            
            if (button == 1)
                
                status_box.String = 'Now choose right end of cut...';
                [tmp2, ~]	= ginput(1);
				tmp2		= interp1((1e-3 .* data_cat.dist_lin), ind_num_trace, tmp2, 'nearest', 'extrap');
                tmp4		= tmp2; % in case of undo
                tmp1		= ind_x_pk:tmp2;
                [pk_ind_z(:, tmp1), pk_ind_z_norm(:, tmp1), pk_ind_z_flat(:, tmp1)] ...
                            = deal(NaN);
                
                update_pk_plot
				status_box.String = ['All layers cut between ' num2str((1e-3 * data_cat.dist_lin(ind_x_pk)), '%3.1f') ' and ' num2str((1e-3 * data_cat.dist_lin(tmp2)), '%3.1f') ' km.'];
                pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
				
            elseif strcmpi(char(button), 'Q')
                status_box.String = 'Quit trimming.';
                pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
                return
            end
        end
    end

%% Merge two layers

    function pk_merge(src, event)
        
        if ~any(strcmp(disp_type, {'twtt' 'depth' 'norm' 'flat'}))
            status_box.String = 'Must be in twtt, depth, norm or flat to merge layers.';
            return
        end
        if ~pk_done
            status_box.String = 'No picked layers to merge yet.';
            return
        end
        if ((pk.num_layer < 2) && ((pk.num_layer <= 1) && ~(surf_avail || bed_avail)))
            status_box.String = 'Not enough layers to merge.';
            return
        end
		
        [pk_gui.KeyPressFcn, pk_gui.WindowButtonDownFcn] = deal('');
        
        if (curr_layer == (pk.num_layer + 1))
            status_box.String = 'Pick layer to merge with surface (Q: cancel)...';
        elseif (curr_layer == (pk.num_layer + 2))
            status_box.String = 'Pick layer to merge with bed (Q: cancel)...';
        else
            status_box.String = ['Pick layer to merge with layer #' num2str(curr_layer) ' (Q: cancel)...'];
        end
        
        % get pick and convert to indices
        [ind_x_pk, ind_z_pk, button] ...
                            = ginput(1);
        if strcmpi(char(button), 'Q')
            pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
            status_box.String = 'Layer merging canceled.';
            return
        end
        
		ind_x_pk    = interp1((1e-3 .* data_cat.dist_lin), ind_num_trace, ind_x_pk, 'nearest', 'extrap');
        switch disp_type
            case {'twtt' 'norm' 'flat'}
                ind_z_pk    = interp1(data_cat.twtt, 1:num_sample_trim, (1e-6 * ind_z_pk), 'nearest', 'extrap');
            case 'depth'
                ind_z_pk    = interp1(data_cat.twtt, 1:num_sample_trim, (1e-6 * (ind_z_pk + (1e6 * (data_cat.twtt_surf(ind_x_pk) - data_cat.twtt(1))))), 'nearest', 'extrap');
        end
        
        % get current layer positions at ind_x_pk, depending on what layers are available
		switch disp_type
			case {'twtt' 'depth'}
				tmp1		= pk_ind_z(:, ind_x_pk);
			case 'norm'
				tmp1		= pk_ind_z_norm(:, ind_x_pk);
			case 'flat'
				tmp1		= pk_ind_z_flat(:, ind_x_pk);
		end
		
        switch disp_type
            case 'twtt'
                if surf_avail
                    tmp1    = [tmp1; ind_surf(ind_x_pk)];
                else
                    tmp1    = [tmp1; NaN];
                end
                if bed_avail
                    tmp1    = [tmp1; ind_bed(ind_x_pk)];
                end
            case 'depth'
				if bed_avail
                    tmp1    = [tmp1; (ind_bed(ind_x_pk) - ind_surf(ind_x_pk) + 1)];
				end
			case 'norm'
				tmp1		= [tmp1; ind_surf_norm(ind_x_pk); ind_bed_norm(ind_x_pk)];
            case 'flat'
                if surf_avail
                    tmp1    = [tmp1; ind_surf_flat(ind_x_pk)];
                else
                    tmp1    = [tmp1; NaN];
                end
                if bed_avail
                    tmp1    = [tmp1; ind_bed_flat(ind_x_pk)];
                end
        end
        tmp2                = find(~isnan(tmp1));
        tmp1                = tmp1(tmp2);
        [tmp1, tmp3]        = unique(tmp1); % test for repeating tmp1 values
        if (length(tmp1) > 1)
            tmp1            = interp1(tmp1, tmp2(tmp3), ind_z_pk, 'nearest', 'extrap');
        else
            tmp1            = tmp2;
        end
        if (tmp1 == curr_layer)
            pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
            status_box.String = 'Aborted merging because picked layer is the same as current layer.';
            return
		elseif isempty(tmp1)
            pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
            status_box.String = 'No layer picked to be merged. Pick more precisely.';
            return
		elseif (((tmp1 == (pk.num_layer + 1)) && (curr_layer == (pk.num_layer + 2))) || ((tmp1 == (pk.num_layer + 2)) && (curr_layer == (pk.num_layer + 1))))
            pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
            status_box.String = 'Cannot merge surface and bed.';
            return
        end
        
        % force curr_layer to be surface or bed if merging with one of those
        if ((curr_layer <= pk.num_layer) && any(tmp1 == (pk.num_layer + [1 2])))
            [curr_layer, tmp1] ...
                            = deal(tmp1, curr_layer);
        end
        
        if (curr_layer == (pk.num_layer + 1))
            
            % replace NaN values in surface with those in the second layer
            ind_surf(isnan(ind_surf)) ...
                            = pk_ind_z(tmp1, isnan(ind_surf));
            data_cat.twtt_surf(isnan(data_cat.twtt_surf) & ~isnan(pk_ind_z(tmp1, :))) ...
                            = data_cat.twtt(pk_ind_z(tmp1, (isnan(data_cat.twtt_surf) & ~isnan(pk_ind_z(tmp1, :)))))';
            if flat_done
                ind_surf_flat= NaN(1, data_cat.num_trace);
                for ii = find(sum(~isnan(ind_z_flat)) > ord_poly)
                    ind_surf_flat(ii) ...
                            = interp1(ind_z_flat(~isnan(ind_z_flat(:, ii)), ii), find(~isnan(ind_z_flat(:, ii))), ind_surf(ii), 'nearest', 'extrap');
                end
            end
			
            amp_depth       = NaN(num_sample_trim, data_cat.num_trace, 'single');
            for ii = find(~isnan(ind_surf))
                amp_depth(1:(num_sample_trim - ind_surf(ii) + 1), ii) ...
                            = amp(ind_surf(ii):num_sample_trim, ii); % shift data up to surface
            end
            disp_check(2).Visible ...
							= 'on';
            depth_avail     = true;
            
			delete(p_pk(tmp1))
            tmp3            = setdiff(1:pk.num_layer, tmp1);
            [pk_ind_z, p_pk, pk_ind_z_norm, pk_ind_z_flat, pk.num_layer] ...
                            = deal(pk_ind_z(tmp3, :), p_pk(tmp3), pk_ind_z_norm(tmp3, :), pk_ind_z_flat(tmp3, :), (pk.num_layer - 1));
			
            status_box.String = ['Surface and layer #' num2str(tmp1) ' merged.'];
            
            if (curr_layer > tmp1)
                curr_layer  = curr_layer - 1;
            end
            if (isempty(curr_layer) || any(curr_layer < 1))
                curr_layer  = 1;
            end
			if (bed_avail && norm_done)
				do_norm
			end
            if strcmp(disp_type, 'depth')
                plot_depth
            end
            
        elseif (curr_layer == (pk.num_layer + 2))
            
            % replace NaN values in bed with those in the second layer
            ind_bed(isnan(ind_bed)) ...
                            = pk_ind_z(tmp1, isnan(ind_bed));
            data_cat.twtt_bed(isnan(data_cat.twtt_bed) & ~isnan(pk_ind_z(tmp1, :))) ...
                            = data_cat.twtt(pk_ind_z(tmp1, (isnan(data_cat.twtt_bed) & ~isnan(pk_ind_z(tmp1, :)))))';
			
            if flat_done
                ind_bed_flat= NaN(1, data_cat.num_trace);
                for ii = find(sum(isnan(ind_z_flat)) > ord_poly)
                    ind_bed_flat(ii) ...
                            = interp1(ind_z_flat(~isnan(ind_z_flat(:, ii)), ii), find(~isnan(ind_z_flat(:, ii))), ind_bed(ii), 'nearest', 'extrap');
                end
            end
			
			delete(p_pk(tmp1))
            tmp3            = setdiff(1:pk.num_layer, tmp1);
            [pk_ind_z, p_pk, pk_ind_z_norm, pk_ind_z_flat, pk.num_layer] ...
                            = deal(pk_ind_z(tmp3, :), p_pk(tmp3), pk_ind_z_norm(tmp3, :), pk_ind_z_flat(tmp3, :), (pk.num_layer - 1));
            
            status_box.String = ['Bed and layer #' num2str(tmp1) ' merged.'];
            
            if (curr_layer > tmp1)
                curr_layer  = curr_layer - 1;
            end
            if (isempty(curr_layer) || any(curr_layer < 1))
                curr_layer  = 1;
            end
			if norm_done
				do_norm
			end
			
        else
            
            % replace NaN values in first layer with those in the second layer
            pk_ind_z(curr_layer, isnan(pk_ind_z(curr_layer, :))) ...
                            = pk_ind_z(tmp1, isnan(pk_ind_z(curr_layer, :)));
            pk_ind_z_flat(curr_layer, isnan(pk_ind_z_flat(curr_layer, :))) ...
                            = pk_ind_z_flat(tmp1, isnan(pk_ind_z_flat(curr_layer, :)));
			
			if norm_done
				pk_ind_z_norm(curr_layer, isnan(pk_ind_z_norm(curr_layer, :))) ...
							= pk_ind_z_norm(tmp1, isnan(pk_ind_z_norm(curr_layer, :)));
			end
			
			delete(p_pk(tmp1))
            tmp2            = setdiff(1:pk.num_layer, tmp1);
            [pk_ind_z, p_pk, pk_ind_z_flat, pk_ind_z_norm, pk.num_layer] ...
                            = deal(pk_ind_z(tmp2, :), p_pk(tmp2), pk_ind_z_flat(tmp2, :), pk_ind_z_norm(tmp2, :), (pk.num_layer - 1));
            
            status_box.String = ['Layers #' num2str(curr_layer) ' and #' num2str(tmp1) ' merged.'];
            
            if (curr_layer > tmp1)
                curr_layer  = curr_layer - 1;
            end
			if (isempty(curr_layer) || any(curr_layer < 1))
                curr_layer  = 1;
			end
        end
        
		update_pk_color
		p_pk(curr_layer).Color = pk_color(curr_layer, :);
		update_pk_plot
        layer_list.String = [num2cell(1:pk.num_layer) 'surface' 'bed']; layer_list.Value = curr_layer;
        pk_select
        pk_cross
        show_surfbed
        show_pk
        pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
    end

%% Split layer

    function pk_split(source, eventdata)
        
		if ~any(strcmp(disp_type, {'twtt' 'depth' 'norm' 'flat'}))
            status_box.String = 'Must be in twtt, depth, norm or flat to split layers.';
            return
		end
        if ~pk_done
            status_box.String = 'No layers to split yet. Load picks first.';
            return
        end
		if any(layer_list.Value == (pk.num_layer + [1 2]))
            status_box.String = 'Cannot split surface or bed.';
            return
		end
		
		[pk_gui.KeyPressFcn, pk_gui.WindowButtonDownFcn] = deal('');
        
        status_box.String = ['Pick location to split layer #' num2str(curr_layer) ' (Q: cancel)...'];
        pause(0.1)
        
        % get pick and convert to indices
        [ind_x_pk, ~, button] ...
                            = ginput(1);
        
        if strcmpi(char(button), 'Q')
			pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
            status_box.String = 'Layer splitting canceled.';
            return
        end
        
        % get index position of ind_x_pk
		ind_x_pk			= interp1((1e-3 .* data_cat.dist_lin), ind_num_trace, ind_x_pk, 'nearest',  'extrap');
		
        % check split is worthwhile
		switch disp_type
			case {'twtt' 'depth'}
				if (isempty(find(~isnan(pk_ind_z(curr_layer, 1:ind_x_pk)), 1)) || isempty(find(~isnan(pk_ind_z(curr_layer, ind_x_pk:end)), 1)))
					pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
                    status_box.String = 'Aborted splitting because current layer is absent on one side of split point.';
                    return
				end
			case 'norm'
				if (isempty(find(~isnan(pk_ind_z_norm(curr_layer, 1:ind_x_pk)), 1)) || isempty(find(~isnan(pk_ind_z_norm(curr_layer, ind_x_pk:end)), 1)))
					pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
                    status_box.String = 'Aborted splitting because current layer is absent on one side of split point.';
                    return
				end
			case 'flat'
                if (isempty(find(~isnan(pk_ind_z_flat(curr_layer, 1:ind_x_pk)), 1)) || isempty(find(~isnan(pk_ind_z_flat(curr_layer, ind_x_pk:end)), 1)))
					pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
                    status_box.String = 'Aborted splitting because current layer is absent on one side of split point.';
                    return
                end
		end		
        
		tmp1				= line((1e-3 .* data_cat.dist_lin(ind_x_pk([1 1]))), (1e6 .* [twtt_min_ref twtt_max_ref]), 'Color', 'w', 'LineWidth', 3, 'LineStyle', '--');
        status_box.String = 'Split current layer across division? Y: yes; otherwise: no.';
        waitforbuttonpress
		delete(tmp1)
		tmp1				= NaN;
        if ~strcmpi(pk_gui.CurrentCharacter, 'Y')
            status_box.String = 'Layer splitting canceled.';
            pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
            return
        end
		
        pk.num_layer        = pk.num_layer + 1; % add a layer
        update_pk_color
        
		% breakout layers
		pk_ind_z(pk.num_layer, :) ...
							= NaN(1, data_cat.num_trace);
		pk_ind_z(pk.num_layer, ind_x_pk:end) ...
							= pk_ind_z(curr_layer, ind_x_pk:end);
		pk_ind_z(curr_layer, ind_x_pk:end) ...
							= NaN;
		pk_ind_z_norm = [pk_ind_z_norm; NaN(1, data_cat.num_trace)];	
		if norm_done
			pk_ind_z_norm(pk.num_layer, ind_x_pk:end) ...
							= pk_ind_z_norm(curr_layer, ind_x_pk:end);
			pk_ind_z_norm(curr_layer, ind_x_pk:end) ...
							= NaN;
		end
		pk_ind_z_flat = [pk_ind_z_flat; NaN(1, data_cat.num_trace)];
		if flat_done
			pk_ind_z_flat(pk.num_layer, ind_x_pk:end) ...
							= pk_ind_z_flat(curr_layer, ind_x_pk:end);
			pk_ind_z_flat(curr_layer, ind_x_pk:end) ...
							= NaN;
		end
		
		% breakout plots
		switch disp_type
			case 'twtt'
				[p_pk(curr_layer).XData, p_pk(curr_layer).YData] ...
							= deal((1e-3 .* data_cat.dist_lin(~isnan(pk_ind_z(curr_layer, :)))), (1e6 .* data_cat.twtt(pk_ind_z(curr_layer, ~isnan(pk_ind_z(curr_layer, :))))'));
                p_pk(pk.num_layer) = line((1e-3 .* data_cat.dist_lin(~isnan(pk_ind_z(pk.num_layer, :)))), (1e6 .* data_cat.twtt(pk_ind_z(pk.num_layer, ~isnan(pk_ind_z(pk.num_layer, :))))'), ...
									    'LineStyle', 'none', 'Marker', '.', 'Color', pk_color(pk.num_layer, :), 'MarkerSize', 8);
			case 'depth'
                tmp1        = find(~isnan(pk_ind_z(curr_layer, :)) & ~isnan(ind_surf));
                tmp2        = pk_ind_z(curr_layer, tmp1) - ind_surf(tmp1) + 1;
                tmp2((tmp2 < 1) | (tmp2 > num_sample_trim)) ...
                            = NaN;
				[p_pk(curr_layer).XData, p_pk(curr_layer).YData] ...
							= deal((1e-3 .* data_cat.dist_lin(tmp1(~isnan(tmp2)))), (1e6 .* data_cat.twtt(tmp2(~isnan(tmp2)))'));
                tmp1        = find(~isnan(pk_ind_z(pk.num_layer, :)) & ~isnan(ind_surf));
                tmp2        = pk_ind_z(pk.num_layer, tmp1) - ind_surf(tmp1) + 1;
                tmp2((tmp2 < 1) | (tmp2 > num_sample_trim)) ...
                            = NaN;
                p_pk(pk.num_layer) = line((1e-3 .* data_cat.dist_lin(tmp1(~isnan(tmp2)))), (1e6 .* data_cat.twtt(tmp2(~isnan(tmp2)))'), 'LineStyle', 'none', 'Marker', '.', 'Color', pk_color(pk.num_layer, :), 'MarkerSize', 8);
			case 'norm'
				[p_pk(curr_layer).XData, p_pk(curr_layer).YData] ...
							= deal((1e-3 .* data_cat.dist_lin(~isnan(pk_ind_z_norm(curr_layer, :)))), (1e6 .* data_cat.twtt(pk_ind_z_norm(curr_layer, ~isnan(pk_ind_z_norm(curr_layer, :))))'));
                p_pk(pk.num_layer) = line((1e-3 .* data_cat.dist_lin(~isnan(pk_ind_z_norm(pk.num_layer, :)))), (1e6 .* data_cat.twtt(pk_ind_z_norm(pk.num_layer, ~isnan(pk_ind_z_norm(pk.num_layer, :))))'), ...
										'LineStyle', 'none', 'Marker', '.', 'Color', pk_color(pk.num_layer, :), 'MarkerSize', 8);
			case 'flat'
				[p_pk(curr_layer).XData, p_pk(curr_layer).YData] ...
							= deal((1e-3 .* data_cat.dist_lin(~isnan(pk_ind_z_flat(curr_layer, :)))), (1e6 .* data_cat.twtt(pk_ind_z_flat(curr_layer, ~isnan(pk_ind_z_flat(curr_layer, :))))'));
                p_pk(pk.num_layer) = line((1e-3 .* data_cat.dist_lin(~isnan(pk_ind_z_flat(pk.num_layer, :)))), (1e6 .* data_cat.twtt(pk_ind_z_flat(pk.num_layer, ~isnan(pk_ind_z_flat(pk.num_layer, :))))'), ...
										'LineStyle', 'none', 'Marker', '.', 'Color', pk_color(pk.num_layer, :), 'MarkerSize', 8);
		end
        
		% pk_sort
		
        layer_list.String = [num2cell(1:pk.num_layer) 'surface' 'bed'];
		status_box.String = ['Layer #' num2str(curr_layer) ' split.'];
		pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
	end

%% Shift z position of a layer (deals with older imported layers that were offset)

    function pk_shift(source, eventdata)
		
        if ~pk_done
            status_box.String = 'No layers to shift yet. Load picks first.';
            return
        end
		if any(layer_list.Value == (pk.num_layer + [1 2]))
            status_box.String = 'Cannot shift surface or bed.';
            return
		end
		
		[pk_gui.KeyPressFcn, pk_gui.WindowButtonDownFcn] = deal('');
		
        status_box.String = ['Vertically shift layer ' num2str(curr_layer) ' up (U) in the ice column or down (D)? (Q: cancel)...'];
        pause(0.1)
		waitforbuttonpress
		
        if strcmpi(pk_gui.CurrentCharacter, 'Q')
			pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
            status_box.String = 'Layer shifting canceled.';
            return
		elseif strcmpi(pk_gui.CurrentCharacter, 'U')
			tmp1			= 'u';
            status_box.String = 'Layer to be shifted UP...';
		elseif strcmpi(pk_gui.CurrentCharacter, 'D')
			tmp1			= 'd';
            status_box.String = 'Layer to be shifted DOWN...';
        end
		pause(0.5)
		
        status_box.String = ['Shift layer ' num2str(curr_layer) ' by how many samples? (1 to 9; Q or 0: cancel)...'];
        pause(0.1)
		waitforbuttonpress
		
		if any(strcmpi(pk_gui.CurrentCharacter, {'Q' '0'}))
			pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
            status_box.String = 'Layer shifting canceled.';
            return
		elseif ~isnan(str2double(pk_gui.CurrentCharacter))
			tmp2			= str2double(pk_gui.CurrentCharacter);
			if strcmp(tmp1, 'u')
				tmp2		= -tmp2;
			end	
			status_box.String = ['Layer to be shifted by ' num2str(tmp2) ' samples...'];
		end
		
		status_box.String = 'Shift ALL layers? (A: yes; otherwise no)...';
		pause(0.1)
		waitforbuttonpress
		
		if strcmpi(pk_gui.CurrentCharacter, 'A')
			tmp3			= 1:pk.num_layer;
			status_box.String = 'Shifting all layers...';
			pause(0.1)
		else
			tmp3			= curr_layer;
		end
		
		% shift layer(s)
		for ii = tmp3
			pk_ind_z(ii, :) = pk_ind_z(ii, :) + tmp2;
		end
		if norm_done
			tmp1			= (1:num_sample_trim)'; % normalization vector
			tmp2			= (tmp1(:, ones(1, data_cat.num_trace)) - ind_surf_smooth(ones(num_sample_trim, 1), :)) ./ ind_thick_smooth(ones(num_sample_trim, 1), :); % normalization matrix
			for ii = tmp3
				tmp4		= find(~isnan(pk_ind_z(ii, :)) & ~isnan(tmp2(1, :)));
				for jj = tmp4
					pk_ind_z_norm(ii, jj) ...
							= interp1(tmp1, tmp2(:, jj), pk_ind_z(ii, jj));
				end
				pk_ind_z_norm(ii, :) ...
							= round(pk_ind_z_norm(ii, :) .* num_sample_trim);
			end
			pk_ind_z_norm((pk_ind_z_norm < 1) | (pk_ind_z_norm > num_sample_trim)) ...
							= NaN;
		end
		if flat_done
			for ii = find(sum(~isnan(ind_z_flat), 1) > ord_poly)
				pk_ind_z_flat(~isnan(pk_ind_z(:, ii)), ii) ...
							= interp1(ind_z_flat(~isnan(ind_z_flat(:, ii)), ii), find(~isnan(ind_z_flat(:, ii))), pk_ind_z(~isnan(pk_ind_z(:, ii)), ii), 'nearest', 'extrap');
			end
		end
		
		% redo plots
		switch disp_type
			case 'twtt'
				for ii = tmp3
					p_pk(ii).YData ...
							= 1e6 .* data_cat.twtt(pk_ind_z(ii, ~isnan(pk_ind_z(ii, :))))';
				end
			case 'depth'
				for ii = tmp3
                	tmp1    = find(~isnan(pk_ind_z(ii, :)) & ~isnan(ind_surf));
                	tmp2    = pk_ind_z(ii, tmp1) - ind_surf(tmp1) + 1;
                	tmp2((tmp2 < 1) | (tmp2 > num_sample_trim)) ...
                            = NaN;
					[p_pk(ii).XData, p_pk(ii).YData] ...
							= deal((1e-3 .* data_cat.dist_lin(tmp1(~isnan(tmp2)))), (1e6 .* data_cat.twtt(tmp2(~isnan(tmp2)))'));
				end
			case 'norm'
				for ii = tmp3
					[p_pk(ii).XData, p_pk(ii).YData] ...
							= deal((1e-3 .* data_cat.dist_lin(~isnan(pk_ind_z_norm(ii, :)))), (1e6 .* data_cat.twtt(pk_ind_z_norm(ii, ~isnan(pk_ind_z_norm(ii, :))))'));
				end
			case 'flat'
				for ii = tmp3
					[p_pk(ii).XData, p_pk(ii).YData] ...
							= deal((1e-3 .* data_cat.dist_lin(~isnan(pk_ind_z_flat(ii, :)))), (1e6 .* data_cat.twtt(pk_ind_z_flat(ii, ~isnan(pk_ind_z_flat(ii, :))))'));
				end
		end
		
		if (length(tmp3) > 1)
			status_box.String = 'All layers shifted.';
		else
			status_box.String = ['Layer #' num2str(curr_layer) ' shifted.'];
		end
		
		pk_cross
		pk_sort
		
		pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
    end

%% Assign current layer to surface or bed

    function pk_surfbed(src, event)
		
        if ~(pk_done && pk.num_layer && curr_layer)
            status_box.String = 'No picked layers to assign to surface or bed.';
            return
        end
		
        [pk_gui.KeyPressFcn, pk_gui.WindowButtonDownFcn] = deal('');
        status_box.String = 'Assign current layer? S: surface; B: bed; otherwise: no...';
        waitforbuttonpress
		
		if strcmpi(pk_gui.CurrentCharacter, 'S')
			
            status_box.String = 'Assigning layer to surface...';
			
            ind_surf        = pk_ind_z(curr_layer, :);
            data_cat.twtt_surf = NaN(1, data_cat.num_trace);
            data_cat.twtt_surf(~isnan(pk_ind_z(curr_layer, :))) ...
                            = data_cat.twtt(pk_ind_z(curr_layer, ~isnan(pk_ind_z(curr_layer, :))))'; % fix concatenated data traveltime
			
            pk_del_breakout
            [depth_avail, surf_avail] ...
							= deal(true);
            amp_depth		= NaN(num_sample_trim, data_cat.num_trace, 'single');
            for ii = find(~isnan(ind_surf))
                amp_depth(1:(num_sample_trim - ind_surf(ii) + 1), ii) ...
							= amp(ind_surf(ii):num_sample_trim, ii); % shift data up to surface
            end
            disp_check(2).Visible = 'on';
			
			if bed_avail
				do_norm
			end
			if flat_done
				ind_surf_flat ...
							= NaN(1, data_cat.num_trace);
				for ii = find(sum(~isnan(ind_z_flat)) > ord_poly)
					ind_surf_flat(ii) ...
                            = interp1(ind_z_flat(~isnan(ind_z_flat(:, ii)), ii), find(~isnan(ind_z_flat(:, ii))), ind_surf(ii), 'nearest', 'extrap');
				end
				ind_surf_flat((ind_surf_flat < 1) | (ind_surf_flat > num_sample_trim)) ...
                            = NaN;
			end
			
            curr_layer      = pk.num_layer + 1;
            layer_list.String = [num2cell(1:pk.num_layer) 'surface' 'bed']; layer_list.Value = (pk.num_layer + 1);
            show_surfbed
			update_pk_plot
            pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
            status_box.String = 'Layer assigned to surface.';
			
        elseif strcmpi(pk_gui.CurrentCharacter, 'B')
			
            status_box.String = 'Assigning layer to bed...';
			
            ind_bed         = pk_ind_z(curr_layer, :);
            data_cat.twtt_bed  = NaN(1, data_cat.num_trace);
            data_cat.twtt_bed(~isnan(pk_ind_z(curr_layer, :))) ...
                            = data_cat.twtt(pk_ind_z(curr_layer, ~isnan(pk_ind_z(curr_layer, :))))';
            pk_del_breakout
            bed_avail       = true;
            if surf_avail
                depth_avail = true;
				do_norm
            end
			if flat_done
				ind_bed_flat ...
							= NaN(1, data_cat.num_trace);
				for ii = find(sum(~isnan(ind_z_flat)) > ord_poly)
					ind_bed_flat(ii) ...
                            = interp1(ind_z_flat(~isnan(ind_z_flat(:, ii)), ii), find(~isnan(ind_z_flat(:, ii))), ind_bed(ii), 'nearest', 'extrap');
				end
				ind_bed_flat((ind_bed_flat < 1) | (ind_bed_flat > num_sample_trim)) ...
                            = NaN;
			end
			
            layer_list.String = [num2cell(1:pk.num_layer) 'surface' 'bed']; layer_list.Value = (pk.num_layer + 2);
            show_surfbed
			update_pk_plot
            pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
            status_box.String = 'Layer assigned to bed.';
			
        else
            pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
            status_box.String = 'Layer not assigned.';
		end
		
		update_pk_color
    end

%% Save layer picks

    function pk_save(src, event)
        
        [pk_gui.KeyPressFcn, pk_gui.WindowButtonDownFcn] = deal('');
		
        % want everything done before saving
        if ~cross_check
			cross_check		= true;
			pk_cross
			if ~cross_pass
				pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
            	status_box.String = 'Crossing layers must be addressed before saving.';
            	return
			end
        end
        
        reset_xz
        pause(0.1)
        
        if isempty(path_pk)
            if ispc
                if (exist([path_data '..\pk'], 'dir') || exist([path_data '..\..\pk'], 'dir'))
                    path_pk = [path_data(1:strfind(path_data, '\cat')) 'pk\'];
                end
            else
                if (exist([path_data '../pk'], 'dir') || exist([path_data '../../pk'], 'dir'))
                    path_pk = [path_data(1:strfind(path_data, '/cat')) 'pk/'];
                end
            end
        end
        
        tmp1                = path_pk;
        if isempty(file_pk)
            file_pk         = [file_data(1:(end - 4)) '_pk.mat'];
        end
        
        if ~isempty(path_pk)
            [file_pk, path_pk] = uiputfile('*.mat', 'Save picks:', [path_pk file_data(1:(end - 4)) '_pk.mat']);
        elseif ~isempty(path_data)
            [file_pk, path_pk] = uiputfile('*.mat', 'Save picks:', [path_data file_data(1:(end - 4)) '_pk.mat']);
        else
            [file_pk, path_pk] = uiputfile('*.mat', 'Save picks:', [file_data(1:(end - 4)) '_pk.mat']);
        end
        
        if ~ischar(file_pk)
            [file_pk, path_pk] = deal('', tmp2);
			pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
            status_box.String = 'Saving canceled.';
            return
        end
        
        status_box.String = 'Saving picks...';
		pause(0.1)
		
		pk_sort
		
        tmp1                = pk;
        
        % save many variables in pk structure for easy reference independent of data later on
        [pk.lat, pk.lon, pk.x, pk.y, pk.num_sample, pk.num_trace, pk.file_in, pk.file_cat, pk.twtt_min_ref, pk.twtt_max_ref, pk.dist, pk.dist_lin, pk.elev_air, pk.time, pk.elev_surf_corr, pk.twtt_surf, pk.twtt_bed, pk.twtt] ...
                            = deal(data_cat.lat, data_cat.lon, data_cat.x, data_cat.y, data_cat.num_sample, data_cat.num_trace, {data_cat.radframeproc_call(:).file_in}, file_data(1:(end - 4)), twtt_min_ref, twtt_max_ref, ...
								   data_cat.dist, data_cat.dist_lin, data_cat.elev_air, data_cat.time, data_cat.elev_surf_corr, data_cat.twtt_surf, data_cat.twtt_bed, data_cat.twtt);
        
        % get traveltimes and echo intensities from indices, and adjust indices as appropriate assuming trimming has occurred
        pk.elev_surf        = pk.elev_air - (data_cat.twtt_surf .* (speed_vacuum / 2)); % ice-sheet surface elevation (used to calculate layer elevations)
        pk.elev_bed         = pk.elev_surf - ((data_cat.twtt_bed - data_cat.twtt_surf) .* (speed_ice / 2)); % bed elevation
        pk.elev_bed_corr    = pk.elev_surf_corr - ((data_cat.twtt_bed - data_cat.twtt_surf) .* (speed_ice / 2));
        
        [pk.int_surf, pk.int_bed] ...
                            = deal(NaN(1, data_cat.num_trace));
        if surf_avail
            pk.int_surf(~isnan(ind_surf)) ...
                            = amp(sub2ind([num_sample_trim data_cat.num_trace], ind_surf(~isnan(ind_surf)), find(~isnan(ind_surf))));
        end
        if bed_avail
            pk.int_bed(~isnan(ind_bed)) ...
                            = amp(sub2ind([num_sample_trim data_cat.num_trace], ind_bed(~isnan(ind_bed)), find(~isnan(ind_bed))));
        end
        
		for ii = 1:pk.num_layer
            [pk.layer(ii).twtt, pk.layer(ii).int] ...
                            = deal(NaN(1, data_cat.num_trace));
            pk.layer(ii).twtt(~isnan(pk_ind_z(ii, :))) ...
                            = data_cat.twtt(pk_ind_z(ii, ~isnan(pk_ind_z(ii, :))))';
            pk.layer(ii).twtt_ice ...
                            = pk.layer(ii).twtt - data_cat.twtt_surf;
            pk.layer(ii).depth ...
                            = pk.layer(ii).twtt_ice .* (speed_ice / 2);
            pk.layer(ii).elev ...
                            = pk.elev_surf - pk.layer(ii).depth;
            pk.layer(ii).elev_corr ...
                            = pk.elev_surf_corr - pk.layer(ii).depth;
            pk.layer(ii).int(~isnan(pk_ind_z(ii, :))) ...
                            = amp(sub2ind([num_sample_trim data_cat.num_trace], pk_ind_z(ii, ~isnan(pk_ind_z(ii, :))), find(~isnan(pk_ind_z(ii, :)))));
            pk.layer(ii).ind_z ...
							= pk_ind_z(ii, :) + pk.ind_trim_start - 1; % adjust due to trimming
		end
		
        pk                  = orderfields(pk);
        pk.layer            = orderfields(pk.layer);
        
        try
            save([path_pk file_pk], '-v7.3', 'pk')
        catch tmp2
            pk              = tmp1;
			pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
            status_box.String = ['Saving pk file unsuccessful. Try again. Error: ' tmp2.message];
            return
        end
        
        % make a simple figure that also gets saved
        set(0, 'DefaultFigureWindowStyle', 'default')
        pk_fig               = figure('Position', [10 10 1600 1000]);
        imagesc((1e-3 .* data_cat.dist_lin), (1e6 .* data_cat.twtt), amp, [db_min db_max])
		
        colormap(bone)
        if surf_avail
            line((1e-3 .* data_cat.dist_lin), (1e6 .* data_cat.twtt_surf), 'LineStyle', 'none', 'Marker', '.', 'Color', 'm', 'MarkerSize', 8)
        end
        if bed_avail
            line((1e-3 .* data_cat.dist_lin), (1e6 .* data_cat.twtt_bed), 'LineStyle', 'none', 'Marker', '.', 'Color', 'm', 'MarkerSize', 8)
        end
        for ii = 1:pk.num_layer
            line((1e-3 .* data_cat.dist_lin), (1e6 .* pk.layer(ii).twtt), 'LineStyle', 'none', 'Marker', '.', 'Color', pk_color(ii, :), 'MarkerSize', 8)
        end
        set(gca, 'FontSize', 18)
        xlabel('Distance (km)')
        ylabel('Traveltime ({\mu}s)')
        title(file_pk(1:(end - 4)), 'FontWeight', 'bold', 'Interpreter', 'none')
        grid on
        box on
        
        pause(0.1)
		warning('off', 'MATLAB:graphics:axestoolbar:PrintWarning')
        print(pk_fig, '-dpng', [path_pk file_pk(1:(end - 4)) '.png'])
        warning('on', 'MATLAB:graphics:axestoolbar:PrintWarning')
		
        pk                  = tmp1;
        set(0, 'DefaultFigureWindowStyle', 'docked')
        
		pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
        status_box.String = [num2str(pk.num_layer) ' picked layer(s) saved as ' file_pk(1:(end - 4)) '.'];
		close(pk_fig)
    end

%% Update minimum twtt

    function slide_twtt_min(src, event)
        if (((1e6 * twtt_max_ref) - (twtt_min_slide.Value - (1e6 * twtt_min_ref))) < (1e6 * twtt_max))
            if twttfix_check
                tmp1        = twtt_max - twtt_min;
            end
            twtt_min        = ((1e6 * twtt_max_ref) - (twtt_min_slide.Value - (1e6 * twtt_min_ref))) * 1e-6;
            if twttfix_check
                twtt_max    = twtt_min + tmp1;
                if (twtt_max > twtt_max_ref)
                    twtt_max= twtt_max_ref;
                    twtt_min= twtt_max - tmp1;
                    if (twtt_min < twtt_min_ref)
                        twtt_min ...
                            = twtt_min_ref;
                    end
                    if ((1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref))) < twtt_min_slide.Min)
                        twtt_min_slide.Value = twtt_min_slide.Min;
                    elseif ((1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref))) > twtt_min_slide.Max)
                        twtt_max_slide.Value = twtt_min_slide.Max;
                    else
                        twtt_min_slide.Value = (1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref)));
                    end
                end
                twtt_max_edit.String = sprintf('%3.1f', (1e6 * twtt_max));
                if ((1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref))) < twtt_max_slide.Min)
                    twtt_max_slide.Value = twtt_max_slide.Min;
                elseif ((1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref))) > twtt_max_slide.Max)
                    twtt_max_slide.Value = twtt_max_slide.Max;
                else
                    twtt_max_slide.Value = (1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref)));
                end
            end
            twtt_min_edit.String = sprintf('%3.1f', (1e6 * twtt_min));
            update_twtt_range
        else
            if ((1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref))) < twtt_min_slide.Min)
                twtt_min_slide.Value = twtt_min_slide.Min;
            elseif ((1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref))) > twtt_min_slide.Max)
                twtt_max_slide.Value = twtt_min_slide.Max;
            else
                twtt_min_slide.Value = (1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref)));
            end
        end
        twtt_min_slide.Enable = 'off';
        drawnow
        twtt_min_slide.Enable = 'on';
    end

%% Update maximum twtt

    function slide_twtt_max(src, event)
        if (((1e6 * twtt_max_ref) - (twtt_max_slide.Value - (1e6 * twtt_min_ref))) > (1e6 * twtt_min))
            if twttfix_check
                tmp1        = twtt_max - twtt_min;
            end
            twtt_max        = ((1e6 * twtt_max_ref) - (twtt_max_slide.Value - (1e6 * twtt_min_ref))) * 1e-6;
            if twttfix_check
                twtt_min    = twtt_max - tmp1;
                if (twtt_min < twtt_min_ref)
                    twtt_min= twtt_min_ref;
                    twtt_max= twtt_min + tmp1;
                    if (twtt_max > twtt_max_ref)
                        twtt_max ...
                            = twtt_max_ref;
                    end
                    if ((1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref))) < twtt_max_slide.Min)
                        twtt_max_slide.Value = twtt_max_slide.Min;
                    elseif ((1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref))) > twtt_max_slide.Max)
                        twtt_max_slide.Value = twtt_max_slide.Max;
                    else
                        twtt_max_slide.Value = (1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref)));
                    end
                end
                twtt_min_edit.String = sprintf('%3.1f', (1e6 * twtt_min));
                if ((1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref))) < twtt_min_slide.Min)
                    twtt_min_slide.Value = twtt_min_slide.Min;
                elseif ((1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref))) > twtt_min_slide.Max)
                    twtt_max_slide.Value = twtt_min_slide.Max;
                else
                    twtt_min_slide.Value = (1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref)));
                end
            end
            twtt_max_edit.String = sprintf('%3.1f', (1e6 * twtt_max));
            update_twtt_range
        else
            if ((1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref))) < twtt_max_slide.Min)
                twtt_max_slide.Value = twtt_max_slide.Min;
            elseif ((1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref))) > twtt_max_slide.Max)
                twtt_max_slide.Value = twtt_max_slide.Max;
            else
                twtt_max_slide.Value = (1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref)));
            end
        end
        twtt_max_slide.Enable = 'off';
		drawnow
        twtt_max_slide.Enable = 'on';
    end

%% Reset minimum twtt

    function reset_twtt_min(src, event)
        if ((1e6 * twtt_max_ref) < twtt_min_slide.Min)
            twtt_min_slide.Value = twtt_min_slide.Min;
        elseif ((1e6 * twtt_max_ref) > twtt_min_slide.Max)
            twtt_max_slide.Value = twtt_min_slide.Max;
        else
            twtt_min_slide.Value = (1e6 * twtt_max_ref);
        end
        twtt_min_edit.String = sprintf('%3.1f', (1e6 * twtt_min_ref));
        twtt_min              = twtt_min_ref;
        update_twtt_range
    end

%% Reset maximum twtt

    function reset_twtt_max(src, event)
        if ((1e6 * twtt_min_ref) < twtt_max_slide.Min)
            twtt_max_slide.Value = twtt_max_slide.Min;
        elseif ((1e6 * twtt_min_ref) > twtt_max_slide.Max)
            twtt_max_slide.Value = twtt_max_slide.Max;
        else
            twtt_max_slide.Value = (1e6 * twtt_min_ref);
        end
        twtt_max_edit.String = sprintf('%3.1f', (1e6 * twtt_max_ref));
        twtt_max              = twtt_max_ref;
        update_twtt_range
    end

%% Update twtt range

    function update_twtt_range(src, event)
        ax_radar.YLim = 1e6 .* [twtt_min twtt_max];
        narrow_cb
    end

%% Update minimum color slider

    function slide_db_min(src, event)
        tmp1				= [db_min db_max];
        tmp2				= [db_min_ref db_max_ref];
        if (cb_min_slide.Value < tmp1(2))
            if cbfix_check1
                tmp3        = diff(tmp1);
            end
            tmp1(1)         = cb_min_slide.Value;
            if cbfix_check1
                tmp1(2)     = tmp1(1) + tmp3;
                if (tmp1(2) > tmp2(2))
                    tmp1(2) = tmp2(2);
                    tmp1(1) = tmp1(2) - tmp3;
                    if (tmp1(1) < tmp2(1))
                        tmp1(1) ...
                            = tmp2(1);
                    end
                    if (tmp1(1) < cb_min_slide.Min)
                        cb_min_slide.Value = cb_min_slide.Min;
                    else
                        cb_min_slide.Value = tmp1(1);
                    end
                end
                cb_max_edit.String = sprintf('%3.0f', tmp1(2));
                if (tmp1(2) > cb_max_slide.Max)
                    cb_max_slide.Value = cb_max_slide.Max;
                else
                    cb_max_slide.Value = tmp1(2);
                end
            end
            cb_min_edit.String = sprintf('%3.0f', tmp1(1));
        else
            if (tmp1(1) < cb_min_slide.Min)
                cb_min_slide.Value = cb_min_slide.Min;
            else
                cb_min_slide.Value = tmp1(1);
            end
        end
		[db_min, db_max]	= deal(tmp1(1), tmp1(2));
        update_db_range
        cb_min_slide.Enable = 'off';
        drawnow
        cb_min_slide.Enable = 'on';
    end

%% Update maximum color slider

    function slide_db_max(src, event)
        tmp1				= [db_min db_max];
        tmp2				= [db_min_ref db_max_ref];
        if (cb_max_slide.Value > tmp1(1))
            if cbfix_check1
                tmp3        = diff(tmp1);
            end
            tmp1(2)         = cb_max_slide.Value;
            if cbfix_check1
                tmp1(1)     = tmp1(2) - tmp1;
                if (tmp1(1) < tmp2(1))
                    tmp1(1) = tmp2(1);
                    tmp1(2) = tmp1(1) + tmp1;
                    if (tmp1(2) > tmp2(2))
                        tmp1(2) ...
                            = tmp2(2);
                    end
                    if (tmp1(2) > cb_max_slide.Max)
                        cb_max_slide.Value = cb_max_slide.Max;
                    else
                        cb_max_slide.Value = tmp1(2);
                    end
                end
                cb_min_edit.String = sprintf('%3.0f', tmp1(1));
                if (tmp1(1) < cb_min_slide.Min)
                    cb_min_slide.Value = cb_min_slide.Min;
                else
                    cb_min_slide.Value = tmp1(1);
                end
            end
            cb_max_edit.String = sprintf('%3.0f', tmp1(2));
        else
            if (tmp1(2) > cb_max_slide.Max)
                cb_max_slide.Value = cb_max_slide.Max;
            else
                cb_max_slide.Value = tmp1(2);
            end
        end
        [db_min, db_max]	= deal(tmp1(1), tmp1(2));
        update_db_range
        cb_max_slide.Enable = 'off';
        drawnow
        cb_max_slide.Enable = 'on';
    end

%% Reset color sliders

    function reset_db(src, event)
		[tmp1, db_min]		= deal(db_min_ref);
        if (tmp1 < cb_min_slide.Min)
            cb_min_slide.Value = cb_min_slide.Min;
        else
            cb_min_slide.Value = tmp1;
        end
        cb_min_edit.String = num2str(tmp1);
		[tmp1, db_max]		= deal(db_max_ref);
        if (tmp1 > cb_max_slide.Max)
            cb_max_slide.Value = cb_max_slide.Max;
        else
            cb_max_slide.Value = tmp1;
        end
        cb_max_edit.String = num2str(tmp1);	
        update_db_range
    end

%% Update color

    function update_db_range(src, event)
		ax_radar.CLim = [db_min db_max];
    end

%% Update minimum distance

    function slide_dist_min(src, event)
        if ((1e3 * dist_min_slide.Value) < dist_max)
            if distfix_check
                tmp1        = dist_max - dist_min;
            end
            dist_min        = (1e3 * dist_min_slide.Value);
            if distfix_check
                dist_max    = dist_min + tmp1;
                if (dist_max > dist_max_ref)
                    dist_max= dist_max_ref;
                    dist_min= dist_max - tmp1;
                    if (dist_min < dist_min_ref)
                        dist_min ...
                            = dist_min_ref;
                    end
                    if (dist_min < (1e3 * dist_min_slide.Min))
                        dist_min_slide.Value = dist_min_slide.Min;
                    else
                        dist_min_slide.Value = 1e-3 * dist_min;
                    end
                end
                if (dist_max > (1e3 * dist_max_slide.Max))
                    dist_max_slide.Value = dist_max_slide.Max;
                end
                dist_max_slide.Value = 1e-3 * dist_max;
                dist_max_edit.String = sprintf('%3.1f', (1e-3 * dist_max));
            end
            dist_min_edit.String = sprintf('%3.1f', (1e-3 * dist_min));
            update_dist_range
        else
            if (dist_min < (1e3 * dist_min_slide.Min))
                dist_min_slide.Value = dist_min_slide.Min;
            else
                dist_min_slide.Value = 1e-3 * dist_min;
            end
        end
        dist_min_slide.Enable = 'off';
        drawnow
        dist_min_slide.Enable = 'on';
    end

%% Update maximum distance

    function slide_dist_max(src, event)
        if (dist_max_slide.Value > (1e-3 * dist_min))
            if distfix_check
                tmp1        = dist_max - dist_min;
            end
            dist_max        = 1e3 * dist_max_slide.Value;
            if distfix_check
                dist_min    = dist_max - tmp1;
                if (dist_min < dist_min_ref)
                    dist_min= dist_min_ref;
                    dist_max= dist_min + tmp1;
                    if (dist_max > dist_max_ref)
                        dist_max ...
                            = dist_max_ref;
                    end
                    if (dist_max > (1e3 * dist_max_slide.Max))
                        dist_max_slide.Value = dist_max_slide.Max;
                    else
                        dist_max_slide.Value = 1e-3 * dist_max;
                    end
                end
                if (dist_min < (1e3 * dist_min_slide.Min))
                    dist_min_slide.Value = dist_min_slide.Min;
                else
                    dist_min_slide.Value = 1e-3 * dist_min;
                end
                dist_min_edit.String = sprintf('%3.1f', (1e-3 * dist_min));
            end
            dist_max_edit.String = sprintf('%3.1f', (1e-3 * dist_max));
            update_dist_range
        else
            if (dist_max > (1e3 * dist_max_slide.Max))
                dist_max_slide.Value = dist_max_slide.Max;
            else
                dist_max_slide.Value = 1e-3 * dist_max;
            end
        end
        dist_max_slide.Enable = 'off';
        drawnow
        dist_max_slide.Enable = 'on';
    end

%% Reset minimum distance

    function reset_dist_min(src, event)
        if (dist_min_ref < (1e3 * dist_min_slide.Min))
            dist_min_slide.Value = dist_min_slide.Min;
        else
            dist_min_slide.Value = 1e-3 * dist_min_ref;
        end
        dist_min_edit.String = sprintf('%3.1f', (1e-3 * dist_min_ref));
        dist_min            = dist_min_ref;
        update_dist_range
    end

%% Reset maximum distance

    function reset_dist_max(src, event)
        if (dist_max_ref > (1e3 * dist_max_slide.Max))
            dist_max_slide.Value = dist_max_slide.Max;
        else
            dist_max_slide.Value = 1e-3 * dist_max_ref;
        end
        dist_max_edit.String = sprintf('%3.1f', (1e-3 * dist_max_ref));
        dist_max            = dist_max_ref;
        update_dist_range
    end

%% Update distance range

    function update_dist_range(src, event)
        ax_radar.XLim = 1e-3 .* [dist_min dist_max];
        narrow_cb
    end

%% Reset both distance (x) and two-way traveltime (y)

    function reset_xz(src, event)
        reset_twtt_min
        reset_twtt_max
        reset_dist_min
        reset_dist_max
    end

%% Zoom on/off

	function zoom_on(src, event)
		zoom(pk_gui, 'on')
	end
	
	function zoom_off(src, event)
		zoom(pk_gui, 'off')
	end

%% Adjust slider limits after panning or zooming

    function pan_zoom(src, event)
        if ~load_done
            status_box.String = 'No data loaded yet.';
            return
        end
        tmp1                = ax_radar.XLim;
        if (tmp1(1) < (1e-3 * dist_min_ref))
            reset_dist_min
        else
            if (tmp1(1) < dist_min_slide.Min)
                tmp1(1)     = 1e-3 * dist_min_ref;
            end
            dist_min_slide.Value = tmp1(1);
            dist_min_edit.String = sprintf('%3.1f', tmp1(1));
            dist_min        = 1e3 * tmp1(1);
        end
        if (tmp1(2) > (1e-3 * dist_max_ref))
            reset_dist_max
        else
            if (tmp1(2) > dist_max_slide.Max)
                tmp1(2)     = 1e-3 * dist_max_ref;
            end
            dist_max_slide.Value = tmp1(2);
            dist_max_edit.String = sprintf('%3.1f', tmp1(2));
            dist_max        = 1e3 * tmp1(2);
        end
        tmp1                = ax_radar.YLim;
        if (tmp1(1) < (1e6 * twtt_min_ref))
            reset_twtt_min
        else
            if (((1e6 * twtt_max_ref) - (tmp1(1) - (1e6 * twtt_min_ref))) < twtt_min_slide.Min)
                twtt_min_slide.Value = twtt_min_slide.Min;
            elseif (((1e6 * twtt_max_ref) - (tmp1(1) - (1e6 * twtt_min_ref))) > twtt_min_slide.Max)
                twtt_max_slide.Value = twtt_min_slide.Max;
            else
                twtt_min_slide.Value = ((1e6 * twtt_max_ref) - (tmp1(1) - (1e6 * twtt_min_ref)));
            end
            twtt_min_edit.String = sprintf('%3.1f', tmp1(1));
            twtt_min        = 1e-6 * tmp1(1);
        end
        if (tmp1(2) > (1e6 * twtt_max_ref))
            reset_twtt_max
        else
            if (((1e6 * twtt_max_ref) - (tmp1(2) - (1e6 * twtt_min_ref))) < twtt_max_slide.Min)
                twtt_max_slide.Value = twtt_max_slide.Min;
            elseif (((1e6 * twtt_max_ref) - (tmp1(2) - (1e6 * twtt_min_ref))) > twtt_max_slide.Max)
                twtt_max_slide.Value = twtt_max_slide.Max;
            else
                twtt_max_slide.Value = ((1e6 * twtt_max_ref) - (tmp1(2) - (1e6 * twtt_min_ref)));
            end
            twtt_max_edit.String = sprintf('%3.1f', tmp1(2));
            twtt_max        = 1e-6 * tmp1(2);
        end
        narrow_cb
    end

%%  Arrow pan/zoom functions

    function pan_left(src, event)
        tmp1                = dist_max - dist_min;
        tmp2                = dist_min - (0.25 * tmp1);
        if (tmp2 < dist_min_ref)
            dist_min        = dist_min_ref;
        else
            dist_min        = tmp2;
        end
        dist_max            = dist_min + tmp1;
        if (dist_max > dist_max_ref)
            dist_max    	= dist_max_ref;
        end
        pan_horiz
    end

    function pan_right(src, event)
        tmp1                = dist_max - dist_min;
        tmp2                = dist_max + (0.25 * tmp1);
        if (tmp2 > dist_max_ref)
            dist_max        = dist_max_ref;
        else
            dist_max        = tmp2;
        end
        dist_min			= dist_max - tmp1;
        if (dist_min < dist_min_ref)
            dist_min        = dist_min_ref;
        end
        pan_horiz
    end

    function pan_horiz(src, event)
        dist_min_edit.String = sprintf('%3.1f', (1e-3 * dist_min));
        dist_max_edit.String = sprintf('%3.1f', (1e-3 * dist_max));
        if (dist_min < (1e3 * dist_min_slide.Min))
            dist_min_slide.Value = dist_min_slide.Min;
        else
            dist_min_slide.Value = 1e-3 * dist_min;
        end
        if (dist_max > (1e3 * dist_max_slide.Max))
            dist_max_slide.Value = dist_max_slide.Max;
        else
            dist_max_slide.Value = 1e-3 * dist_max;
        end
		update_dist_range
    end

    function pan_up(src, event)
        tmp1                = twtt_max - twtt_min;
        tmp2                = twtt_min - (0.25 * tmp1);
        if (tmp2 < twtt_min_ref)
            twtt_min        = twtt_min_ref;
        else
            twtt_min        = tmp2;
        end
        twtt_max            = twtt_min + tmp1;
        pan_vert
    end

    function pan_down(src, event)
        tmp1                = twtt_max - twtt_min;
        tmp2                = twtt_max + (0.25 * tmp1);
        if (tmp2 > twtt_max_ref)
            twtt_max        = twtt_max_ref;
        else
            twtt_max        = tmp2;
        end
        twtt_min            = twtt_max - tmp1;
        pan_vert
    end

    function pan_vert(src, event)
        twtt_min_edit.String = sprintf('%3.1f', (1e6 * twtt_min));
        twtt_max_edit.String = sprintf('%3.1f', (1e6 * twtt_max));
        if ((1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref))) < twtt_min_slide.Min)
            twtt_min_slide.Value = twtt_min_slide.Min;
        elseif ((1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref))) > twtt_min_slide.Max)
            twtt_max_slide.Value = twtt_min_slide.Max;
        else
            twtt_min_slide.Value = (1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref)));
        end
        if ((1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref))) < twtt_max_slide.Min)
            twtt_max_slide.Value = twtt_max_slide.Min;
        elseif ((1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref))) > twtt_max_slide.Max)
            twtt_max_slide.Value = twtt_max_slide.Max;
        else
            twtt_max_slide.Value = (1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref)));
        end
        update_twtt_range
    end

    function zoom_in(src, event)
        ax_radar.XLim = 1e-3 .* [(dist_min + (0.25 * (dist_max - dist_min))) (dist_max - (0.25 * (dist_max - dist_min)))];
		ax_radar.YLim = (1e6 .* [(twtt_min + (0.25 * (twtt_max - twtt_min))) (twtt_max - (0.25 * (twtt_max - twtt_min)))]);
        pan_zoom
    end

    function zoom_out(src, event)
        ax_radar.XLim = 1e-3 .* [(dist_min - (0.25 * (dist_max - dist_min))) (dist_max + (0.25 * (dist_max - dist_min)))];
		ax_radar.YLim = (1e6 .* [(twtt_min - (0.25 * (twtt_max - twtt_min))) (twtt_max + (0.25 * (twtt_max - twtt_min)))]);
        pan_zoom
    end

%% Switch display type

    function disp_radio(~, event)
        disp_type           = event.NewValue.String;
        switch disp_type
            case 'twtt'
                plot_twtt
            case 'depth'
                plot_depth
			case 'norm'
				plot_norm
			case 'clutter'
				plot_clutter
            case 'flat'
                plot_flat
        end
        [disp_check(:).Enable] = deal('off');
        drawnow
        [disp_check(:).Enable] = deal('on');
    end

%% Plot traveltime

    function plot_twtt(src, event)
        if ~load_done
            status_box.String = 'Data not yet loaded.';
            return
        end
		if (cmap_curr ~= 1)
            cmap_curr = 1;
            change_cmap
		end
		if isgraphics(p_data)
			[p_data.XData, p_data.CData] = deal((1e-3 .* data_cat.dist_lin), amp);
		else
			p_data			= imagesc((1e-3 .* data_cat.dist_lin), (1e6 .* data_cat.twtt), amp, [db_min db_max]);
		end
        disp_type           = 'twtt';
		update_pk_plot
        cb_def
    end

%% Plot adjusted traveltime (depth)

    function plot_depth(src, event)
        if ~load_done
            status_box.String = 'Data not yet loaded.';
            return
        end
		if (cmap_curr ~= 1)
            cmap_curr = 1;
            change_cmap
		end
        [p_data.XData, p_data.CData] = deal((1e-3 .* data_cat.dist_lin), amp_depth);
        disp_type           = 'depth';
		update_pk_plot
        cb_def
    end

%% Plot thickness-normalized radargram

	function plot_norm(src, event)
        if ~(surf_avail && bed_avail)
            disp_group.SelectedObject = disp_check(1);
            disp_type       = 'twtt';
            plot_twtt
			status_box.String = 'Cannot show thickness-normalized view without surface and bed picks.';
            return
        end
		if (cmap_curr ~= 1)
            cmap_curr = 1;
            change_cmap
		end
		if ~norm_done
			do_norm
		end
        [p_data.XData, p_data.CData] = deal((1e-3 .* data_cat.dist_lin), amp_norm);
        disp_type           = 'norm';
		update_pk_plot
		cb_def
	end

%% Plot cluttergram

    function plot_clutter(src, event)
		if (cmap_curr ~= 1)
            cmap_curr = 1;
            change_cmap
		end
        [p_data.XData, p_data.CData] = deal((1e-3 .* data_cat.dist_lin), data_cat.clutter);
        disp_type           = 'clutter';
		update_pk_plot
		cb_def
    end

%% Plot layer-flattened radargram

    function plot_flat(src, event)
        if ~flat_done
            disp_group.SelectedObject = disp_check(1);
            disp_type       = 'twtt';
            plot_twtt
            return
        end
		if (cmap_curr ~= 1)
            cmap_curr = 1;
            change_cmap
		end
        [p_data.XData, p_data.CData] = deal((1e-3 .* data_cat.dist_lin), amp_flat);
        disp_type           = 'flat';
		update_pk_plot
		cb_def
    end

%% dB defaults for colorbar

	function cb_def(src, event)
        narrow_cb
        cbl.String = '(dB)';
        cb_min_slide.Min = db_min_ref; cb_min_slide.Max = db_max_ref; cb_min_slide.Value = db_min;
        cb_max_slide.Min = db_min_ref; cb_max_slide.Max = db_max_ref; cb_max_slide.Value = db_max;
        cb_min_edit.String = sprintf('%3.0f', db_min);
        cb_max_edit.String = sprintf('%3.0f', db_max);
	end

%% Update picks display

	function update_pk_plot(src, event)
		
		switch disp_type
			
			case {'twtt' 'clutter'}
				
				for ii = 1:pk.num_layer
					tmp1	= find(~isnan(pk_ind_z(ii, :)));
					[p_pk(ii).XData, p_pk(ii).YData] = deal((1e-3 .* data_cat.dist_lin(tmp1)), (1e6 .* data_cat.twtt(pk_ind_z(ii, tmp1))'));
				end
				if surf_avail
                	[p_surf.XData, p_surf.YData] = deal((1e-3 .* data_cat.dist_lin(~isnan(ind_surf))), (1e6 .* data_cat.twtt(ind_surf(~isnan(ind_surf)))));
				end
				if bed_avail
                	[p_bed.XData, p_bed.YData] = deal((1e-3 .* data_cat.dist_lin(~isnan(ind_bed))), (1e6 .* data_cat.twtt(ind_bed(~isnan(ind_bed)))));
				end
			
			case 'depth'
				
				for ii = 1:pk.num_layer
					tmp1	= find(~isnan(pk_ind_z(ii, :)) & ~isnan(ind_surf));
					[p_pk(ii).XData, p_pk(ii).YData] = deal((1e-3 .* data_cat.dist_lin(tmp1)), (1e6 .* (data_cat.twtt(pk_ind_z(ii, tmp1))' - data_cat.twtt_surf(tmp1) + data_cat.twtt(1))));
				end
				if surf_avail
					[p_surf.XData, p_surf.YData] = deal(NaN);
				end
				if bed_avail
                    tmp1    = find(~isnan(ind_bed) & ~isnan(ind_surf));
                    [p_bed.XData, p_bed.YData] = deal((1e-3 .* data_cat.dist_lin(tmp1)), (1e6 .* (data_cat.twtt_bed(tmp1) - data_cat.twtt_surf(tmp1) + data_cat.twtt(1))));
				end
			
			case 'norm'
				
				for ii = 1:pk.num_layer
					tmp1	= find(~isnan(pk_ind_z_norm(ii, :)));
					[p_pk(ii).XData, p_pk(ii).YData] = deal((1e-3 .* data_cat.dist_lin(tmp1)), (1e6 .* data_cat.twtt(pk_ind_z_norm(ii, tmp1))'));
				end
				[p_surf.XData, p_surf.YData] = deal((1e-3 .* data_cat.dist_lin(~isnan(ind_surf_norm))), (1e6 .* data_cat.twtt(ind_surf_norm(~isnan(ind_surf_norm)))'));
				[p_bed.XData, p_bed.YData] = deal((1e-3 .* data_cat.dist_lin(~isnan(ind_bed_norm))), (1e6 .* data_cat.twtt(ind_bed_norm(~isnan(ind_bed_norm)))'));
			
			case 'flat'
				
				for ii = 1:pk.num_layer
					if isempty(find(~isnan(pk_ind_z_flat(ii, :)), 1))
						[p_pk(ii).XData, p_pk(ii).YData] = deal(NaN);
					else
						tmp1= find(~isnan(pk_ind_z_flat(ii, :)));
						[p_pk(ii).XData, p_pk(ii).YData] = deal((1e-3 .* data_cat.dist_lin(tmp1)), (1e6 .* data_cat.twtt(pk_ind_z_flat(ii, tmp1))'));
					end
				end
				if surf_avail
                    tmp1    = ind_surf_flat;
                    [p_surf.XData, p_surf.YData] = deal((1e-3 .* data_cat.dist_lin(~isnan(tmp1))), (1e6 .* data_cat.twtt(tmp1(~isnan(tmp1)))'));
				end
				if bed_avail
                    tmp1    = ind_bed_flat;
                    [p_bed.XData, p_bed.YData] = deal((1e-3 .* data_cat.dist_lin(~isnan(tmp1))), (1e6 .* data_cat.twtt(tmp1(~isnan(tmp1)))'));
				end
		end
		show_int_core
		show_surfbed
		show_pk
	end

%% Calculate thickness-normalized radargram

	function do_norm(src, event)
		if ~isempty(find(isnan(ind_surf), 1))
			ind_surf_smooth = fillmissing(ind_surf, 'spline', 'EndValues', 'extrap', 'MaxGap', 10);
		else
			ind_surf_smooth	= ind_surf;
		end
		if ~isempty(find(isnan(ind_bed), 1))
			ind_bed_smooth = fillmissing(ind_bed, 'spline', 'EndValues', 'extrap', 'MaxGap', 40);
		else
			ind_bed_smooth	= ind_bed;
		end
		
		ind_surf_smooth		= smoothdata(ind_surf_smooth, 'rlowess', (length_smooth / median(diff(1e-3 .* data_cat.dist_lin), 'omitnan'))); % smoothed surface index
		ind_bed_smooth		= smoothdata(ind_bed_smooth, 'rlowess', (length_smooth / median(diff(1e-3 .* data_cat.dist_lin), 'omitnan'))); % smoothed bed index
		ind_thick_smooth	= ind_bed_smooth - ind_surf_smooth; % thickness index
		tmp1				= (1:num_sample_trim)'; % normalization vector
		tmp2				= (tmp1(:, ones(1, data_cat.num_trace)) - ind_surf_smooth(ones(num_sample_trim, 1), :)) ./ ind_thick_smooth(ones(num_sample_trim, 1), :); % normalization matrix
		[ind_bed_norm, ind_surf_norm] ...
							= deal(NaN(1, data_cat.num_trace));
		for ii = find(~isnan(ind_surf_smooth))
			ind_surf_norm(ii) ...
							= interp1(tmp1, tmp2(:, ii), ind_surf(ii));
		end
		ind_surf_norm		= round(ind_surf_norm .* num_sample_trim);
		ind_surf_norm((ind_surf_norm < 1) | (ind_surf_norm > num_sample_trim)) ...
							= NaN;
		for ii = find(~isnan(ind_bed_smooth))
			ind_bed_norm(ii)= round(interp1(tmp1, tmp2(:, ii), ind_bed(ii)));
		end
		ind_bed_norm		= round(ind_bed_norm .* num_sample_trim);
		ind_bed_norm((ind_bed_norm < 1) | (ind_bed_norm > num_sample_trim)) ...
							= NaN;
		
		% rescaled amplitude at normalized thickness
		tmp3				= linspace(0, 1, num_sample_trim)';
		tmp4				= find(~isnan(ind_surf_smooth) & ~isnan(ind_thick_smooth) & (ind_thick_smooth > 0));
		
		amp_norm			= NaN(num_sample_trim, data_cat.num_trace, 'single');
		if parallel_check
			tmp2			= tmp2(:, tmp4);
			tmp5			= amp(:, tmp4);
			parfor ii = 1:length(tmp4)
				tmp5(:, ii)	...
							= interp1(tmp2(:, ii), tmp5(:, ii), tmp3);
			end
			amp_norm(:, tmp4) ...
							= tmp5;
			tmp5			= 0;
		else
			for ii = tmp4
				amp_norm(:, ii)	...
							= interp1(tmp2(:, ii), amp(:, ii), tmp3);
			end
		end
		norm_done			= true;
	end

%% Show intersections

	 function show_int_core(src, event)
		if int_core_check
			if any(isgraphics(p_int))
				[p_int(isgraphics(p_int)).Visible] = deal('on');
				uistack(p_int(isgraphics(p_int)), 'top')
			end
			if any(isgraphics(p_core))
				[p_core(isgraphics(p_core)).Visible] = deal('on');
				uistack(p_core(isgraphics(p_core)), 'top')
			end
		else
			if any(isgraphics(p_int))
				[p_int(isgraphics(p_int)).Visible] = deal('off');
			end
			if any(isgraphics(p_core))
				[p_core(isgraphics(p_core)).Visible] = deal('off');
			end
		end
	end

%% Show surface/bed

    function show_surfbed(src, event)
        if (surf_avail || bed_avail)
			if surfbed_check
                switch disp_type
					case {'twtt' 'norm' 'clutter' 'flat'}
                        if isgraphics(p_bed)
                            p_bed.Visible = 'on';
                            uistack(p_bed, 'top')
                        end
                        if isgraphics(p_surf)
                            p_surf.Visible = 'on';
                            uistack(p_surf, 'top')
                        end
                    case 'depth'
						if isgraphics(p_bed)
                            p_bed.Visible = 'on';
                            uistack(p_bed, 'top')
						end
						if isgraphics(p_surf)
                            p_surf.Visible = 'off';
						end
                end
            else
                if isgraphics(p_bed)
                    p_bed.Visible = 'off';
                end
                if isgraphics(p_surf)
                    p_surf.Visible = 'off';
                end
			end
        elseif surfbed_check
            surfbed_check = false;
        end
	end

%% Show picked layers

    function show_pk(src, event)
        if pk_done
			if pk_check
				[p_pk(isgraphics(p_pk)).Visible] = deal('on');
            else
                [p_pk(isgraphics(p_pk)).Visible] = deal('off');
			end
        elseif pk_check
            pk_check = false;
        end
	end

%% Switch to a chunk

    function plot_chunk(src, event)
        curr_chunk          = chunk_list.Value;
        if (curr_chunk <= num_chunk)
            ax_radar.XLim = 1e-3 .* dist_chunk([curr_chunk (curr_chunk + 1)]);
        else
            ax_radar.XLim = 1e-3 .* data_cat.dist_lin([1 end]);
        end
        tmp1                = 1e3 .* ax_radar.XLim;
        [dist_min, dist_max] = deal(tmp1(1), tmp1(2));
        if (dist_min < (1e3 * dist_min_slide.Min))
            dist_min_slide.Value = dist_min_slide.Min;
        elseif (dist_min > (1e3 * dist_min_slide.Max))
            dist_min_slide.Value = dist_min_slide.Max;
        else
            dist_min_slide.Value = 1e-3 * dist_min;
        end
        if (dist_max < (1e3 * dist_max_slide.Min))
            dist_max_slide.Value = dist_max_slide.Min;
        elseif (dist_max > (1e3 * dist_max_slide.Max))
            dist_max_slide.Value = dist_max_slide.Max;
        else
            dist_max_slide.Value = 1e-3 * dist_max;
        end
        dist_min_edit.String = sprintf('%3.1f', (1e-3 * dist_min));
        dist_max_edit.String = sprintf('%3.1f', (1e-3 * dist_max));
    end

%% Adjust length of each chunk

    function adj_length_chunk(src, event)
        length_chunk        = abs(ceil(str2double(length_chunk_edit.String)));
        % make the horizontal chunks of data to pick (too big is unwieldy)
        num_chunk           = floor((1e-3 * (data_cat.dist_lin(end) - data_cat.dist_lin(1))) ./ length_chunk) + 1;
        dist_chunk          = data_cat.dist_lin(1):(1e3 * length_chunk):data_cat.dist_lin(end);
		dist_chunk(end + 1)	= data_cat.dist_lin(end);
        chunk_list.String = [num2cell(1:num_chunk) 'full']; chunk_list.Value = (num_chunk + 1);
        status_box.String = ['Display chunk length adjusted to ' num2str(length_chunk) ' km.'];
        length_chunk_edit.String = num2str(length_chunk);
    end

%% Adjust length of maximum picking length

    function adj_length_pk_max(src, event)
        tmp1				= abs(ceil(str2double(length_pk_max_edit.String)));
		if ~isnan(tmp1)
			length_pk_max	= tmp1;
        	status_box.String = ['Maximum picking length each side adjusted to ' num2str(length_pk_max) ' km.'];
        	length_pk_max_edit.String = num2str(length_pk_max);
		else
        	length_pk_max_edit.String = num2str(length_pk_max);
		end
    end

%% Adjust number of indices above/below layer to search for peak

    function adj_num_win(src, event)
        if ~ceil(str2double(num_win_edit.String))
            num_win_edit.String = num2str(pk.num_win);
            status_box.String = 'Search window must be non-zero.';
            return
        end
        pk.num_win          = abs(ceil(str2double(num_win_edit.String)));
        num_win_edit.String = num2str(pk.num_win);
        status_box.String = ['Vertical search window adjusted to +/- ' num2str(pk.num_win) ' sample(s).'];
    end

%% Change colormap

    function change_cmap(src, event)
        pk_gui.Colormap	= cmaps{cmap_curr};
    end

%% Pop out figure

    function pop_fig(src, event)
        if ~load_done
            status_box.String = 'Data must be loaded prior to popping out a figure.';
            return
        end
        % make a simple figure that also gets saved
        set(0, 'DefaultFigureWindowStyle', 'default')
        figure('Position', [100 100 1200 800])
        axis ij tight
        colormap(cmaps{cmap_curr})
        switch disp_type
            case 'twtt'
				imagesc((1e-3 .* data_cat.dist_lin), (1e6 .* data_cat.twtt), amp, [db_min db_max])				
				axis([(1e-3 .* [dist_min dist_max]) (1e6 .* [twtt_min twtt_max])])
                clim([db_min db_max])
                if surfbed_check
                    if surf_avail
                        line((1e-3 .* data_cat.dist_lin), (1e6 .* data_cat.twtt_surf), 'LineStyle', 'none', 'Marker', '.', 'Color', 'm', 'MarkerSize', 8)
                    end
                    if bed_avail
                        line((1e-3 .* data_cat.dist_lin), (1e6 .* data_cat.twtt_bed), 'LineStyle', 'none', 'Marker', '.', 'Color', 'm', 'MarkerSize', 8)
                    end
                end
				if pk_check
                    for ii = 1:pk.num_layer
                        if ~isempty(find(~isnan(pk_ind_z(ii, :)), 1))
                            tmp1 = line((1e-3 .* data_cat.dist_lin(~isnan(pk_ind_z(ii, :)))), (1e6 .* data_cat.twtt(pk_ind_z(ii, ~isnan(pk_ind_z(ii, :))))'), ...
								'LineStyle', 'none', 'Marker', '.', 'Color', pk_color(ii, :), 'MarkerSize', 8);
                            if (ii == curr_layer)
                                tmp1.MarkerSize = 24;
                            end
                        end
                    end
				end
                text((1e-3 * (dist_max + (0.015 * (dist_max - dist_min)))), (1e6 * (twtt_min - (0.04 * (twtt_max - twtt_min)))), '(dB)', 'Color', 'k', 'FontSize', 20)
            case 'depth'
                imagesc((1e-3 .* data_cat.dist_lin), (1e6 .* data_cat.twtt), amp_depth, [db_min db_max])
				axis([(1e-3 .* [dist_min dist_max]) (1e6 .* [twtt_min twtt_max])])
                clim([db_min db_max])
                if surfbed_check
                    if bed_avail
                        tmp1= find(~isnan(ind_bed) & ~isnan(ind_surf));
                        line((1e-3 .* data_cat.dist_lin(tmp1)), (1e6 .* (data_cat.twtt_bed(tmp1) - data_cat.twtt_surf(tmp1))), 'LineStyle', 'none', 'Marker', '.', 'Color', 'm', 'MarkerSize', 8)
                    end
                end
				if pk_check
                    for ii = 1:pk.num_layer
                        if ~isempty(find(~isnan(pk_ind_z(ii, :)), 1))
                            tmp1 = find(~isnan(pk_ind_z(ii, :)) & ~isnan(ind_surf));
                            tmp2 = line((1e-3 .* data_cat.dist_lin(tmp1)), (1e6 .* data_cat.twtt(pk_ind_z(ii, tmp1) - ind_surf(tmp1) + 1)'), 'LineStyle', 'none', 'Marker', '.', 'Color', pk_color(ii, :), 'MarkerSize', 8);
                            if (ii == curr_layer)
                                tmp2.MarkerSize = 24;
                            end
                        end
                    end
				end
                text((1e-3 * (dist_max + (0.015 * (dist_max - dist_min)))), (1e6 * (twtt_min - (0.04 * (twtt_max - twtt_min)))), '(dB)', 'Color', 'k', 'FontSize', 20)
			case 'norm'
                imagesc((1e-3 .* data_cat.dist_lin), (1e6 .* data_cat.twtt), amp_norm, [db_min db_max])
                if surfbed_check
                    if surf_avail
                        line((1e-3 .* data_cat.dist_lin(~isnan(ind_surf_norm))), (1e6 .* data_cat.twtt(ind_surf_norm(~isnan(ind_surf_norm)))'), 'LineStyle', 'none', 'Marker', '.', 'Color', 'g', 'MarkerSize', 8)
                    end
                    if bed_avail
                        line((1e-3 .* data_cat.dist_lin(~isnan(ind_bed_norm))), (1e6 .* data_cat.twtt(ind_bed_norm(~isnan(ind_bed_norm)))'), 'LineStyle', 'none', 'Marker', '.', 'Color', 'g', 'MarkerSize', 8)
                    end
                end
				if pk_check
                    for ii = 1:pk.num_layer
                        tmp1 = line((1e-3 .* data_cat.dist_lin(~isnan(pk_ind_z_norm(ii, :)))), (1e6 .* data_cat.twtt(pk_ind_z_norm(ii, ~isnan(pk_ind_z_norm(ii, :))))'), ...
									'LineStyle', 'none', 'Marker', '.', 'Color', pk_color(ii, :), 'MarkerSize', 8);
                        if (ii == curr_layer)
                            tmp1.MarkerSize = 24;
                        end
                    end
				end
                text((1e-3 * (dist_max + (0.015 * (dist_max - dist_min)))), (1e6 * (twtt_min - (0.04 * (twtt_max - twtt_min)))), '(dB)', 'Color', 'k', 'FontSize', 20)
			case 'clutter'
				imagesc((1e-3 .* data_cat.dist_lin), (1e6 .* data_cat.clutter), amp, [db_min db_max])				
				axis([(1e-3 .* [dist_min dist_max]) (1e6 .* [twtt_min twtt_max])])
                clim([db_min db_max])
                if surfbed_check
                    if surf_avail
                        line((1e-3 .* data_cat.dist_lin), (1e6 .* data_cat.twtt_surf), 'LineStyle', 'none', 'Marker', '.', 'Color', 'm', 'MarkerSize', 8)
                    end
                    if bed_avail
                        line((1e-3 .* data_cat.dist_lin), (1e6 .* data_cat.twtt_bed), 'LineStyle', 'none', 'Marker', '.', 'Color', 'm', 'MarkerSize', 8)
                    end
                end
				if pk_check
                    for ii = 1:pk.num_layer
                        if ~isempty(find(~isnan(pk_ind_z(ii, :)), 1))
                            tmp1 = line((1e-3 .* data_cat.dist_lin(~isnan(pk_ind_z(ii, :)))), (1e6 .* data_cat.twtt(pk_ind_z(ii, ~isnan(pk_ind_z(ii, :))))'), ...
								'LineStyle', 'none', 'Marker', '.', 'Color', pk_color(ii, :), 'MarkerSize', 8);
                            if (ii == curr_layer)
                                tmp1.MarkerSize = 24;
                            end
                        end
                    end
				end
            case 'flat'
                imagesc((1e-3 .* data_cat.dist_lin), (1e6 .* data_cat.twtt), amp_flat, [db_min db_max])
                if surfbed_check
                    if surf_avail
                        tmp1= ind_surf_flat;
                        line((1e-3 .* data_cat.dist_lin(~isnan(tmp1))), (1e6 .* data_cat.twtt(tmp1(~isnan(tmp1)))'), 'LineStyle', 'none', 'Marker', '.', 'Color', 'g', 'MarkerSize', 8)
                    end
                    if bed_avail
                        tmp1= ind_bed_flat;
                        line((1e-3 .* data_cat.dist_lin(~isnan(tmp1))), (1e6 .* data_cat.twtt(tmp1(~isnan(tmp1)))'), 'LineStyle', 'none', 'Marker', '.', 'Color', 'g', 'MarkerSize', 8)
                    end
                end
				if pk_check
                    for ii = 1:pk.num_layer
                        tmp1 = line((1e-3 .* data_cat.dist_lin(~isnan(pk_ind_z_flat(ii, :)))), (1e6 .* data_cat.twtt(pk_ind_z_flat(ii, ~isnan(pk_ind_z_flat(ii, :))))'), ...
									'LineStyle', 'none', 'Marker', '.', 'Color', pk_color(ii, :), 'MarkerSize', 8);
                        if (ii == curr_layer)
                            tmp1.MarkerSize = 24;
                        end
                    end
				end
                text((1e3 * (dist_max + (0.015 * (dist_max - dist_min)))), (1e6 * (twtt_min - (0.04 * (twtt_max - twtt_min)))), '(dB)', 'Color', 'k', 'FontSize', 20)
        end
        set(gca, 'FontSize', 20, 'Layer', 'top')
        xlabel('Distance (km)')
        ylabel('Traveltime ({\mu}s)')
        title(file_data(1:(end - 4)), 'FontWeight', 'bold', 'Interpreter', 'none')
        colorbar('FontSize', 20)
        if grid_check
            grid on
        end
        box on
        set(0, 'DefaultFigureWindowStyle', 'docked')
    end

%% Toggle gridlines

	function toggle_grid(src, event)
        if grid_check
            [ax_radar.XGrid, ax_radar.YGrid] = deal('on');
        else
            [ax_radar.XGrid, ax_radar.YGrid] = deal('off');
        end
    end

%% Narrow color axis to +/- 2 standard deviations of current mean value

    function narrow_cb(src, event)
        if (cbfix_check2 && ~isequal(narrow_ax, [ax_radar.XLim ax_radar.YLim]))
            tmp1            = zeros(2);
            tmp1(1, :)      = interp1(data_cat.dist_lin, ind_num_trace, [dist_min dist_max], 'nearest', 'extrap');
            tmp1(2, :)      = interp1(data_cat.twtt, 1:num_sample_trim, [twtt_min twtt_max], 'nearest', 'extrap');
            switch disp_type
                case 'twtt'
                    tmp1    = amp(tmp1(2, 1):10:tmp1(2, 2), tmp1(1, 1):10:tmp1(1, 2));
                case 'depth'
                    tmp1    = amp_depth(tmp1(2, 1):10:tmp1(2, 2), tmp1(1, 1):10:tmp1(1, 2));
				case 'norm'
                    tmp1    = amp_norm(tmp1(2, 1):10:tmp1(2, 2), tmp1(1, 1):10:tmp1(1, 2));
				case 'clutter'
                    tmp1    = data_cat.clutter(tmp1(2, 1):10:tmp1(2, 2), tmp1(1, 1):10:tmp1(1, 2));			
                case 'flat'
                    tmp1    = amp_flat(tmp1(2, 1):10:tmp1(2, 2), tmp1(1, 1):10:tmp1(1, 2));
            end
            tmp2            = NaN(1, 2);
            [tmp2(1), tmp2(2)] ...
                            = deal(mean(tmp1(~isinf(tmp1)), 'omitnan'), std(tmp1(~isinf(tmp1)), 'omitnan'));
            tmp1			= zeros(1, 2);
			tmp4			= [db_min_ref db_max_ref];
            if ((tmp2(1) - (2 * tmp2(2))) < tmp4(1))
                tmp1(1)     = tmp4(1);
            else
                tmp1(1)     = tmp2(1) - (2 * tmp2(2));
            end
            if ((tmp2(1) + (2 * tmp2(2))) > tmp4(2))
                tmp1(2)     = tmp4(2);
            else
                tmp1(2)     = tmp2(1) + (2 * tmp2(2));
            end
            [tmp4(1), tmp4(2)] ...
                            = deal(tmp1(1), tmp1(2));
            if (tmp4(1) < cb_min_slide.Min)
                tmp4(1)     = db_min_ref;
            end
            cb_min_slide.Value = tmp4(1);
            if (tmp4(2) > cb_max_slide.Max)
                tmp4(2)     = db_max_ref;
            end
            if any(isnan(tmp4))
                status_box.String = 'Color scale could not be changed due to NaN color values.';
                return
            end
            if ~issorted(tmp4)
                tmp4        = sort(tmp4);
            end
            cb_max_slide.Value = tmp4(2);
            cb_min_edit.String = sprintf('%3.0f', tmp4(1));
            cb_max_edit.String = sprintf('%3.0f', tmp4(2));
            [db_min, db_max]= deal(tmp4(1), tmp4(2));
			ax_radar.CLim = tmp4;
			narrow_ax		= [ax_radar.XLim ax_radar.YLim];
        end
    end

%% Update layer colors

	function update_pk_color(src, event)
		pk_color			= repmat(pk_color_def, ceil(pk.num_layer / size(pk_color_def, 1)), 1); % extend predefined color pattern
		pk_color			= pk_color(1:pk.num_layer, :); % trim to current range
	end

%% Keyboard shortcuts for various functions

    function keypress(src, event)
		switch event.Key
			case '1'
				chunk_list.Value = 1;
				plot_chunk
            case '3'
                if pk_check
                    pk_check = false;
                else
                    pk_check = true;
                end
                show_pk
            case '4'
                if surfbed_check
                    surfbed_check = false;
                else
                    surfbed_check = true;
                end
                show_surfbed
            case '5'
                if int_core_check
                    int_core_check = false;
                else
                    int_core_check = true;
                end
                show_int_core
			case '7'
				pk_trimx
            case 'a'
                pk_adj
            case 'b'
				pk_surfbed
            case 'c'
				pk_focus
			case 'd'
                pk_del
            case 'e'
                reset_xz
            case 'f'
                flatten
            case 'g'
                if grid_check 
                    grid_check = false;
                else
                    grid_check = true;
                end
                toggle_grid
			case 'h'
				pk_shift
            case 'i'
                pk_load
			case 'j'
				pk_last
            case 'l'
                load_data
            case 'm'
                pk_merge
            case 'n'
                pk_next
			case 'o'
				if (pk.num_layer > 1)
					switch disp_type
						case {'twtt' 'depth' 'clutter'}
							twtt_min = data_cat.twtt(min(pk_ind_z, [], 'all', 'omitnan'));
							twtt_max = data_cat.twtt(max(pk_ind_z, [], 'all', 'omitnan'));
						case 'norm'
							twtt_min = data_cat.twtt(min(pk_ind_z_norm, [], 'all', 'omitnan'));
							twtt_max = data_cat.twtt(max(pk_ind_z_norm, [], 'all', 'omitnan'));
						case 'flat'
							twtt_min = data_cat.twtt(min(pk_ind_z_flat, [], 'all', 'omitnan'));
							twtt_max = data_cat.twtt(max(pk_ind_z_flat, [], 'all', 'omitnan'));
					end
					twtt_min_slide.Value = (1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref)));
					twtt_max_slide.Value = (1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref)));
					twtt_min_edit.String = sprintf('%3.1f', (1e6 * twtt_min));
					twtt_max_edit.String = sprintf('%3.1f', (1e6 * twtt_max));
					ax_radar.YLim = 1e6 .* [twtt_min twtt_max];
					switch disp_type
						case {'twtt' 'depth' 'clutter'}
							tmp1 = sum(~isnan(pk_ind_z));
						case 'norm'
							tmp1 = sum(~isnan(pk_ind_z_norm));
						case 'flat'
							tmp1 = sum(~isnan(pk_ind_z_flat));
					end
					tmp1	= [find(tmp1, 1) find(tmp1, 1, 'last')];
					dist_min = data_cat.dist_lin(tmp1(1));
					dist_min_slide.Value = 1e-3 * dist_min;
					dist_min_edit.String = sprintf('%3.1f', (1e-3 * dist_min));
					dist_max = data_cat.dist_lin(tmp1(2));
					dist_max_slide.Value = 1e-3 * dist_max;
					dist_max_edit.String = sprintf('%3.1f', (1e-3 * dist_max));
					if strcmp(disp_type, 'depth')
						tmp1= data_cat.twtt(min(ind_surf(tmp1(1):tmp1(2))));
						[twtt_min, twtt_max] ...
							= deal((twtt_min - tmp1), (twtt_max - tmp1));
						twtt_min_slide.Value = (1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref)));
						twtt_max_slide.Value = (1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref)));
						twtt_min_edit.String = sprintf('%3.1f', (1e6 * twtt_min));
						twtt_max_edit.String = sprintf('%3.1f', (1e6 * twtt_max));
						ax_radar.YLim = 1e6 .* [twtt_min twtt_max];
					end
					update_dist_range
				end
			case 'p'
                pk_auto
            case 'q'
                pop_fig
            case 's'
                pk_save
			case 't'
				misctest
			case 'v'
				if cross_check
					cross_check = false;
					status_box.String = 'Crossing layers checking turned off.';
				else
					cross_check = true;
					status_box.String = 'Crossing layers checking turned on.';
					pk_cross
				end
            case 'w'
                if (cmap_curr == 1)
                    cmap_curr = 2;
                else
                    cmap_curr = 1;
                end
                change_cmap
            case 'x'
                pk_split
            case 'y'
				if twttfix_check
                    twttfix_check = false;
					status_box.String = 'Vertical scale free.';
                else
                    twttfix_check = true;
					status_box.String = 'Vertical scale fixed.';
				end
            case 'z'
				zoom_on
            case 'slash'
				switch ord_poly
                    case 2
                        ord_poly ...
                            = 3;
                        status_box.String = 'Flattening polynomial switched to 3rd order';
                    case 3
                        ord_poly ...
                            = 2;
                        status_box.String = 'Flattening polynomial switched to 2nd order.';
				end
			case 'equal'
				tmp1		= interp1(data_cat.dist_lin, ind_num_trace, [dist_min dist_max], 'nearest');
				if (~isempty(find(~isnan(ind_surf(tmp1(1):tmp1(2))), 1)) && ~isempty(find(~isnan(ind_bed(tmp1(1):tmp1(2))), 1)) && any(strcmp(disp_type, {'twtt' 'depth' 'norm' 'clutter' 'flat'})))
					switch disp_type
						case {'twtt' 'depth' 'clutter'}
							tmp2	= [mean(ind_surf(tmp1(1):tmp1(2)), 'omitnan') mean(ind_bed(tmp1(1):tmp1(2)), 'omitnan')];
						case 'norm'
							tmp2	= [mean(ind_surf_norm(tmp1(1):tmp1(2)), 'omitnan') mean(ind_bed_norm(tmp1(1):tmp1(2)), 'omitnan')];
							if isnan(tmp2(1))
								tmp2(1) = 1;
							end
							if isnan(tmp2(2))
								tmp2(2) = num_sample_trim;
							end
						case 'flat'
							tmp2	= [mean(ind_surf_flat(tmp1(1):tmp1(2)), 'omitnan') mean(ind_bed_flat(tmp1(1):tmp1(2)), 'omitnan')];
							if isnan(tmp2(1))
								tmp2(1) = 1;
							end
							if isnan(tmp2(2))
								tmp2(2) = num_sample_trim;
							end
					end
					tmp2	= round([(tmp2(1) + (0.05 * diff(tmp2))) (tmp2(2) - (0.05 * diff(tmp2)))]);
					twtt_min = data_cat.twtt(tmp2(1));
					twtt_max = data_cat.twtt(tmp2(2));
					twtt_min_slide.Value = (1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref)));
					twtt_max_slide.Value = (1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref)));
					twtt_min_edit.String = sprintf('%3.1f', (1e6 * twtt_min));
					twtt_max_edit.String = sprintf('%3.1f', (1e6 * twtt_max));
					update_twtt_range
				end
			case 'quote'
				if distfix_check
                    distfix_check = false;
					status_box.String = 'Distance scale free.';
                else
                    distfix_check = true;
					status_box.String = 'Distance scale fixed.';
				end
			case 'leftbracket'
				if cbfix_check1
                    cbfix_check1 = false;
					status_box.String = 'Color scale free.';
                else
                    cbfix_check1 = true;
					status_box.String = 'Color scale fixed.';
				end
			case 'rightbracket'
				if cbfix_check2
                    cbfix_check2 = false;
					status_box.String = 'Color scale free.';
                else
                    cbfix_check2 = true;
					status_box.String = 'Color scale fixed to mean +/- 1 standard deviation of current view.';
					narrow_cb
				end
            case 'downarrow'
                pan_down
            case 'leftarrow'
                pan_left
            case 'rightarrow'
                pan_right
            case 'uparrow'
                pan_up
            case 'space'
                if ~load_done
                    return
                end
                switch disp_type
                    case 'twtt'
                        if flat_done
                            disp_group.SelectedObject = disp_check(5);
                            disp_type = 'flat';
                            plot_flat
						elseif norm_done
                            disp_group.SelectedObject = disp_check(3);
                            disp_type = 'norm';
                            plot_norm
						elseif clutter_avail
                            disp_group.SelectedObject = disp_check(4);
                            disp_type = 'clutter';
                            plot_clutter
						else
                            disp_group.SelectedObject = disp_check(2);
                            disp_type = 'depth';
                            plot_depth
                        end
                    case 'depth'
						if flat_done
                            disp_group.SelectedObject = disp_check(5);
                            disp_type = 'flat';
                            plot_flat
                        else
                            disp_group.SelectedObject = disp_check(1);
                            disp_type = 'twtt';
                            plot_twtt
						end
                    otherwise
                        disp_group.SelectedObject = disp_check(1);
                        disp_type = 'twtt';
                        plot_twtt
                end
		end
    end

%% Mouse click

    function mouse_click(src, event)
		if ((logical(pk.num_layer) || surf_avail || bed_avail))
			tmp1            = src.CurrentPoint;
			tmp2            = pk_gui.Position;
			tmp3            = ax_radar.Position;
			tmp4            = [(tmp2(1) + (tmp2(3) * tmp3(1))) (tmp2(1) + (tmp2(3) * (tmp3(1) + tmp3(3)))); (tmp2(2) + (tmp2(4) * tmp3(2))) (tmp2(2) + (tmp2(4) * (tmp3(2) + tmp3(4))))];
			if ((tmp1(1) > (tmp4(1, 1))) && (tmp1(1) < (tmp4(1, 2))) && (tmp1(2) > (tmp4(2, 1))) && (tmp1(2) < (tmp4(2, 2))))
				tmp1        = [((tmp1(1) - tmp4(1, 1)) / diff(tmp4(1, :))) ((tmp4(2, 2) - tmp1(2)) / diff(tmp4(2, :)))];
				tmp2        = [ax_radar.XLim; ax_radar.YLim];
				[ind_x_pk, ind_z_pk] ...
							= deal(((tmp1(1) * diff(tmp2(1, :))) + tmp2(1, 1)), ((tmp1(2) * diff(tmp2(2, :))) + tmp2(2, 1)));
				switch src.SelectionType
					case 'normal'
						pk_select_gui
					case 'extend' % shift / left cut
						button ...
							= 'L1';
						pk_adj
						button ...
							= NaN;
					case 'alt' % ctrl / right cut
						button ...
							= 'R1';
						pk_adj
						button ...
							= NaN;
				end
			end
		end
    end

%% Test something

    function misctest(src, event)
		
		pk_gui.KeyPressFcn = @keypress; pk_gui.WindowButtonDownFcn = @mouse_click;
		status_box.String = 'Test done.';
	end
%%
end