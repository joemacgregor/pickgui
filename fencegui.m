function fencegui
% FENCEGUI Interactive comparison of layer picks between intersecting radar segments.
%   
%   FENCEGUI loads a pair of GUIs (3D and 2D) for interactive comparison
%   of a primary segment's radar layers with those from segments that
%   intersect it. These layers must have been traced previously across
%   concatenated datasets using PICKGUI.
% 
% Joe MacGregor (NASA)
% Last updated: 12 November 2024

%% Intialize variables

% initial paths, update as needed
[path_int, path_match]		= deal('/Users/jamacgre/OneDrive - NASA/research/matlab/pickgui_v2/mat/');
path_pk						= '/Users/jamacgre/OneDrive - NASA/data/gris_strat2/arctic_cat_pk_rev_r1/';
path_data					= '/Users/jamacgre/OneDrive - NASA/data/gris_strat2/arctic_cat_r1/';
[file_int, file_match]		= deal('merge_xy_int.mat', 'layer_match_list.mat');

load([path_int file_int], 'int_all', 'name_pk');
load([path_match file_match], 'layer_match_list');

% depth defaults
[depth_min_ref, depth_max_ref] ...
                            = deal(0, 3500);
[depth_min, depth_max]      = deal(depth_min_ref, depth_max_ref);

% distance defaults
dist_int_half				= 10e3;
[dist_min_ref, dist_max_ref]= deal(zeros(1, 2), (1e3 .* ones(1, 2)));
[dist_min, dist_max]        = deal(dist_min_ref, dist_max_ref);

% dB default
[db_min_ref, db_max_ref]    = deal(repmat(-130, 1, 2), zeros(1, 2));
[db_min, db_max]            = deal(repmat(-80, 1, 2), repmat(-20, 1, 2));

% some default values
speed_vacuum                = 299792458; % speed of light in the vacuum, m/s
permitt_ice                 = 3.15; % real part of relative permittivity of pure ice at RF
speed_ice                   = speed_vacuum / sqrt(permitt_ice); % RF ice speed, m/s
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
                               0.75 0       0]; % RGB layer colors

% allocate a bunch of variables
[data_check, data_done, distfix_check, grid_check, pk_check, pk_done, zfix_check] ...
							= deal(false(1, 2));
[amp_depth, colors, depth_bed, depth_layer, depth_ref, dist_lin, file_data, file_pk, file_pk_short, ind_layer, layer_str, narrow_ax, p_core, p_core_name, p_int1, p_int2, p_pk, twtt] ...
                            = deal(cell(1, 2));
[curr_layer, curr_pk, num_layer, num_sample, num_trace] ...
                            = deal(zeros(1, 2));
[curr_ind_int, ii, ind_x_pk, ind_y_pk, jj, num_int, tmp1, tmp2, tmp3, tmp4, tmp5] ...
                            = deal(0);
[curr_int, curr_rad]		= deal(1);
curr_rad_alt                = 2;
[ax, cb_min_slide, cb_min_edit, cb_max_slide, cb_max_edit, dist_min_slide, dist_min_edit, dist_max_slide, dist_max_edit, int_check, layer_list, p_bed, p_data, rad_check] ...
							= deal(gobjects(1, 2));
match_done					= false;

size_font					= 18;
width_slide					= 0.02;

tmp_list					= [];

%% draw GUI

set(0,'DefaultFigureWindowStyle','docked')
fgui						= figure('Toolbar', 'figure', 'Name', 'FENCEGUI', 'MenuBar', 'none', 'KeyPressFcn', @keypress, 'WindowButtondownFcn', @mouse_click);
ax(1)                       = subplot('Position', [0.065 0.06 0.41 0.85]);
ax(2)                       = subplot('Position', [0.55 0.06 0.41 0.85]);

[ax(:).FontSize] = deal(size_font); [ax(:).Layer] = deal('top'); [ax(:).YDir] = deal('reverse'); [ax(:).GridLineWidth] = deal(2);

axes(ax(1))
hold on
colormap(bone)
clim([db_min(1) db_max(1)])
axis([(1e-3 .* [dist_min_ref(1) dist_max_ref(1)]) depth_min_ref depth_max_ref])
box on
grid off
ylabel('(m)');
colorbar('FontSize', size_font)
% pan/zoom callbacks
h_pan(1)                    = pan;
h_pan(1).ActionPostCallback = @panzoom;
h_zoom(1)                   = zoom;
h_zoom(1).ActionPostCallback= @panzoom;

axes(ax(2))
hold on
colormap(bone)
clim([db_min(2) db_max(2)])
axis([(1e-3 .* [dist_min_ref(2) dist_max_ref(2)]) depth_min_ref depth_max_ref])
box on
grid off
colorbar('FontSize', size_font)
% pan/zoom callbacks
h_pan(2)                    = pan;
h_pan(2).ActionPostCallback = @panzoom;
h_zoom(2)                   = zoom;
h_zoom(2).ActionPostCallback= @panzoom;

% sliders
cb_min_slide(1)             = uicontrol(fgui, 'Style', 'slider', 'Units', 'normalized', 'Position', [0.475 0.07 width_slide 0.32], 'Callback', @slide_db_min1, 'Min', db_min_ref(1), 'Max', db_max_ref(1), 'Value', db_min(1), 'SliderStep', [0.01 0.1]);
cb_max_slide(1)             = uicontrol(fgui, 'Style', 'slider', 'Units', 'normalized', 'Position', [0.475 0.50 width_slide 0.32], 'Callback', @slide_db_max1, 'Min', db_min_ref(1), 'Max', db_max_ref(1), 'Value', db_max(1), 'SliderStep', [0.01 0.1]);
cb_min_slide(2)             = uicontrol(fgui, 'Style', 'slider', 'Units', 'normalized', 'Position', [0.97 0.07 width_slide 0.32], 'Callback', @slide_db_min2, 'Min', db_min_ref(2), 'Max', db_max_ref(2), 'Value', db_min(2), 'SliderStep', [0.01 0.1]);
cb_max_slide(2)             = uicontrol(fgui, 'Style', 'slider', 'Units', 'normalized', 'Position', [0.97 0.50 width_slide 0.32], 'Callback', @slide_db_max2, 'Min', db_min_ref(2), 'Max', db_max_ref(2), 'Value', db_max(2), 'SliderStep', [0.01 0.1]);
dist_min_slide(1)           = uicontrol(fgui, 'Style', 'slider', 'Units', 'normalized', 'Position', [0.10 0.005 0.12 0.02], 'Callback', @slide_dist_min1, 'Min', (1e-3 * dist_min_ref(1)), 'Max', (1e-3 * dist_max_ref(1)), 'Value', (1e-3 * dist_min_ref(1)), 'SliderStep', [0.01 0.1]);
dist_max_slide(1)           = uicontrol(fgui, 'Style', 'slider', 'Units', 'normalized', 'Position', [0.33 0.005 0.12 0.02], 'Callback', @slide_dist_max1, 'Min', (1e-3 * dist_min_ref(1)), 'Max', (1e-3 * dist_max_ref(1)), 'Value', (1e-3 * dist_max_ref(1)), 'SliderStep', [0.01 0.1]);
dist_min_slide(2)           = uicontrol(fgui, 'Style', 'slider', 'Units', 'normalized', 'Position', [0.58 0.005 0.12 0.02], 'Callback', @slide_dist_min2, 'Min', (1e-3 * dist_min_ref(2)), 'Max', (1e-3 * dist_max_ref(2)), 'Value', (1e-3 * dist_min_ref(2)), 'SliderStep', [0.01 0.1]);
dist_max_slide(2)           = uicontrol(fgui, 'Style', 'slider', 'Units', 'normalized', 'Position', [0.81 0.005 0.12 0.02], 'Callback', @slide_dist_max2, 'Min', (1e-3 * dist_min_ref(2)), 'Max', (1e-3 * dist_max_ref(2)), 'Value', (1e-3 * dist_max_ref(2)), 'SliderStep', [0.01 0.1]);
z_min_slide					= uicontrol(fgui, 'Style', 'slider', 'Units', 'normalized', 'Position', [0.005 0.07 width_slide 0.32], 'Callback', @slide_z_min, 'Min', depth_min_ref, 'Max', depth_max_ref, 'Value', depth_min_ref, 'SliderStep', [0.01 0.1]);
z_max_slide					= uicontrol(fgui, 'Style', 'slider', 'Units', 'normalized', 'Position', [0.005 0.50 width_slide 0.32], 'Callback', @slide_z_max, 'Min', depth_min_ref, 'Max', depth_max_ref, 'Value', depth_max_ref, 'SliderStep', [0.01 0.1]);

% slider values
cb_min_edit(1)              = annotation('textbox', [0.48 0.39 0.04 0.03], 'String', num2str(db_min(1)), 'FontSize', size_font, 'Color', 'k', 'EdgeColor', 'none');
cb_min_edit(2)              = annotation('textbox', [0.965 0.39 0.04 0.03], 'String', num2str(db_min(2)), 'FontSize', size_font, 'Color', 'k', 'EdgeColor', 'none');
cb_max_edit(1)              = annotation('textbox', [0.48 0.82 0.04 0.03], 'String', num2str(db_max(1)), 'FontSize', size_font, 'Color', 'k', 'EdgeColor', 'none');
cb_max_edit(2)              = annotation('textbox', [0.965 0.82 0.04 0.03], 'String', num2str(db_max(2)), 'FontSize', size_font, 'Color', 'k', 'EdgeColor', 'none');
dist_min_edit(1)            = annotation('textbox', [0.07 0.005 0.04 0.03], 'String', num2str(1e-3 * dist_min_ref(1)), 'FontSize', size_font, 'Color', 'k', 'EdgeColor', 'none');
dist_min_edit(2)            = annotation('textbox', [0.55 0.005 0.04 0.03], 'String', num2str(1e-3 * dist_min_ref(2)), 'FontSize', size_font, 'Color', 'k', 'EdgeColor', 'none');
dist_max_edit(1)            = annotation('textbox', [0.295 0.005 0.04 0.03], 'String', num2str(1e-3 * dist_max_ref(1)), 'FontSize', size_font, 'Color', 'k', 'EdgeColor', 'none');
dist_max_edit(2)            = annotation('textbox', [0.775 0.005 0.04 0.03], 'String', num2str(1e-3 * dist_max_ref(2)), 'FontSize', size_font, 'Color', 'k', 'EdgeColor', 'none');
z_min_edit					= annotation('textbox', [0.005 0.39 0.04 0.03], 'String', num2str(depth_min_ref), 'FontSize', size_font, 'Color', 'k', 'EdgeColor', 'none');
z_max_edit					= annotation('textbox', [0.005 0.82 0.04 0.03], 'String', num2str(depth_min_ref), 'FontSize', size_font, 'Color', 'k', 'EdgeColor', 'none');
int_count_edit				= annotation('textbox', [0.63 0.95 0.04 0.03], 'String', '0 /', 'FontSize', size_font, 'FontWeight', 'bold', 'Color', 'b', 'EdgeColor', 'none');
int_tot1_edit				= annotation('textbox', [0.66 0.95 0.04 0.03], 'String', '0', 'FontSize', size_font, 'FontWeight', 'bold', 'Color', 'b', 'EdgeColor', 'none');
int_tot2_edit				= annotation('textbox', [0.21 0.92 0.04 0.03], 'String', '0', 'FontSize', size_font, 'FontWeight', 'bold', 'Color', 'b', 'EdgeColor', 'none');

% push buttons
button_param2				= {'Load picks & data'		[0.005 0.925 0.08 0.03]		@load_pk1			'b';
							   'Match'					[0.72 0.925 0.04 0.03]		@pk_match			'm';
							   'Unmatch'				[0.765 0.925 0.05 0.03]		@pk_unmatch			'm';
							   'Test'					[0.96 0.925 0.03 0.03]		@misctest			'r';
							   'Reset x/z'				[0.41 0.925 0.05 0.03]		@reset_xz1			'r';
							   'Reset x/z'				[0.91 0.925 0.04 0.03]		@reset_xz2			'r';
							   'Reset'					[0.49 0.03 0.03 0.03]		@reset_db_min1		'r';
							   'Reset'					[0.49 0.46 0.03 0.03]		@reset_db_max1		'r';
							   'Reset'					[0.965 0.03 0.03 0.03]		@reset_db_min2		'r';
							   'Reset'					[0.965 0.46 0.03 0.03]		@reset_db_max2		'r'};
for ii  = 1:size(button_param2, 1)
	uicontrol(fgui, 'Style', 'pushbutton', 'String', button_param2{ii, 1}, 'Units', 'normalized', 'Position', button_param2{ii, 2}, 'Callback', button_param2{ii, 3}, 'FontSize', size_font, 'ForegroundColor', button_param2{ii, 4})
end

% fixed text annotations        position                string              color
text_param2                 = {[0.49 0.42 0.03 0.03]    'dB_{min}'          'k';
                               [0.49 0.85 0.03 0.03]    'dB_{max}'          'k';
                               [0.965 0.42 0.03 0.03]   'dB_{min}'          'k';
                               [0.965 0.85 0.03 0.03]   'dB_{max}'          'k';
                               [0.18 0.965 0.03 0.03]   'Layer'             'm';
                               [0.52 0.965 0.03 0.03]   'Layer'             'm';
                               [0.035 0.005 0.03 0.03]  'dist_{min}'        'k';
                               [0.515 0.005 0.03 0.03]  'dist_{min}'        'k';
                               [0.26 0.005 0.03 0.03]   'dist_{max}'        'k';
                               [0.74 0.005 0.03 0.03]   'dist_{max}'        'k';
                               [0.005 0.42 0.03 0.03]   'depth_{max}'       'k';
                               [0.005 0.85 0.03 0.03]   'depth_{min}'       'k';
                               [0.10 0.92 0.08 0.03]    'Intersection #'    'b';
                               [0.54 0.92 0.08 0.03]    'Intersections'     'b';
                               [0.815 0.925 0.08 0.03]  'Nearest'           'm';
                               [0.865 0.925 0.08 0.03]  'Match'             'm'};
for ii = 1:size(text_param2, 1)
	if ~ispc
		annotation('textbox', text_param2{ii, 1}, 'String', text_param2{ii, 2}, 'FontSize', size_font, 'Color', text_param2{ii, 3}, 'EdgeColor', 'none', 'FontWeight', 'bold')
	else
		annotation('textbox', text_param2{ii, 1}, 'String', text_param2{ii, 2}, 'FontSize', size_font, 'Color', text_param2{ii, 3}, 'EdgeColor', 'none')
	end
end

% variable text annotations
file_box(1)                 = annotation('textbox', [0.005 0.965 0.16 0.03], 'String', '', 'Color', 'k', 'FontSize', size_font, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Interpreter', 'none', 'LineWidth', 1);
file_box(2)                 = annotation('textbox', [0.345 0.965 0.16 0.03], 'String', '', 'Color', 'k', 'FontSize', size_font, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Interpreter', 'none', 'LineWidth', 1);
status_box					= annotation('textbox', [0.69 0.965 0.30 0.03], 'String', '', 'Color', 'k', 'FontSize', size_font, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Interpreter', 'none', 'LineWidth', 1);

rad_group                   = uibuttongroup('Position', [0.48 0.925 0.05 0.03], 'SelectionChangeFcn', @rad_radio);
uicontrol(fgui, 'Style', 'text', 'parent', rad_group, 'Units', 'normalized', 'Position', [0 0.6 0.9 0.3], 'FontSize', size_font)
rad_check(1)                = uicontrol(fgui, 'Style', 'radio', 'String', 'P', 'Units', 'normalized', 'Position', [0.01 0.1 0.45 0.8], 'parent', rad_group, 'FontSize', size_font, 'HandleVisibility', 'off');
rad_check(2)                = uicontrol(fgui, 'Style', 'radio', 'String', 'I', 'Units', 'normalized', 'Position', [0.51 0.1 0.45 0.8], 'parent', rad_group, 'FontSize', size_font, 'HandleVisibility', 'off');
rad_group.SelectedObject = rad_check(1);

% menus
layer_list(1)				= uicontrol(fgui, 'Style', 'popupmenu', 'String', 'N/A', 'Value', 1, 'Units', 'normalized', 'Position', [0.21 0.955 0.04 0.04], 'FontSize', size_font, 'ForegroundColor', 'k', 'Callback', @pk_select1);
layer_list(2)				= uicontrol(fgui, 'Style', 'popupmenu', 'String', 'N/A', 'Value', 1, 'Units', 'normalized', 'Position', [0.55 0.955 0.04 0.04], 'FontSize', size_font, 'ForegroundColor', 'k', 'Callback', @pk_select2);
int_list                    = uicontrol(fgui, 'Style', 'popupmenu', 'String', 'N/A', 'Value', 1, 'Units', 'Normalized', 'Position', [0.60 0.91 0.12 0.04], 'FontSize', size_font, 'ForegroundColor', 'k', 'Callback', @load_pk2);
intnum_list                 = uicontrol(fgui, 'Style', 'popupmenu', 'String', 'N/A', 'Value', 1, 'Units', 'normalized', 'Position', [0.15 0.91 0.04 0.04], 'FontSize', size_font, 'ForegroundColor', 'k', 'Callback', @change_int);

% check boxes
int_check(1)                = uicontrol(fgui, 'Style', 'checkbox', 'Units', 'normalized', 'Position', [0.20 0.92 0.01 0.02], 'Callback', @show_int1, 'FontSize', size_font, 'Value', 0, 'BackgroundColor', fgui.Color);
int_check(2)                = uicontrol(fgui, 'Style', 'checkbox', 'Units', 'normalized', 'Position', [0.59 0.92 0.01 0.02], 'Callback', @show_int2, 'FontSize', size_font, 'Value', 0, 'BackgroundColor', fgui.Color);
nearest_check               = uicontrol(fgui, 'Style', 'checkbox', 'Units', 'normalized', 'Position', [0.85 0.925 0.01 0.02], 'FontSize', size_font, 'Value', 1, 'BackgroundColor', fgui.Color);
match_check                 = uicontrol(fgui, 'Style', 'checkbox', 'Units', 'normalized', 'Position', [0.895 0.925 0.01 0.02], 'FontSize', size_font, 'Value', 1, 'BackgroundColor', fgui.Color);

linkaxes(ax, 'y')

%% Clear plots

    function clear_plots(src, event)
        delete([p_bed(curr_rad); p_core{curr_rad}(isgraphics(p_core{curr_rad}))'; p_core_name{curr_rad}(isgraphics(p_core_name{curr_rad}))'; p_data(curr_rad); p_pk{curr_rad}(isgraphics(p_pk{curr_rad}))'])
        for ii = 1:2
			delete(p_int1{ii}(isgraphics(p_int1{ii})))
			delete(p_int2{ii}(isgraphics(p_int2{ii})))
        end
        file_box(curr_rad).String = '';
        [int_check(curr_rad).Value] = deal(0);
        [layer_list(curr_rad).String, intnum_list(:).String] = deal('N/A'); [layer_list(curr_rad).Value, intnum_list(:).Value] = deal(1); int_tot2_edit.String = '0';
    end

%% Clear data and picks

    function clear_data(src, event)
        [data_done(curr_rad), pk_done(curr_rad), match_done] ...
                            = deal(false);
        [amp_depth{curr_rad}, colors{curr_rad}, depth_bed{curr_rad}, depth_layer{curr_rad}, depth_ref{curr_rad}, dist_lin{curr_rad}, file_pk_short{curr_rad}, ind_layer{curr_rad}, layer_str{curr_rad}, narrow_ax{curr_rad}, twtt{curr_rad}] ...
							= deal([]);
        [curr_layer(curr_rad), curr_pk(curr_rad), num_layer(curr_rad), num_sample(curr_rad), num_trace(curr_rad), curr_ind_int, ii, ind_x_pk, ind_y_pk, jj, num_int, tmp1, tmp2, tmp3, tmp4, tmp5] ...
                            = deal(0);
        [p_bed(curr_rad), p_core{curr_rad}, p_core_name{curr_rad}, p_data(curr_rad), p_int1{1}, p_int1{2}, p_int2{1}, p_int2{2}, p_pk{curr_rad}] ...
							= deal(gobjects(1));
	end

%% Load picks for this segment

    function load_pk1(src, event)
        [curr_rad, curr_rad_alt] ...
                            = deal(1, 2);
		if (exist([path_pk file_pk{curr_rad}], 'file') == 2) % already have a picked file, want to move onto next one
			tmp1			= dir([path_pk '*.mat']);
			tmp2			= find(contains({tmp1.name}, file_pk{curr_rad}));
			if (tmp2 < length(tmp1))
				file_pk{curr_rad} ...
							= tmp1(tmp2 + 1).name;
			else
				status_box.String = 'Reached end of directory';
				return
			end
		else
			file_pk{curr_rad} ...
							= '';
		end
        [fgui.KeyPressFcn, fgui.WindowButtonDownFcn] = deal('');
        load_pk_breakout
        fgui.KeyPressFcn = @keypress; fgui.WindowButtonDownFcn = @mouse_click;
    end

    function load_pk2(src, event)
        [curr_rad, curr_rad_alt] ...
                            = deal(2, 1);
		if ~strcmp(int_list.String, 'none')
			file_pk{curr_rad} ...
							= ['Data_' int_list.String{int_list.Value} '_pk.mat'];
		else
			status_box.String = 'No later intersections!';
			return
		end
        [fgui.KeyPressFcn, fgui.WindowButtonDownFcn] = deal('');
        load_pk_breakout
        fgui.KeyPressFcn = @keypress; fgui.WindowButtonDownFcn = @mouse_click;
	end

    function load_pk_breakout(src, event)
		
        rad_group.SelectedObject = rad_check(curr_rad);
        
        if ((curr_rad == 2) && ~pk_done(1))
            status_box.String = 'Load primary segment before intersecting segment.';
            return
        end
        
        % dialog box to choose picks file to load
		if (exist([path_pk file_pk{curr_rad}], 'file') ~= 2)
			if ~isempty(path_pk)
            	[file_pk{curr_rad}, path_pk] ...
 							= uigetfile('*.mat', 'Load merged picks:', path_pk);
        	elseif ~isempty(path_data)
            	[file_pk{curr_rad}, path_pk] ...
							= uigetfile('*.mat', 'Load merged picks:', path_data);
        	else
            	[file_pk{curr_rad}, path_pk] ...
							= uigetfile('*.mat', 'Load merged picks:');
			end
		end
		
		if isnumeric(file_pk{curr_rad})
            file_pk{curr_rad} ...
							= '';
		end
        if isempty(file_pk{curr_rad})
            status_box.String = 'No picks loaded.';
            return
        end
		
		% get rid of intersecting picks/data if moving on to another primary segment
		if (curr_rad == 1)
			curr_rad		= 2;
			clear_plots
			clear_data
			curr_rad		= 1;
		end
		clear_plots
		clear_data
		
        file_pk_short{curr_rad} ...
							= file_pk{curr_rad}(6:(end - 7));
		
        % current index of subsegment within complete list
		curr_pk(curr_rad)	= find(ismember(name_pk, file_pk_short(curr_rad)));
        
        % load picks file
        status_box.String = ['Loading ' file_pk{curr_rad}(1:(end - 4)) '...'];
        pause(0.05)
        try
            tmp1            = load([path_pk file_pk{curr_rad}], 'pk');
            tmp1			= tmp1.pk;
			[dist_lin{curr_rad}, file_data{curr_rad}, num_layer(curr_rad), num_trace(curr_rad)] ...
							= deal(tmp1.dist_lin, tmp1.file_cat, tmp1.num_layer, tmp1.num_trace);
			depth_layer{curr_rad} ...
							= reshape([tmp1.layer(:).depth], num_trace(curr_rad), num_layer(curr_rad))';
			depth_bed{curr_rad}	...
							= (speed_ice / 2) .* (tmp1.twtt_bed - tmp1.twtt_surf);
			tmp1			= 0;
        catch % give up, force restart
            status_box.String = [file_pk{curr_rad} ' does not contain a pk structure. Try again.'];
            return
        end
        
		file_box(curr_rad).String = file_pk_short{curr_rad};
        
        % list matches
        tmp1                = find(ismember(layer_match_list(:, 1), curr_pk(curr_rad), 'rows'));
        tmp2                = find(ismember(layer_match_list(:, 3), curr_pk(curr_rad), 'rows'));
        if (~isempty(tmp1) || ~isempty(tmp2))
            ind_layer{curr_rad} ...
                            = unique([layer_match_list(tmp1, 2:4); layer_match_list(tmp2, [4 1:2])], 'rows');
        else
            ind_layer{curr_rad} ...
                            = NaN(0, 3);
        end
        
        % find intersections for primary segment
        switch curr_rad
            
            case 1 % primary
                
				if isempty(tmp_list)
					tmp1    = find((int_all(:, 3) == curr_pk(1)));
					tmp2    = find((int_all(:, 1) == curr_pk(1)));
                	tmp3    = cell((length(tmp1) + length(tmp2)), 1);
					for ii = 1:length(tmp1)
						tmp3{ii} ...
							= name_pk{int_all(tmp1(ii), 1)};
					end
					for ii = (length(tmp1) + 1):(length(tmp1) + length(tmp2))
						tmp3{ii} ...
							= name_pk{int_all(tmp2(ii - length(tmp1)), 3)};
					end
					tmp2	= tmp3;
				else
					tmp2	= cell(1, length(tmp_list));
					for ii = 1:length(tmp2)
						tmp2{ii} ...
							= tmp_list{ii}(6:(end - 7));
					end
				end
				if ~isempty(unique(tmp2))
					int_list.String = unique(tmp2);
					int_tot1_edit.String = length(unique(tmp2));
				else
					int_list.String = 'none';
					int_count_edit.String = '0 /';
				end
				int_list.Value = 1; int_count_edit.String = [num2str(int_list.Value) ' /'];
				
            case 2 % figure out intersections
				
                tmp1        = find(ismember(int_all(:, [1 3]), curr_pk(1:2), 'rows'));
                tmp2        = find(ismember(int_all(:, [3 1]), curr_pk(1:2), 'rows'));
                curr_ind_int= [int_all(tmp1, [2 4]); int_all(tmp2, [4 2])];
				
				% check intersections for layers
				tmp5		= [];
				for ii = 1:size(curr_ind_int, 1)
					[tmp1, tmp2, tmp3, tmp4] ...
                            = deal((dist_lin{1}(curr_ind_int(ii, 1)) - dist_int_half), (dist_lin{1}(curr_ind_int(ii, 1)) + dist_int_half), ...
								   (dist_lin{2}(curr_ind_int(ii, 2)) - dist_int_half), (dist_lin{2}(curr_ind_int(ii, 2)) + dist_int_half));
					if (isempty(find(~isnan(depth_layer{1}(:, interp1(dist_lin{1}, 1:num_trace(1), tmp1, 'nearest', 'extrap'):interp1(dist_lin{1}, 1:num_trace(1), tmp2, 'nearest', 'extrap'))), 1)) || ...
						isempty(find(~isnan(depth_layer{2}(:, interp1(dist_lin{2}, 1:num_trace(2), tmp3, 'nearest', 'extrap'):interp1(dist_lin{2}, 1:num_trace(2), tmp4, 'nearest', 'extrap'))), 1)))
						tmp5=[tmp5 ii]; %#ok<AGROW>
					end
				end
				
				curr_ind_int= curr_ind_int(setdiff(1:size(curr_ind_int, 1), tmp5), :);
                num_int     = size(curr_ind_int, 1);
				
				if ~num_int
					status_box.String = 'No traced intersections worth comparing.';
					pause(0.5)
					if (int_list.Value < length(int_list.String))
						int_list.Value = int_list.Value + 1;
						int_count_edit.String = [num2str(int_list.Value) ' /'];
						load_pk2
						if num_int
							load_data
						end
					else
						pk_save
						status_box.String = 'End of intersecting list.';
					end
					return
				else
					intnum_list.String = num2cell(1:num_int); intnum_list.Value = 1; int_tot2_edit.String = num2str(num_int);
				end
        end
        
        % display merged picks
        colors{curr_rad}    = repmat(colors_def, ceil(num_layer(curr_rad) / size(colors_def, 1)), 1); % extend predefined color pattern
        colors{curr_rad}    = colors{curr_rad}(1:num_layer(curr_rad), :);
        p_pk{curr_rad}		= gobjects(1, num_layer(curr_rad));
        layer_str{curr_rad} = num2cell(1:num_layer(curr_rad));
        curr_layer(curr_rad) = 1;
        layer_list(curr_rad).String = layer_str{curr_rad}; layer_list(curr_rad).Value = 1;
		
        % display picks, surface and bed
        for ii = 1:num_layer(curr_rad) %#ok<*FXUP>
			if ~isempty(find(~isnan(depth_layer{curr_rad}(ii, :)), 1))
				p_pk{curr_rad}(ii) ...
                            = line(ax(curr_rad), (1e-3 .* dist_lin{curr_rad}(~isnan(depth_layer{curr_rad}(ii, :)))), depth_layer{curr_rad}(ii, ~isnan(depth_layer{curr_rad}(ii, :))), ...
								   'Marker', '.', 'Color', colors{curr_rad}(ii, :), 'MarkerSize', 12, 'LineStyle', 'none', 'Visible', 'off');
			else
				p_pk{curr_rad}(ii) ...
                            = line(ax(curr_rad), NaN, NaN, 'Marker', '.', 'Color', colors{curr_rad}(ii, :), 'MarkerSize', 12, 'LineStyle', 'none', 'Visible', 'off');
			end
        end
        
		switch curr_rad
			case 1
				depth_max	= max(depth_bed{curr_rad}, [], 'all', 'omitnan');
			case 2
				depth_max	= max([depth_bed{:}], [], 'all', 'omitnan');
		end
		if (isinf(depth_max) || isnan(depth_max))
			depth_max		= depth_max_ref;
		end
		[z_max_slide.Max, z_max_slide.Value] = deal(depth_max_ref);
		z_min_slide.Max = depth_max_ref;
		
		p_bed(curr_rad)		= line(ax(curr_rad), (1e-3 .* dist_lin{curr_rad}), depth_bed{curr_rad}, 'Marker', '.', 'Color', 'm', 'MarkerSize', 12, 'LineStyle', 'none', 'Visible', 'off');
		
        % get new min/max limits
        [dist_min_ref(curr_rad), dist_max_ref(curr_rad), dist_min(curr_rad), dist_max(curr_rad)] ...
                            = deal(min(dist_lin{curr_rad}), max(dist_lin{curr_rad}), min(dist_lin{curr_rad}), max(dist_lin{curr_rad}));
        
		dist_min_slide(curr_rad).Min = 1e-3 * dist_min_ref(curr_rad); dist_min_slide(curr_rad).Max = 1e-3 * dist_max_ref(curr_rad); dist_min_slide(curr_rad).Value = 1e-3 * dist_min_ref(curr_rad);
		dist_max_slide(curr_rad).Min = 1e-3 * dist_min_ref(curr_rad); dist_max_slide(curr_rad).Max = 1e-3 * dist_max_ref(curr_rad); dist_max_slide(curr_rad).Value = 1e-3 * dist_max_ref(curr_rad);
		dist_min_edit(curr_rad).String = sprintf('%4.1f', (1e-3 * dist_min_ref(curr_rad)));
		dist_max_edit(curr_rad).String = sprintf('%4.1f', (1e-3 * dist_max_ref(curr_rad)));
        
        [pk_check(curr_rad), pk_done(curr_rad)] ...
							= deal(true);
		
        % intersection plots
        if (curr_rad == 2)
			for ii = 1:2
				p_int1{ii}	= gobjects(1, num_int);
			end
            p_int2{1}		= gobjects(1, num_layer(2));
			p_int2{2}		= gobjects(1, num_layer(1));
            for ii = 1:num_int
				for jj = 1:2
					p_int1{jj}(ii) = line(ax(jj), (1e-3 .* dist_lin{jj}(curr_ind_int(ii, [jj jj]))), [depth_min_ref depth_max_ref], 'LineStyle', '--', 'Color', 'm', 'LineWidth', 2, 'Visible', 'off');
				end
            end
            for ii = 1:num_layer(2)
                p_int2{1}(ii) = line(ax(1), (1e-3 * dist_lin{1}(curr_ind_int(:, 1))), depth_layer{2}(ii, curr_ind_int(:, 2)), 'Marker', 'o', 'Color', 'k', 'LineWidth', 1, 'MarkerSize', 12, 'MarkerFaceColor', colors{2}(ii, :), 'LineStyle', 'none', 'Visible', 'off');
            end
			for ii = 1:num_layer(1)
                p_int2{2}(ii) = line(ax(2), (1e-3 * dist_lin{2}(curr_ind_int(:, 2))), depth_layer{1}(ii, curr_ind_int(:, 1)), 'Marker', 'o', 'Color', 'k', 'LineWidth', 1, 'MarkerSize', 12, 'MarkerFaceColor', colors{1}(ii, :), 'LineStyle', 'none', 'Visible', 'off');
			end
            if ~isempty(ind_layer{1})
                [p_int2{2}(ind_layer{1}(:, 1)).Marker] = deal('s');
            end
			if ~isempty(ind_layer{2})
                [p_int2{1}(ind_layer{2}(:, 1)).Marker] = deal('s');
			end
            for ii = find(ind_layer{1}(:, 2) == curr_pk(2))' % match to current segment
                % if (length(find((ind_layer{1}(:, end) == ind_layer{1}(ii, end)))) > 1)
                %     tmp1= find((ind_layer{1}(:, end) == ind_layer{1}(ii, end)));
                %     colors{1}(ind_layer{1}(tmp1(2:end), 1), :) ...
                %         = repmat(colors{1}(ind_layer{1}(tmp1(1), 1), :), length(tmp1(2:end)), 1);
                %     p_pk{1}(ind_layer{1}(ii, 1)).Color = colors{1}(ind_layer{1}(tmp1(1), 1), :);
                % end
                [p_int2{2}(ind_layer{1}(ii, 1)).Marker, p_int2{1}(ind_layer{1}(ii, 3)).Marker] = deal('^');
				[p_int2{2}(ind_layer{1}(ii, 1)).MarkerFaceColor, p_int2{1}(ind_layer{1}(ii, 3)).MarkerFaceColor] = deal(colors{1}(ind_layer{1}(ii, 1), :));
                p_pk{2}(ind_layer{1}(ii, 3)).Color = colors{1}(ind_layer{1}(ii, 1), :);
            end
            [int_check(:).Value] = deal(1);
        end
		
        ax(curr_rad).XLim = 1e-3 .* [dist_min(curr_rad) dist_max(curr_rad)];
		ax(curr_rad).YLim = [depth_min depth_max];
		
        if all(pk_done)
            intnum_list.Value = 1;
            change_int
        end
        status_box.String = ['Loaded ' file_pk_short{curr_rad} '.'];
		
		pause(0.05)
		
		if all(pk_done)
			curr_rad		= 2;
			[nearest_check.Value, match_check.Value] = deal(1); int_count_edit.String = [num2str(int_list.Value) ' /'];
		end
		load_data
    end

%% Load radar data

    function load_data(src, event) %#ok<*INUSD>
        
        rad_group.SelectedObject = rad_check(curr_rad);
		
        if (exist([path_data file_data{curr_rad} '.mat'], 'file') ~= 2)
			[file_data{curr_rad}, path_data] ...
                            = uigetfile('*.mat', 'Load radar data:');
        end
        
		% quick checks
		if isnumeric(file_data{curr_rad})
            [file_data{curr_rad}, path_data] ...
                            = deal('');
		end
        if isempty(file_data{curr_rad})
            status_box.String = 'No radar data loaded.';
            return
        end
		
        status_box.String = ['Loading ' file_data{curr_rad} '...'];
        pause(0.05)
        tmp1            = load([path_data file_data{curr_rad}], 'amp', 'twtt', 'twtt_surf', 'num_sample');
		if ~isfield(tmp1, 'amp')
            status_box.String = [file_data{curr_rad} ' does not contain expected fields. Try again.'];
            return
		end
        
        [num_sample(curr_rad), twtt{curr_rad}] ...
							= deal(tmp1.num_sample, tmp1.twtt);
		tmp2				= diff(twtt{curr_rad}(1:2));
        depth_ref{curr_rad} = (speed_ice / 2) .* (0:tmp2:((num_sample(curr_rad) - 1) * tmp2))'; % simple monotonically increasing depth vector
		
		tmp3				= interp1(twtt{curr_rad}, 1:num_sample(curr_rad), tmp1.twtt_surf, 'nearest', 'extrap'); % surface traveltime indices
        if ~isempty(find(isnan(tmp3), 1))
            tmp3(isnan(tmp3)) ...
							= round(interp1(find(~isnan(tmp3)), tmp3(~isnan(tmp3)), find(isnan(tmp3)), 'linear', 'extrap'));
        end
        tmp3(tmp3 < 1)		= 1;
        tmp3(tmp3 > num_sample(curr_rad)) ...
                            = num_sample(curr_rad);
		
    	amp_depth{curr_rad} = NaN(num_sample(curr_rad), num_trace(curr_rad), 'single');
    	for ii = find(~isnan(tmp3))
        	amp_depth{curr_rad}(1:(num_sample(curr_rad) - tmp3(ii) + 1), ii) ...
							= tmp1.amp(tmp3(ii):num_sample(curr_rad), ii); % shift data up to surface
    	end
		
		% done with data file
		tmp1				= 0;
		
        % assign traveltime and distance reference values/sliders based on data
        [db_min_ref(curr_rad), db_max_ref(curr_rad), db_min(curr_rad), db_max(curr_rad)] ...
                            = deal(min(amp_depth{curr_rad}, [], 'all', 'omitnan'), max(amp_depth{curr_rad}, [], 'all', 'omitnan'), min(amp_depth{curr_rad}, [], 'all', 'omitnan'), max(amp_depth{curr_rad}, [], 'all', 'omitnan'));
        cb_min_slide(curr_rad).Min = db_min_ref(curr_rad); cb_min_slide(curr_rad).Max = db_max_ref(curr_rad); cb_min_slide(curr_rad).Value = db_min_ref(curr_rad);
        cb_max_slide(curr_rad).Min = db_min_ref(curr_rad); cb_max_slide(curr_rad).Max = db_max_ref(curr_rad); cb_max_slide(curr_rad).Value = db_max_ref(curr_rad);
        cb_min_edit(curr_rad).String = sprintf('%3.0f', db_min_ref(curr_rad));
        cb_max_edit(curr_rad).String = sprintf('%3.0f', db_max_ref(curr_rad));
        
        % plot data
        [data_check(curr_rad), data_done(curr_rad)] ...
							= deal(true);
        plot_data
		if (curr_rad == 2)
			change_int
			show_int1
			show_int2
		else
			load_pk2
		end
        status_box.String = 'Segment radar data loaded.';
    end

%% Select current layer

    function pk_select1(src, event)
        [curr_rad, curr_rad_alt] ...
                        = deal(1, 2);
        pk_select
    end

    function pk_select2(src, event)
        [curr_rad, curr_rad_alt] ...
                        = deal(2, 1);
        pk_select
	end

    function pk_select(src, event)
		
        rad_group.SelectedObject = rad_check(curr_rad);
        curr_layer(curr_rad)= layer_list(curr_rad).Value;
		
		if pk_done(curr_rad)
            layer_list(curr_rad).Value = curr_layer(curr_rad);
            [p_pk{curr_rad}(isgraphics(p_pk{curr_rad})).MarkerSize] = deal(12);
            [p_int2{curr_rad_alt}(isgraphics(p_int2{curr_rad_alt})).MarkerSize] = deal(12);
            if isgraphics(p_pk{curr_rad}(curr_layer(curr_rad)))
                p_pk{curr_rad}(curr_layer(curr_rad)).MarkerSize = 24;
            end
            if all(pk_done)
                if isgraphics(p_int2{curr_rad_alt}(curr_layer(curr_rad)))
                    p_int2{curr_rad_alt}(curr_layer(curr_rad)).MarkerSize = 18;
                end
            end
            [p_pk{curr_rad}(isgraphics(p_pk{curr_rad})).MarkerSize] = deal(12);
            p_pk{curr_rad}(curr_layer(curr_rad)).MarkerSize = 24;
		end
		
        if ischar(tmp1)
            if strcmp(tmp1, 'reselect')
                return
            end
        end
        tmp1                = '';
		if (all(pk_done) && ~isempty(ind_layer{2}) && match_check.Value)
            switch curr_rad
                case 1
                    if ~isempty(find(((ind_layer{2}(:, 3) == curr_layer(1)) & (ind_layer{2}(:, 2) == curr_pk(1))), 1))
                        curr_layer(2) = ind_layer{2}(find(((ind_layer{2}(:, 3) == curr_layer(1)) & (ind_layer{2}(:, 2) == curr_pk(1))), 1), 1);
                        for ii = 1:2
                            [p_int2{ii}(isgraphics(p_int2{ii})).MarkerSize] = deal(12);
                        end
                        if isgraphics(p_int2{1}(curr_layer(2)))
                            p_int2{1}(curr_layer(2)).MarkerSize = 18;
                        end
                        if isgraphics(p_int2{2}(curr_layer(1)))
                            p_int2{2}(curr_layer(1)).MarkerSize = 18;
                        end
                        [p_pk{2}(isgraphics(p_pk{2})).MarkerSize] = deal(12);
                        if isgraphics(p_pk{2}(curr_layer(2)))
                            p_pk{2}(curr_layer(2)).MarkerSize = 24;
                        end
                        tmp1 = 'matched';
                    end
                    layer_list(2).Value = curr_layer(2);
                case 2
                    if ~isempty(find(((ind_layer{2}(:, 1) == curr_layer(2)) & (ind_layer{2}(:, 2) == curr_pk(1))), 1))
                        curr_layer(1) = ind_layer{2}(find(((ind_layer{2}(:, 1) == curr_layer(2)) & (ind_layer{2}(:, 2) == curr_pk(1))), 1), 3);
						for ii = 1:2
							[p_int2{ii}(isgraphics(p_int2{ii})).MarkerSize] = deal(12);
						end
                        if isgraphics(p_int2{1}(curr_layer(2)))
                            p_int2{1}(curr_layer(2)).MarkerSize = 18;
                        end
                        if isgraphics(p_int2{2}(curr_layer(1)))
                            p_int2{2}(curr_layer(1)).MarkerSize = 18;
                        end
                        [p_pk{1}(isgraphics(p_pk{1})).MarkerSize] = deal(12);
                        if isgraphics(p_pk{1}(curr_layer(1)))
                            p_pk{1}(curr_layer(1)).MarkerSize = 24;
                        end
                        tmp1 = 'matched';
                    end
                    layer_list(1).Value = curr_layer(1);
            end
		end

        if (nearest_check.Value && all(pk_done) && isempty(tmp1))
			if any(~isnan(depth_layer{curr_rad_alt}(:, curr_ind_int(curr_int, curr_rad_alt))))
                [~, curr_layer(curr_rad_alt)] ...
                            = min(abs(depth_layer{curr_rad}(curr_layer(curr_rad), curr_ind_int(curr_int, curr_rad)) - depth_layer{curr_rad_alt}(:, curr_ind_int(curr_int, curr_rad_alt))));
                [curr_rad, curr_rad_alt] ...
                            = deal(curr_rad_alt, curr_rad);
                layer_list(curr_rad).Value = curr_layer(curr_rad);
                tmp1        = 'reselect';
                pk_select
                [curr_rad, curr_rad_alt] ...
                            = deal(curr_rad_alt, curr_rad);
        		rad_group.SelectedObject = rad_check(curr_rad);
        		curr_layer(curr_rad) ...
							= layer_list(curr_rad).Value;
			end
        else
            status_box.String = deal(['Layer #' num2str(curr_layer(curr_rad)) ' selected.']);
        end
    end

%% Focus on a layer

    function pk_focus(src, event)
        if isempty(find(~isnan(depth_layer{curr_rad}(curr_layer(curr_rad), :)), 1))
            status_box.String = 'Cannot focus on a layer that is hidden.';
            return
        end
        ax(curr_rad).XLim = 1e-3 .* [dist_lin{curr_rad}(find(~isnan(depth_layer{curr_rad}(curr_layer(curr_rad), :)), 1)) dist_lin{curr_rad}(find(~isnan(depth_layer{curr_rad}(curr_layer(curr_rad), :)), 1, 'last'))];
        tmp1                = 1e3 .* ax(curr_rad).XLim;
        [tmp1(1), tmp1(2)]  = deal((tmp1(1) - diff(tmp1)), (tmp1(2) + diff(tmp1)));
        if (tmp1(1) < dist_min_ref(curr_rad))
            tmp1(1)         = dist_min_ref(curr_rad);
        end
        if (tmp1(2) > dist_max_ref(curr_rad))
            tmp1(2)         = dist_max_ref(curr_rad);
        end
        ax(curr_rad).XLim = 1e-3 .* tmp1;
        [dist_min(curr_rad), dist_max(curr_rad)] ...
                            = deal(tmp1(1), tmp1(2));
        if (dist_min(curr_rad) < (1e3 * dist_min_slide(curr_rad).Min))
            dist_min_slide(curr_rad).Value = dist_min_slide(curr_rad).Min;
        elseif (dist_min(curr_rad) > dist_min_slide(curr_rad).Max)
            dist_min_slide(curr_rad).Value = dist_min_slide(curr_rad).Max;
        else
            dist_min_slide(curr_rad).Value = 1e-3 * dist_min(curr_rad);
        end
        if (dist_max(curr_rad) < (1e3 * dist_max_slide(curr_rad).Min))
            dist_max_slide(curr_rad).Value = dist_max_slide(curr_rad).Min;
        elseif (dist_max(curr_rad) > (1e3 * dist_max_slide(curr_rad).Max))
            dist_max_slide(curr_rad).Value = dist_max_slide(curr_rad).Max;
        else
            dist_max_slide(curr_rad).Value = 1e-3 * dist_max(curr_rad);
        end
        dist_min_edit(curr_rad).String = sprintf('%3.1f', (1e-3 * dist_min(curr_rad)));
        dist_max_edit(curr_rad).String = sprintf('%3.1f', (1e-3 * dist_max(curr_rad)));
        ax(curr_rad).YLim = [min(depth_layer{curr_rad}(curr_layer(curr_rad), :), [], 'omitnan') max(depth_layer{curr_rad}(curr_layer(curr_rad), :), [], 'omitnan')];
        tmp1				= ax(curr_rad).YLim;
        [tmp1(1), tmp1(2)]	= deal((tmp1(1) - diff(tmp1)), (tmp1(2) + diff(tmp1)));
        if (tmp1(1) < depth_min_ref)
            tmp1(1)			= depth_min_ref;
        end
        if (tmp1(2) > depth_max_ref)
            tmp1(2)			= depth_max_ref;
        end
        ax(curr_rad).YLim = tmp1;
        [depth_min, depth_max] ...
							= deal(tmp1(1), tmp1(2));
        if ((depth_max_ref - (depth_max - depth_min_ref)) < z_min_slide.Min)
            z_min_slide.Value = z_min_slide.Min;
        elseif ((depth_max_ref - (depth_max - depth_min_ref)) > z_min_slide.Max)
            z_min_slide.Value = z_min_slide.Max;
        else
            z_min_slide.Value = (depth_max_ref - (depth_max - depth_min_ref));
        end
        if ((depth_max_ref - (depth_min - depth_min_ref)) < z_max_slide.Min)
            z_max_slide.Value = z_max_slide.Min;
        elseif ((depth_max_ref - (depth_min - depth_min_ref)) > z_max_slide.Max)
            z_max_slide.Value = z_max_slide.Max;
        else
            z_max_slide.Value = (depth_max_ref - (depth_min - depth_min_ref));
        end
        z_min_edit.String = sprintf('%4.0f', depth_max);
        z_max_edit.String = sprintf('%4.0f', depth_min);
        narrow_cb
        status_box.String = ['Focused on layer #' num2str(curr_layer(curr_rad)) '.'];
    end

%% Switch to previous layer in list

    function pk_last(src, event)
        rad_group.SelectedObject = rad_check(curr_rad);
        if (pk_done(curr_rad) && (curr_layer(curr_rad) > 1))
            curr_layer(curr_rad) ...
                            = curr_layer(curr_rad) - 1;
            layer_list(curr_rad).Value = curr_layer(curr_rad);
            pk_select
        end
    end

%% Switch to next layer in the list

    function pk_next(src, event)
        rad_group.SelectedObject = rad_check(curr_rad);
        if (pk_done(curr_rad) && (curr_layer(curr_rad) < num_layer(curr_rad)))
            curr_layer(curr_rad) ...
                            = curr_layer(curr_rad) + 1;
            layer_list(curr_rad).Value = curr_layer(curr_rad);
            pk_select
        end
    end

%% Match two intersecting layers

    function pk_match(src, event)
        
        [fgui.KeyPressFcn, fgui.WindowButtonDownFcn] = deal('');
        
        if ~all(pk_done)
            status_box.String = 'Picks must be loaded for both primary and intersecting segments.';
            fgui.KeyPressFcn = @keypress; fgui.WindowButtonDownFcn = @mouse_click;
            return
        end
        
        % check if the intersecting layer is already matched to a primary layer
        if ~isempty(ind_layer{1})
            if ~isempty(find(ismember(ind_layer{1}(:, 2:3), [curr_pk(2) curr_layer(2)], 'rows'), 1))
                tmp1        = find(ismember(ind_layer{1}(:, 2:3), [curr_pk(2) curr_layer(2)], 'rows'));
                tmp2        = colors{1}(ind_layer{1}(tmp1(1), 1), :);
            else
                tmp2        = colors{1}(curr_layer(1), :);
            end
        else
            tmp2            = colors{1}(curr_layer(1), :);
        end
        
        % check for existing match
        if ~isempty(ind_layer{2})
            if ~isempty(find(ismember(ind_layer{2}, [curr_layer(2) curr_pk(1) curr_layer(1)], 'rows'), 1))
                status_box.String = 'This layer pair is already matched.';
                fgui.KeyPressFcn = @keypress; fgui.WindowButtonDownFcn = @mouse_click;
                return
            end
        end
        
        % colorize then verify
		if isgraphics(p_pk{2}(curr_layer(2)))
			p_pk{2}(curr_layer(2)).Color = tmp2;
		end
        if isgraphics(p_int2{1}(curr_layer(2)))
            p_int2{1}(curr_layer(2)).MarkerFaceColor = tmp2;
        end
        if isgraphics(p_int2{2}(curr_layer(1)))
            p_int2{2}(curr_layer(1)).MarkerFaceColor = tmp2;
        end
        
        status_box.String = ['Matching primary segment layer #' num2str(curr_layer(1)) ' with intersecting segment layer #' num2str(curr_layer(2)) '...'];
        
        pause(0.05)
        
        % reassign primary color to first primary color if a primary layer is matched to multiple intersecting segments
        if ~isempty(ind_layer{1})
            if ~isempty(find(ismember(ind_layer{1}(:, 2:3), [curr_pk(2) curr_layer(2)], 'rows'), 1))
                colors{1}(ind_layer{1}(tmp1, 1), :) ...
                            = repmat(tmp2, length(tmp1), 1);
				[p_pk{1}(ind_layer{1}(tmp1, 1)).Color, p_int2{2}(ind_layer{1}(tmp1, 1)).MarkerFaceColor] = deal(tmp2);
            end
        end
        
        % add each layer to the other's list
        ind_layer{1}		= [ind_layer{1}; [curr_layer(1) curr_pk(2) curr_layer(2)]];
        ind_layer{2}		= [ind_layer{2}; [curr_layer(2) curr_pk(1) curr_layer(1)]];
        
        if isgraphics(p_int2{1}(curr_layer(2)))
            p_int2{1}(curr_layer(2)).Marker = '^';
        end
        if isgraphics(p_int2{2}(curr_layer(1)))
            p_int2{2}(curr_layer(1)).Marker = '^';
        end
        
		match_done			= true;
        status_box.String = ['Intersecting segment layer #' num2str(curr_layer(2)) ' matched to primary segment layer #' num2str(curr_layer(1)) '.'];
        fgui.KeyPressFcn = @keypress; fgui.WindowButtonDownFcn = @mouse_click;
    end

%% Unmatch two intersecting layers

    function pk_unmatch(src, event)
        
        [fgui.KeyPressFcn, fgui.WindowButtonDownFcn] = deal('');
		
        if ~all(pk_done)
            status_box.String = 'Picks must be loaded for both primary and intersecting segments.';
            fgui.KeyPressFcn = @keypress; fgui.WindowButtonDownFcn = @mouse_click;
            return
        end
        
        status_box.String = 'Unmatch current pair? (Y: yes; otherwise: cancel)...';
        waitforbuttonpress
        
        if ~strcmpi(fgui.CurrentCharacter, 'Y')
            status_box.String = 'Unmatching canceled.';
            fgui.KeyPressFcn = @keypress; fgui.WindowButtonDownFcn = @mouse_click;
            return
        end
        
        status_box.String = ['Unmatching primary segment layer #' num2str(curr_layer(1)) ' from intersecting segment layer # ' num2str(curr_layer(2)) '...'];
        pause(0.05)
        
        if (~isempty(find(ismember(ind_layer{1}(:, [1 3]), [curr_layer(1) curr_layer(2)], 'rows'), 1)) && ~isempty(find(ismember(ind_layer{2}(:, [1 3]), [curr_layer(2) curr_layer(1)], 'rows'), 1)))
            ind_layer{1} = ind_layer{1}(setdiff(1:size(ind_layer{1}, 1), find(((ind_layer{1}(:, 1) == curr_layer(1)) & (ind_layer{1}(:, 3) == curr_layer(2))))), :);
            ind_layer{2} = ind_layer{2}(setdiff(1:size(ind_layer{2}, 1), find(((ind_layer{2}(:, 1) == curr_layer(2)) & (ind_layer{2}(:, 3) == curr_layer(1))))), :);
        else
            status_box.String = 'Current layer pair not matched. Unmatching canceled.';
            fgui.KeyPressFcn = @keypress; fgui.WindowButtonDownFcn = @mouse_click;
            return
        end
        
        if isgraphics(p_pk{2}(curr_layer(2)))
            p_pk{2}(curr_layer(2)).Color = colors{2}(curr_layer(2), :);
        end
        if isgraphics(p_int2{1}(curr_layer(2)))
            p_int2{1}(curr_layer(2)).MarkerFaceColor = colors{2}(curr_layer(2), :);
        end
        if isgraphics(p_int2{2}(curr_layer(1)))
            p_int2{2}(curr_layer(1)).MarkerFaceColor = colors{1}(curr_layer(1), :);
        end
        
        if ~isempty(find((ind_layer{1}(:, 1) == curr_layer(1)), 1))
            if isgraphics(p_int2{2}(curr_layer(1)))
                p_int2{2}(curr_layer(1)).Marker = 's';
            end
        else
            if isgraphics(p_int2{2}(curr_layer(1)))
                p_int2{2}(curr_layer(1)).Marker = 'o';
            end
        end
        if isempty(find((ind_layer{1}(:, 3) == curr_layer(2)), 1))
            if isgraphics(p_int2{1}(curr_layer(2)))
                p_int2{1}(curr_layer(2)).Marker = 'o';
            end
        end
        if ~isempty(find((ind_layer{2}(:, 1) == curr_layer(2)), 1))
            if isgraphics(p_int2{1}(curr_layer(2)))
                p_int2{1}(curr_layer(2)).Marker = 's';
            end
        else
            if isgraphics(p_int2{1}(curr_layer(2)))
                p_int2{1}(curr_layer(2)).Marker = 'o';
            end
        end
        if isempty(find((ind_layer{2}(:, 3) == curr_layer(1)), 1))
            if isgraphics(p_int2{1}(curr_layer(2)))
                p_int2{1}(curr_layer(2)).Marker = 'o';
            end
        end
        
		match_done			= true;
        status_box.String = ['Intersecting layer # ' num2str(curr_layer(2)) ' unmatched from primary layer #' num2str(curr_layer(1)) '.'];
        fgui.KeyPressFcn = @keypress; fgui.WindowButtonDownFcn = @mouse_click;
    end

%% Save intersection assignments

    function pk_save(src, event)
        
		if ~match_done
			return
		end
		
        [fgui.KeyPressFcn, fgui.WindowButtonDownFcn] = deal('');
        
        % want everything done before saving
        if ~all(pk_done)
            status_box.String = 'Intersecting picks not loaded yet.';
            fgui.KeyPressFcn = @keypress; fgui.WindowButtonDownFcn = @mouse_click;
            return
        end
        
        if (~ischar(file_pk{1}) || ~ischar(file_pk{2}))
            status_box.String = 'Picks'' locations not all assigned. Saving canceled.';
            fgui.KeyPressFcn = @keypress; fgui.WindowButtonDownFcn = @mouse_click;
            return
        end
        
        status_box.String = 'Preparing matches for saving...';
        pause(0.05)
        
        % remove earlier set for these sub-segments
        layer_match_list    = layer_match_list((~ismember(layer_match_list(:, 1), curr_pk(1)) & ~ismember(layer_match_list(:, 3), curr_pk(1)) & ~ismember(layer_match_list(:, 1), curr_pk(2)) & ~ismember(layer_match_list(:, 3), curr_pk(2))), :);
        
        % add new match set for these sub-segments, then reorder and remove any repeated rows
        layer_match_list    = [layer_match_list; repmat(curr_pk(1), size(ind_layer{1}, 1), 1) ind_layer{1}; repmat(curr_pk(2), size(ind_layer{2}, 1), 1) ind_layer{2}];
		
		tmp1				= find(layer_match_list(:, 3) < layer_match_list(:, 1));
		layer_match_list(tmp1, :) ...
							= layer_match_list(tmp1, [3 4 1 2]);
		
        layer_match_list    = sortrows(unique(layer_match_list, 'rows'), [1 3 2 4]);


        status_box.String = 'Saving matching layer list...';
        pause(0.05)
        try
            save([path_match file_match], '-v7.3', 'layer_match_list')
            status_box.String = 'Saved matching layer list.';
        catch
            status_box.String = 'MATCHING LAYER LIST DID NOT SAVE. Try saving again shortly. Do not perform any other operation.';
            fgui.KeyPressFcn = @keypress; fgui.WindowButtonDownFcn = @mouse_click;
            return
        end
        
        % move on to next intersecting file
        if (int_list.Value < length(int_list.String))
            int_list.Value = int_list.Value + 1;
			int_count_edit.String = [num2str(int_list.Value) ' /'];			
            load_pk2
			load_data
        else
            fgui.KeyPressFcn = @keypress; fgui.WindowButtonDownFcn = @mouse_click;
        end
    end

%% Update minimum dB/z/dist

	function slide_db_min1(src, event)
        [curr_rad, curr_rad_alt] ...
                            = deal(1, 2);
        slide_db_min
    end

	function slide_db_min2(src, event)
        [curr_rad, curr_rad_alt] ...
                            = deal(2, 1);
        slide_db_min
    end

    function slide_db_min(src, event)
        if (cb_min_slide(curr_rad).Value < db_max(curr_rad))
            db_min(curr_rad) = cb_min_slide(curr_rad).Value;
            cb_min_edit(curr_rad).String = sprintf('%3.0f', db_min(curr_rad));
            update_db_range
        else
            if (db_min(curr_rad) < cb_min_slide(curr_rad).Min)
                cb_min_slide(curr_rad).Value = cb_min_slide(curr_rad).Min;
            else
                cb_min_slide(curr_rad).Value = db_min(curr_rad);
            end
        end
        cb_min_slide(curr_rad).Enable = 'off';
        drawnow
        cb_min_slide(curr_rad).Enable = 'on';
    end

    function slide_dist_min1(src, event)
        [curr_rad, curr_rad_alt] ...
                            = deal(1, 2);
        slide_dist_min
    end

    function slide_dist_min2(src, event)
        [curr_rad, curr_rad_alt] ...
                            = deal(2, 1);
        slide_dist_min
    end

    function slide_dist_min(src, event)
        if ((1e3 * dist_min_slide(curr_rad).Value) < dist_max(curr_rad))
            if distfix_check(curr_rad)
                tmp1        = dist_max(curr_rad) - dist_min(curr_rad);
            end
            dist_min(curr_rad) = 1e3 * dist_min_slide(curr_rad).Value;
            if distfix_check(curr_rad)
                dist_max(curr_rad) = dist_min(curr_rad) + tmp1;
                if (dist_max(curr_rad) > dist_max_ref(curr_rad))
                    dist_max(curr_rad) = dist_max_ref(curr_rad);
                    dist_min(curr_rad) = dist_max(curr_rad) - tmp1;
                    if (dist_min(curr_rad) < dist_min_ref(curr_rad))
                        dist_min(curr_rad) = dist_min_ref(curr_rad);
                    end
                    if (dist_min(curr_rad) < (1e3 * dist_min_slide(curr_rad).Min))
                        dist_min_slide(curr_rad).Value = dist_min_slide(curr_rad).Min;
                    else
                        dist_min_slide(curr_rad).Value = 1e-3 * dist_min(curr_rad);
                    end
                end
                dist_max_edit(curr_rad).String = sprintf('%3.0f', (1e-3 * dist_max(curr_rad)));
                if (dist_max(curr_rad) > (1e3 * dist_max_slide(curr_rad).Max))
                    dist_max_slide(curr_rad).Value = dist_max_slide(curr_rad).Max;
                else
                    dist_max_slide(curr_rad).Value = 1e-3 * dist_max(curr_rad);
                end
            end
            dist_min_edit(curr_rad).String = sprintf('%3.0f', (1e-3 * dist_min(curr_rad)));
            update_dist_range
        else
            if (dist_min(curr_rad) < (1e3 * dist_min_slide(curr_rad).Min))
                dist_min_slide(curr_rad).Value = dist_min_slide(curr_rad).Min;
            else
                dist_min_slide(curr_rad).Value = 1e-3 * dist_min(curr_rad);
            end
        end
        dist_min_slide(curr_rad).Enable = 'off';
        drawnow
        dist_min_slide(curr_rad).Enable = 'on';
    end

    function slide_z_min(src, event)
        if((depth_max_ref - (z_min_slide.Value - depth_min_ref)) > depth_min)
            if zfix_check
                tmp1        = depth_max - depth_min;
            end
            depth_max = depth_max_ref - (z_min_slide.Value - depth_min_ref);
            if zfix_check
                depth_min = depth_max - tmp1;
                if (depth_min < depth_min_ref)
                    depth_min = depth_min_ref;
                    depth_max = depth_min + tmp1;
                    if (depth_max > depth_max_ref)
                        depth_max = depth_max_ref;
                    end
                    if ((depth_max_ref - (depth_max - depth_min_ref)) < z_min_slide.Min)
                        z_min_slide.Value = z_min_slide.Min;
                    elseif ((depth_max_ref - (depth_max - depth_min_ref)) > z_min_slide.Max)
                        z_min_slide.Value = z_min_slide.Max;
                    else
                        z_min_slide.Value = (depth_max_ref - (depth_max - depth_min_ref));
                    end
                end
                z_min_edit.String = sprintf('%4.0f', depth_max);
                if ((depth_max_ref - (depth_min - depth_min_ref)) < z_max_slide.Min)
                    z_max_slide.Value = z_max_slide.Min;
                elseif ((depth_max_ref - (depth_min - depth_min_ref)) > z_max_slide.Max)
                    z_max_slide.Value = z_max_slide.Max;
                else
                    z_max_slide.Value = (depth_max_ref - (depth_min - depth_min_ref));
                end
            end
            z_min_edit.String = sprintf('%4.0f', depth_max);
            update_z_range
        else
            if ((depth_max_ref - (depth_max - depth_min_ref)) < z_min_slide.Min)
                z_min_slide.Value = z_min_slide.Min;
            elseif ((depth_max_ref - (depth_max - depth_min_ref)) > z_min_slide.Max)
                z_min_slide.Value = z_min_slide.Max;
            else
                z_min_slide.Value = (depth_max_ref - (depth_max - depth_min_ref));
            end
        end
        z_min_slide.Enable = 'off';
        drawnow
        z_min_slide.Enable = 'on';
    end

%% Update maximum dB/z/dist

    function slide_db_max1(src, event)
        [curr_rad, curr_rad_alt] ...
                            = deal(1, 2);
        slide_db_max
    end

    function slide_db_max2(src, event)
        [curr_rad, curr_rad_alt] ...
                            = deal(2, 1);
        slide_db_max
    end

    function slide_db_max(src, event)
        if (cb_max_slide(curr_rad).Value > db_min(curr_rad))
            db_max(curr_rad) = cb_max_slide(curr_rad).Value;
            cb_max_edit(curr_rad).String = sprintf('%3.0f', db_max(curr_rad));
            update_db_range
        else
            if (db_max(curr_rad) > cb_max_slide(curr_rad).Max)
                cb_max_slide(curr_rad).Value = cb_max_slide(curr_rad).Max;
            else
                cb_max_slide(curr_rad).Value = db_max(curr_rad);
            end
        end
        cb_max_slide(curr_rad).Enable = 'off';
        drawnow
        cb_max_slide(curr_rad).Enable = 'on';
    end

    function slide_dist_max1(src, event)
        [curr_rad, curr_rad_alt] ...
                            = deal(1, 2);
        slide_dist_max
    end

    function slide_dist_max2(src, event)
        [curr_rad, curr_rad_alt] ...
                            = deal(2, 1);
        slide_dist_max
    end

    function slide_dist_max(src, event)
        if ((1e3 * dist_max_slide(curr_rad).Value) > dist_min(curr_rad))
            if distfix_check(curr_rad)
                tmp1        = dist_max(curr_rad) - dist_min(curr_rad);
            end
            dist_max(curr_rad) = 1e3 * dist_max_slide(curr_rad).Value;
            if distfix_check(curr_rad)
                dist_min(curr_rad) = dist_max(curr_rad) - tmp1;
                if (dist_min(curr_rad) < dist_min_ref(curr_rad))
                    dist_min(curr_rad) = dist_min_ref(curr_rad);
                    dist_max(curr_rad) = dist_min(curr_rad) + tmp1;
                    if (dist_max(curr_rad) > dist_max_ref(curr_rad))
                        dist_max(curr_rad) = dist_max_ref(curr_rad);
                    end
                    if (dist_max(curr_rad) > (1e3 * dist_max_slide(curr_rad).Max))
                        dist_max_slide(curr_rad).Value = dist_max_slide(curr_rad).Max;
                    else
                        dist_max_slide(curr_rad).Value = 1e-3 * dist_max(curr_rad);
                    end
                end
                dist_min_edit(curr_rad).String = sprintf('%3.0f', (1e-3 * dist_min(curr_rad)));
                if (dist_min(curr_rad) < dist_min_slide(curr_rad).Min)
                    dist_min_slide(curr_rad).Value = dist_min_slide(curr_rad).Min;
                else
                    dist_min_slide(curr_rad).Value = 1e-3 * dist_min(curr_rad);
                end
            end
            dist_max_edit(curr_rad).String = sprintf('%3.0f', (1e-3 * dist_max(curr_rad)));
            update_dist_range
        else
            if (dist_max(curr_rad) > (1e3 * dist_max_slide(curr_rad).Max))
                dist_max_slide(curr_rad).Value = dist_max_slide(curr_rad).Max;
            else
                dist_max_slide(curr_rad).Value = 1e-3 * dist_max(curr_rad);
            end
        end
        dist_max_slide(curr_rad).Enable = 'off';
        drawnow
        dist_max_slide(curr_rad).Enable = 'on';
    end

    function slide_z_max(src, event)
        if ((depth_max_ref - (z_max_slide.Value - depth_min_ref)) < depth_max)
            if zfix_check
                tmp1        = depth_max - depth_min;
            end
            depth_min = depth_max_ref - (z_max_slide.Value - depth_min_ref);
            if zfix_check
                depth_max = depth_min + tmp1;
                if (depth_max > depth_max_ref)
                    depth_max = depth_max_ref;
                    depth_min = depth_max - tmp1;
                    if (depth_min < depth_min_ref)
                        depth_min = depth_min_ref;
                    end
                    if ((depth_max_ref - (depth_min - depth_min_ref)) < z_max_slide.Min)
                        z_max_slide.Value = z_max_slide.Min;
                    elseif ((depth_max_ref - (depth_min - depth_min_ref)) > z_max_slide.Max)
                        z_max_slide.Value = z_max_slide.Max;
                    else
                        z_max_slide.Value = (depth_max_ref - (depth_min - depth_min_ref));
                    end
                end
                z_min_edit.String = sprintf('%4.0f', depth_max);
                if ((depth_max_ref - (depth_max - depth_min_ref)) < z_min_slide.Min)
                    z_min_slide.Value = z_min_slide.Min;
                elseif ((depth_max_ref - (depth_max - depth_min_ref)) > z_min_slide.Max)
                    z_min_slide.Value = z_min_slide.Max;
                else
                    z_min_slide.Value = (depth_max_ref - (depth_max - depth_min_ref));
                end
            end
            z_max_edit.String = sprintf('%4.0f', depth_min);
            update_z_range
        else
            if ((depth_max_ref - (depth_min - depth_min_ref)) < z_max_slide.Min)
                z_max_slide.Value = z_max_slide.Min;
            elseif ((depth_max_ref - (depth_min - depth_min_ref)) > z_max_slide.Max)
                z_max_slide.Value = z_max_slide.Max;
            else
                z_max_slide.Value = (depth_max_ref - (depth_min - depth_min_ref));
            end
        end
        z_max_slide.Enable = 'off';
        drawnow
        z_max_slide.Enable = 'on';
    end

%% Reset minimum dB/x/y/z/dist

	function reset_db_min1(src, event)
        [curr_rad, curr_rad_alt] ...
                            = deal(1, 2);
        rad_group.SelectedObject = rad_check(curr_rad);
        reset_db_min
    end

	function reset_db_min2(src, event)
        [curr_rad, curr_rad_alt] ...
                            = deal(2, 1);
        rad_group.SelectedObject = rad_check(curr_rad);
        reset_db_min
    end

    function reset_db_min(src, event)
        if (db_min_ref(curr_rad) < cb_min_slide(curr_rad).Min)
            cb_min_slide(curr_rad).Value = cb_min_slide(curr_rad).Min;
        else
            cb_min_slide(curr_rad).Value = db_min_ref(curr_rad);
        end
        cb_min_edit(curr_rad).String = num2str(db_min_ref(curr_rad));
        db_min(curr_rad)     = db_min_ref(curr_rad);
        update_db_range
    end

    function reset_dist_min(src, event)
        if (dist_min_ref(curr_rad) < (1e3 * dist_min_slide(curr_rad).Min))
            dist_min_slide(curr_rad).Value = dist_min_slide(curr_rad).Min;
        else
            dist_min_slide(curr_rad).Value = 1e-3 * dist_min_ref(curr_rad);
        end
        dist_min_edit(curr_rad).String = sprintf('%3.1f', (1e-3 * dist_min_ref(curr_rad)));
        dist_min(curr_rad)  = dist_min_ref(curr_rad);
        update_dist_range
    end

    function reset_z_min(src, event)
    	depth_max			= depth_max_ref;
    	if (depth_min_ref < z_min_slide.Min)
        	z_min_slide.Value = z_min_slide.Min;
    	elseif (depth_min_ref > z_min_slide.Max)
        	z_min_slide.Value = z_min_slide.Max;
    	else
        	z_min_slide.Value = depth_min_ref;
    	end
    	z_min_edit.String = sprintf('%4.0f', depth_max_ref);
        update_z_range
    end

%% Reset maximum dB/x/y/z

    function reset_db_max1(src, event)
        [curr_rad, curr_rad_alt] ...
                            = deal(1, 2);
        reset_db_max
    end

    function reset_db_max2(src, event)
        [curr_rad, curr_rad_alt] ...
                            = deal(2, 1);
        reset_db_max
    end

    function reset_db_max(src, event)
        if (db_max_ref(curr_rad) > cb_max_slide(curr_rad).Max)
            cb_max_slide(curr_rad).Value = cb_max_slide(curr_rad).Max;
        else
            cb_max_slide(curr_rad).Value = db_max_ref(curr_rad);
        end
        cb_max_edit(curr_rad).String = num2str(db_max_ref(curr_rad));
        db_max(curr_rad)     = db_max_ref(curr_rad);
        update_db_range
    end

    function reset_dist_max(src, event)
        if (dist_max_ref(curr_rad) > (1e3 * dist_max_slide(curr_rad).Max))
            dist_max_slide(curr_rad).Value = dist_max_slide(curr_rad).Max;
        else
            dist_max_slide(curr_rad).Value = 1e-3 * dist_max_ref(curr_rad);
        end
        dist_max_edit(curr_rad).String = sprintf('%3.1f', (1e-3 * dist_max_ref(curr_rad)));
        dist_max(curr_rad)  = dist_max_ref(curr_rad);
        update_dist_range
    end

    function reset_z_max(src, event)
        depth_min			= depth_min_ref;
        if (depth_min_ref < z_min_slide.Min)
            z_min_slide.Value = z_min_slide.Min;
        elseif (depth_min_ref > z_min_slide.Max)
            z_min_slide.Value = z_min_slide.Max;
        else
            z_min_slide.Value = depth_min_ref;
        end
        z_min_edit.String = sprintf('%4.0f', depth_max_ref);
        update_z_range
    end

%% Reset all x/z

    function reset_xz1(src, event)
        [curr_rad, curr_rad_alt] ...
                            = deal(1, 2);
        reset_xz
    end

    function reset_xz2(src, event)
        [curr_rad, curr_rad_alt] ...
                            = deal(2, 1);
        reset_xz
    end

    function reset_xz(src, event)
        reset_dist_min
        reset_dist_max
        reset_z_min
        reset_z_max
    end

%% Update dB/x/y/z/depth range

    function update_db_range(src, event)
        rad_group.SelectedObject = rad_check(curr_rad);        
        ax(curr_rad).CLim = [db_min(curr_rad) db_max(curr_rad)];
    end

    function update_dist_range(src, event)
		rad_group.SelectedObject = rad_check(curr_rad);
        ax(curr_rad).XLim = 1e-3 .* [dist_min(curr_rad) dist_max(curr_rad)];
        narrow_cb
    end

    function update_z_range(src, event)
        rad_group.SelectedObject = rad_check(curr_rad);
		ax(curr_rad).YLim = [depth_min depth_max];
        narrow_cb
        switch curr_rad
            case 1
                curr_rad	= 2;
                narrow_cb
                curr_rad	= 1;
            case 2
                curr_rad	= 1;
                narrow_cb
                curr_rad	= 2;
        end
    end

%% Adjust slider limits after panning or zooming

    function panzoom(src, event)
        rad_group.SelectedObject = rad_check(curr_rad);
        tmp1                = 1e3 .* ax(curr_rad).XLim;
        if (tmp1(1) < dist_min_ref(curr_rad))
            reset_dist_min
        else
            if (tmp1(1) < (1e3 * dist_min_slide(curr_rad).Min))
                dist_min_slide(curr_rad).Value = dist_min_slide(curr_rad).Min;
            else
                dist_min_slide(curr_rad).Value = 1e-3 * tmp1(1);
            end
            dist_min_edit(curr_rad).String = sprintf('%3.1f', (1e-3 * tmp1(1)));
            dist_min(curr_rad) ...
                            = tmp1(1);
        end
        if (tmp1(2) > dist_max_ref(curr_rad))
            reset_dist_max
        else
            if (tmp1(2) > (1e3 * dist_max_slide(curr_rad).Max))
                dist_max_slide(curr_rad).Value = dist_max_slide(curr_rad).Max;
            else
                dist_max_slide(curr_rad).Value = 1e-3 * tmp1(2);
            end
            dist_max_edit(curr_rad).String = sprintf('%3.1f', (1e-3 * tmp1(2)));
            dist_max(curr_rad) ...
                            = tmp1(2);
        end
        tmp1                = ax(curr_rad).YLim;
        tmp2				= [depth_min depth_max];
        if (tmp1(1) < depth_min_ref)
            reset_z_min
        else
            if (tmp1(1) < z_min_slide.Min)
                z_min_slide.Value = z_min_slide.Min;
            elseif (tmp1(1) > z_min_slide.Max)
                z_min_slide.Value = z_min_slide.Max;
            else
                z_min_slide.Value = tmp1(1);
            end
            z_min_edit.String = sprintf('%4.0f', tmp1(1));
            depth_min		= tmp1(1);
        end
		if (tmp1(2) > depth_max_ref)
            reset_z_max
        else
            if (tmp1(2) < z_max_slide.Min)
                z_max_slide.Value = z_max_slide.Min;
            elseif (tmp1(2) > z_max_slide.Max)
                z_max_slide.Value = z_max_slide.Max;
            else
                z_max_slide.Value = tmp1(2);
            end
            z_max_edit.String = sprintf('%4.0f', tmp1(2));
            depth_max		= tmp1(2);
		end
        narrow_cb
		switch curr_rad
            case 1
                curr_rad	= 2;
                narrow_cb
                curr_rad	= 1;
            case 2
                curr_rad	= 1;
                narrow_cb
                curr_rad	= 2;
		end
    end

%% Plot data (depth)

    function plot_data(src, event)
        if ~data_done(curr_rad)
            status_box.String = 'Data not loaded yet.';
            return
        end
		z_min_slide.Min = depth_min_ref; z_min_slide.Max = depth_max_ref; z_min_slide.Value = (depth_max_ref - (depth_max - depth_min_ref));
		z_max_slide.Min = depth_min_ref; z_max_slide.Max = depth_max_ref; z_max_slide.Value = (depth_max_ref - (depth_min - depth_min_ref));
		z_min_edit.String = sprintf('%4.0f', depth_max);
		z_max_edit.String = sprintf('%4.0f', depth_min);
		if isgraphics(p_data(curr_rad))
			[p_data(curr_rad).XData, p_data(curr_rad).YData, p_data(curr_rad).CData] ...
						= deal((1e-3 .* dist_lin{curr_rad}), depth_ref{curr_rad}, amp_depth{curr_rad});
		else
			p_data(curr_rad) = imagesc(ax(curr_rad), (1e-3 .* dist_lin{curr_rad}), depth_ref{curr_rad}, amp_depth{curr_rad}, [db_min(curr_rad) db_max(curr_rad)]);
		end
		ax(curr_rad).YLim = [depth_min depth_max];
        show_data
        show_pk
    end

%% Show radar data

    function show_data(src, event)
        rad_group.SelectedObject = rad_check(curr_rad);
        if data_done(curr_rad)
            if (data_check(curr_rad) && isgraphics(p_data(curr_rad)))
                p_data(curr_rad).Visible = 'on';
            elseif isgraphics(p_data(curr_rad))
                p_data(curr_rad).Visible = 'off';
            end
        elseif data_check(curr_rad)
            data_check(curr_rad) = false;
        end
    end

%% Show picked layers

    function show_pk(src, event)
        rad_group.SelectedObject = rad_check(curr_rad);
        if pk_done(curr_rad)
            if pk_check(curr_rad)
                [p_pk{curr_rad}(isgraphics(p_pk{curr_rad})).Visible, p_bed(curr_rad).Visible] = deal('on');
                uistack([p_pk{curr_rad}(isgraphics(p_pk{curr_rad})) p_bed(curr_rad)], 'top')
				show_int
            else
                [p_pk{curr_rad}(isgraphics(p_pk{curr_rad})).Visible, p_bed(curr_rad).Visible] = deal('off');
            end
        elseif pk_check(curr_rad)
            pk_check(curr_rad) = true;
        end
    end

%% Show intersections

    function show_int1(src, event)
        [curr_rad, curr_rad_alt] ...
                            = deal(1, 2);
        show_int
    end

    function show_int2(src, event)
        [curr_rad, curr_rad_alt] ...
                            = deal(2, 1);
        show_int
    end

    function show_int(src, event)
        rad_group.SelectedObject = rad_check(curr_rad);
        if all(pk_done)
            if int_check(curr_rad).Value
                [p_int1{curr_rad}(isgraphics(p_int1{curr_rad})).Visible, p_int2{curr_rad}(isgraphics(p_int2{curr_rad})).Visible] = deal('on');
                uistack(p_int1{curr_rad}(isgraphics(p_int1{curr_rad})), 'top')
                uistack(p_int2{curr_rad}(isgraphics(p_int2{curr_rad})), 'top')
            else
				[p_int1{curr_rad}(isgraphics(p_int1{curr_rad})).Visible, p_int2{curr_rad}(isgraphics(p_int2{curr_rad})).Visible] = deal('off');
            end
        elseif int_check(curr_rad).Value
            int_check(curr_rad).Value = 0;
        end
	end

%% Change intersection

    function change_int(src, event)
		if ~(num_int && (curr_ind_int(1) ~= 0))
            return
		end
        curr_int            = intnum_list.Value;
        for ii = 1:2
            [p_int1{ii}(isgraphics(p_int1{ii})).LineWidth] = deal(2);
            if any(isgraphics(p_int1{ii}(curr_int)))
                p_int1{ii}(curr_int).LineWidth = 4;
            end
        end
        [dist_min(1), dist_max(1), dist_min(2), dist_max(2)] ...
                            = deal((dist_lin{1}(curr_ind_int(curr_int, 1)) - dist_int_half), (dist_lin{1}(curr_ind_int(curr_int, 1)) + dist_int_half), ...
								   (dist_lin{2}(curr_ind_int(curr_int, 2)) - dist_int_half), (dist_lin{2}(curr_ind_int(curr_int, 2)) + dist_int_half));
        if (dist_min(1) < dist_min_ref(1))
            dist_min(1)		= dist_min_ref(1);
        end
        if (dist_min(2) < dist_min_ref(2))
            dist_min(2)		= dist_min_ref(2);
        end
        if (dist_max(1) > dist_max_ref(1))
            dist_max(1)		= dist_max_ref(1);
        end
		if (dist_max(2) > dist_max_ref(2))
            dist_max(2)		= dist_max_ref(2);
		end
        ax(1).XLim = 1e-3 .* [dist_min(1) dist_max(1)];
        ax(2).XLim = 1e-3 .* [dist_min(2) dist_max(2)];
        dist_min_slide(1).Value = 1e-3 * dist_min(1);
        dist_min_slide(2).Value = 1e-3 * dist_min(2);
        dist_max_slide(1).Value = 1e-3 * dist_max(1);
        dist_max_slide(2).Value = 1e-3 * dist_max(2);
        dist_min_edit(1).String = sprintf('%3.0f', (1e-3 * dist_min(1)));
        dist_min_edit(2).String = sprintf('%3.0f', (1e-3 * dist_min(2)));
        dist_max_edit(1).String = sprintf('%3.0f', (1e-3 * dist_max(1)));
        dist_max_edit(2).String = sprintf('%3.0f', (1e-3 * dist_max(2)));
    	depth_min			= 0;
    	z_max_slide.Value = (depth_max_ref - (depth_min - depth_min_ref));
    	z_max_edit.String = sprintf('%4.0f', depth_min);
    	if any(~isnan([depth_bed{1}(curr_ind_int(curr_int, 1)) depth_bed{2}(curr_ind_int(curr_int, 2))]))
        	depth_max		= mean([depth_bed{1}(curr_ind_int(curr_int, 1)) depth_bed{2}(curr_ind_int(curr_int, 2))], 'omitnan');
        	z_min_slide.Value = (depth_max_ref - (depth_max - depth_min_ref));
        	z_min_edit.String = sprintf('%4.0f', depth_max);
    	end
        update_z_range
    end

%% Toggle gridlines
	
    function toggle_grid
        rad_group.SelectedObject = rad_check(curr_rad);
        if grid_check(curr_rad)
            [ax(:).YGrid] = deal('on');
        else
            [ax(:).YGrid] = deal('off');
        end
    end

%% Narrow color axis to +/- 2 standard deviations of current mean value

    function narrow_cb(src, event)
		if isequal(narrow_ax{curr_rad}, [ax(curr_rad).XLim ax(curr_rad).YLim])
			return
		end
        rad_group.SelectedObject = rad_check(curr_rad);
        if ~data_done(curr_rad)
            return
        end
        tmp1                = zeros(2);
        tmp1(1, :)          = interp1(dist_lin{curr_rad}, 1:num_trace(curr_rad), [dist_min(curr_rad) dist_max(curr_rad)], 'nearest', 'extrap');
        tmp1(2, :)			= interp1(depth_ref{curr_rad}, 1:num_sample(curr_rad), [depth_min depth_max], 'nearest', 'extrap');
        tmp2                = NaN(1, 2);
        [tmp2(1), tmp2(2)]  = deal(mean(amp_depth{curr_rad}(tmp1(2, 1):10:tmp1(2, 2), tmp1(1, 1):10:tmp1(1, 2)), 'all', 'omitnan'), ...
								   std(amp_depth{curr_rad}(tmp1(2, 1):10:tmp1(2, 2), tmp1(1, 1):10:tmp1(1, 2)), 0, 'all', 'omitnan'));
        if any(isnan(tmp2))
            return
        end
        tmp1                = zeros(1, 2);
        tmp4                = 2 * tmp2(2);
        if ((tmp2(1) - tmp4) < db_min_ref(curr_rad))
            tmp1(1)         = db_min_ref(curr_rad);
        else
            tmp1(1)         = tmp2(1) - tmp4;
        end
        if ((tmp2(1) + tmp4) > db_max_ref(curr_rad))
            tmp1(2)         = db_max_ref(curr_rad);
        else
            tmp1(2)         = tmp2(1) + tmp4;
        end
        [db_min(curr_rad), db_max(curr_rad)] ...
                            = deal(tmp1(1), tmp1(2));
        if (db_min(curr_rad) < cb_min_slide(curr_rad).Min)
            cb_min_slide(curr_rad).Value = cb_min_slide(curr_rad).Min;
        else
            cb_min_slide(curr_rad).Value = db_min(curr_rad);
        end
        if (db_max(curr_rad) > cb_max_slide(curr_rad).Max)
            cb_max_slide(curr_rad).Value = cb_max_slide(curr_rad).Max;
        else
            cb_max_slide(curr_rad).Value = db_max(curr_rad);
        end
        cb_min_edit(curr_rad).String = sprintf('%3.0f', db_min(curr_rad));
        cb_max_edit(curr_rad).String = sprintf('%3.0f', db_max(curr_rad));
        ax(curr_rad).CLim = [db_min(curr_rad) db_max(curr_rad)];
		narrow_ax{curr_rad}	= [ax(curr_rad).XLim ax(curr_rad).YLim];
    end

%% Switch active segment

    function rad_radio(~, event)
        switch event.NewValue.String
            case 'P'
                [curr_rad, curr_rad_alt] ...
                            = deal(1, 2);
                status_box.String = 'Primary segment now selected/active.';
            case 'I'
                [curr_rad, curr_rad_alt] ...
                            = deal(2, 1);
                status_box.String = 'Intersecting segment now selected/active.';
        end
    end

%% Keyboard shortcuts

    function keypress(~, event)
        switch event.Key
            case '1'
				if any(pk_check)
					pk_check			= false(1, 2);
				else
					pk_check			= true(1, 2);
				end
				show_pk
				[curr_rad, curr_rad_alt] ...
										= deal(curr_rad_alt, curr_rad);
				show_pk
				[curr_rad, curr_rad_alt] ...
										= deal(curr_rad_alt, curr_rad);
            case '2'
                rad_group.SelectedObject = rad_check(1);
                [curr_rad, curr_rad_alt] ...
                            = deal(1, 2);
				if pk_check(curr_rad)
					pk_check(curr_rad) = false;
				else
					pk_check(curr_rad) = true;
				end
				show_pk
            case '3'
                rad_group.SelectedObject = rad_check(2);
                [curr_rad, curr_rad_alt] ...
                            = deal(2, 1);
				if pk_check(curr_rad)
					pk_check(curr_rad) = false;
				else
					pk_check(curr_rad) = true;
				end
				show_pk
            case '4'
				rad_group.SelectedObject = rad_check(1);
                [curr_rad, curr_rad_alt] ...
                            = deal(1, 2);
                if data_check(curr_rad)
                    data_check(curr_rad) = false;
                else
                    data_check(curr_rad) = true;
                end
                show_data
            case '5'
				rad_group.SelectedObject = rad_check(2);
                [curr_rad, curr_rad_alt] ...
                            = deal(2, 1);
                if data_check(curr_rad)
                    data_check(curr_rad) = false;
                else
                    data_check(curr_rad) = true;
                end
                show_data
            case '6'
                [curr_rad, curr_rad_alt] ...
                            = deal(1, 2);
                if int_check(curr_rad).Value
                    int_check(curr_rad).Value = 0;
                else
                    int_check(curr_rad).Value = 1;
                end
                show_int
            case '7'
                [curr_rad, curr_rad_alt] ...
                            = deal(2, 1);
                if int_check(curr_rad).Value
                    int_check(curr_rad).Value = 0;
                else
                    int_check(curr_rad).Value = 1;
                end
                show_int
            case 'a'
                pk_last
            case 'e'
                reset_xz
            case 'g'
                if grid_check(curr_rad)
                    grid_check(curr_rad) = false;
                else
                    grid_check(curr_rad) = true;
                end
                toggle_grid
            case 'h'
                if match_check.Value
                    match_check.Value = 0;
                else
                    match_check.Value = 1;
                end
            case 'm'
                pk_match
            case 'n'
                pk_next
            case 'o'
                pk_focus
			case 'p'
				if (curr_int < num_int)
					status_box.String = 'Review all intersections first.';
				else
					if (int_list.Value < length(int_list.String))
						if match_done
							pk_save
						else
							int_list.Value = int_list.Value + 1;
							int_count_edit.String = [num2str(int_list.Value) ' /'];
							load_pk2
							load_data
						end
					else
						pk_save
						status_box.String = 'End of intersecting list.';
					end
				end
            case 't'
                if (curr_int > 1)
                    intnum_list.Value = (curr_int - 1);
                    change_int
                end
            case 'u'
                pk_unmatch
            case 'v'
				if nearest_check.Value
                    nearest_check.Value = 0;
                else
                    nearest_check.Value = 1;
				end
            case 'y'
                if (curr_int < num_int)
                    intnum_list.Value = (curr_int + 1);
                    change_int
                end
			case 'z'
				if zfix_check
                    zfix_check = false;
					status_box.String = 'Depth scale free.';
                else
                    zfix_check = true;
					status_box.String = 'Depth scale fixed.';
				end				
            case 'downarrow'
                tmp1 = depth_max - depth_min;
                tmp2 = depth_max + (0.25 * tmp1);
                if (tmp2 > depth_max_ref)
                    depth_max = depth_max_ref;
                else
                    depth_max = tmp2;
                end
                depth_min = depth_max - tmp1;
                if (depth_min < depth_min_ref)
                    depth_min = depth_min_ref;
                end
                z_min_edit.String = sprintf('%4.0f', depth_max);
                z_max_edit.String = sprintf('%4.0f', depth_min);
                if ((depth_max_ref - (depth_max - depth_min_ref)) < z_min_slide.Min)
                    z_min_slide.Value = z_min_slide.Min;
                elseif ((depth_max_ref - (depth_max - depth_min_ref)) > z_min_slide.Max)
                    z_min_slide.Value = z_min_slide.Max;
                else
                    z_min_slide.Value = (depth_max_ref - (depth_max - depth_min_ref));
                end
                if ((depth_max_ref - (depth_min - depth_min_ref)) < z_max_slide.Min)
                    z_max_slide.Value = z_max_slide.Min;
                elseif ((depth_max_ref - (depth_min - depth_min_ref)) > z_max_slide.Max)
                    z_max_slide.Value = z_max_slide.Max;
                else
                    z_max_slide.Value = (depth_max_ref - (depth_min - depth_min_ref));
                end
                update_z_range
            case 'leftarrow'
                tmp1        = dist_max(curr_rad) - dist_min(curr_rad);
                tmp2        = dist_min(curr_rad) - (0.25 * tmp1);
                if (tmp2 < dist_min_ref(curr_rad))
                    dist_min(curr_rad) = dist_min_ref(curr_rad);
                else
                    dist_min(curr_rad) = tmp2;
                end
                dist_max(curr_rad) = dist_min(curr_rad) + tmp1;
                if (dist_max(curr_rad) > dist_max_ref(curr_rad))
                    dist_max(curr_rad) = dist_max_ref(curr_rad);
                end
                dist_min_edit(curr_rad).String = sprintf('%3.1f', (1e-3 * dist_min(curr_rad)));
                dist_max_edit(curr_rad).String = sprintf('%3.1f', (1e-3 * dist_max(curr_rad)));
                if (dist_min(curr_rad) < (1e3 * dist_min_slide(curr_rad).Min))
                    dist_min_slide(curr_rad).Value = dist_min_slide(curr_rad).Min;
                else
                    dist_min_slide(curr_rad).Value = 1e-3 * dist_min(curr_rad);
                end
                if (dist_max(curr_rad) > (1e3 * dist_max_slide(curr_rad).Max))
                    dist_max_slide(curr_rad).Value = dist_max_slide(curr_rad).Max;
                else
                    dist_max_slide(curr_rad).Value = 1e-3 * dist_max(curr_rad);
                end
                update_dist_range
            case 'rightarrow'
                tmp1        = dist_max(curr_rad) - dist_min(curr_rad);
                tmp2        = dist_max(curr_rad) + (0.25 * tmp1);
                if (tmp2 > dist_max_ref(curr_rad))
                    dist_max(curr_rad) ...
                            = dist_max_ref(curr_rad);
                else
                    dist_max(curr_rad) ...
                            = tmp2;
                end
                dist_min(curr_rad) ...
                            = dist_max(curr_rad) - tmp1;
                if (dist_min(curr_rad) < dist_min_ref(curr_rad))
                    dist_min(curr_rad) ...
                            = dist_min_ref(curr_rad);
                end
                dist_min_edit(curr_rad).String = sprintf('%3.1f', (1e-3 * dist_min(curr_rad)));
                dist_max_edit(curr_rad).String = sprintf('%3.1f', (1e-3 * dist_max(curr_rad)));
                if (dist_min(curr_rad) < (1e3 * dist_min_slide(curr_rad).Min))
                    dist_min_slide(curr_rad).Value = dist_min_slide(curr_rad).Min;
                else
                    dist_min_slide(curr_rad).Value = 1e-3 * dist_min(curr_rad);
                end
                if (dist_max(curr_rad) > (1e3 * dist_max_slide(curr_rad).Max))
                    dist_max_slide(curr_rad).Value = dist_max_slide(curr_rad).Max;
                else
                    dist_max_slide(curr_rad).Value = 1e-3 * dist_max(curr_rad);
                end
                update_dist_range
            case 'uparrow'
                tmp1 = depth_max - depth_min;
                tmp2 = depth_min - (0.25 * tmp1);
                if (tmp2 < depth_min_ref)
                    depth_min = depth_min_ref;
                else
                    depth_min = tmp2;
                end
                depth_max = depth_min + tmp1;
                depth_min = depth_max - tmp1;
                z_min_edit.String = sprintf('%4.0f', depth_max);
                z_max_edit.String = sprintf('%4.0f', depth_min);
                if ((depth_max_ref - (depth_max - depth_min_ref)) < z_min_slide.Min)
                    z_min_slide.Value = z_min_slide.Min;
                elseif ((depth_max_ref - (depth_max - depth_min_ref)) > z_min_slide.Max)
                    z_min_slide.Value = z_min_slide.Max;
                else
                    z_min_slide.Value = (depth_max_ref - (depth_max - depth_min_ref));
                end
                if ((depth_max_ref - (depth_min - depth_min_ref)) < z_max_slide.Min)
                    z_max_slide.Value = z_max_slide.Min;
                elseif ((depth_max_ref - (depth_min - depth_min_ref)) > z_max_slide.Max)
                    z_max_slide.Value = z_max_slide.Max;
                else
                    z_max_slide.Value = (depth_max_ref - (depth_min - depth_min_ref));
                end
                update_z_range
			case 'quote'
				if distfix_check
                    distfix_check = false;
					status_box.String = 'Distance scale free.';
                else
                    distfix_check = true;
					status_box.String = 'Distance scale fixed.';
				end
        end
    end

%% Mouse click

    function mouse_click(src, event)
		if any(~pk_done)
			return
		end
        tmp1                = src.CurrentPoint; % picked x/y pixels
        tmp2                = fgui.Position; % figure x0/y0/w/h normalized
        tmp3                = [ax(1).Position; ax(2).Position]; % axis 1 x0/y0/w/h; 2 x0/y0/w/h normalized
        tmp4                = [(tmp2(1) + (tmp2(3) * tmp3(1, 1))) (tmp2(1) + (tmp2(3) * (tmp3(1, 1) + tmp3(1, 3)))); (tmp2(2) + (tmp2(4) * tmp3(1, 2))) (tmp2(2) + (tmp2(4) * (tmp3(1, 2) + tmp3(1, 4))))]; % 1 x0/y1;y0/y1 pixels
        tmp5                = [(tmp2(1) + (tmp2(3) * tmp3(2, 1))) (tmp2(1) + (tmp2(3) * (tmp3(2, 1) + tmp3(2, 3)))); (tmp2(2) + (tmp2(4) * tmp3(2, 2))) (tmp2(2) + (tmp2(4) * (tmp3(2, 2) + tmp3(2, 4))))]; % 2 x0/y1;y0/y1 pixels
        if ((tmp1(1) > (tmp4(1, 1))) && (tmp1(1) < (tmp4(1, 2))) && (tmp1(2) > (tmp4(2, 1))) && (tmp1(2) < (tmp4(2, 2))))
            [curr_rad, curr_rad_alt] ...
                            = deal(1, 2);
			tmp1    = [((tmp1(1) - tmp4(1, 1)) / diff(tmp4(1, :))) ((tmp4(2, 2) - tmp1(2)) / diff(tmp4(2, :)))];
            tmp2            = [ax(curr_rad).XLim; ax(curr_rad).YLim];
            [ind_x_pk, ind_y_pk] ...
                            = deal(((tmp1(1) * diff(tmp2(1, :))) + tmp2(1, 1)), ((tmp1(2) * diff(tmp2(2, :))) + tmp2(2, 1)));
            pk_select_gui_breakout
        elseif ((tmp1(1) > (tmp5(1, 1))) && (tmp1(1) < (tmp5(1, 2))) && (tmp1(2) > (tmp5(2, 1))) && (tmp1(2) < (tmp5(2, 2))))
            [curr_rad, curr_rad_alt] ...
                            = deal(2, 1);
			tmp1			= [((tmp1(1) - tmp5(1, 1)) / diff(tmp5(1, :))) ((tmp5(2, 2) - tmp1(2)) / diff(tmp5(2, :)))];
            tmp2            = [ax(curr_rad).XLim; ax(curr_rad).YLim];
            [ind_x_pk, ind_y_pk] ...
                            = deal(((tmp1(1) * diff(tmp2(1, :))) + tmp2(1, 1)), ((tmp1(2) * diff(tmp2(2, :))) + tmp2(2, 1)));
            pk_select_gui_breakout
        end
    end

%% select layer interactively

    function pk_select_gui_breakout(src, event)
		[tmp1, tmp2]		= unique(depth_layer{curr_rad}(:, interp1((1e-3 .* dist_lin{curr_rad}), 1:num_trace(curr_rad), ind_x_pk, 'nearest', 'extrap')));
        if (length(tmp1(~isnan(tmp1))) > 1)
            curr_layer(curr_rad) ...
                            = interp1(tmp1(~isnan(tmp1)), tmp2(~isnan(tmp1)), ind_y_pk, 'nearest', 'extrap');
        elseif isscalar(tmp1(~isnan(tmp1)))
            curr_layer(curr_rad) ...
                            = tmp2(find(~isnan(tmp1), 1));
        else
            status_box.String = 'Layer choice unclear.';
            return
        end
        layer_list(curr_rad).Value = curr_layer(curr_rad);
        pk_select
        status_box.String = ['Layer #' num2str(curr_layer(curr_rad)) ' selected.'];
    end

%% Test something

    function misctest(src, event)
		
		fgui.KeyPressFcn = @keypress; fgui.WindowButtonDownFcn = @mouse_click;
        status_box.String = 'Test done.';
    end
%%
end