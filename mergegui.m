function mergegui
% MERGEGUI Interactive merging of layer picks from consecutive, overlapping radar data blocks.
%   
%   MERGEGUI loads a GUI for merging a transect's radar layers traced
%   previously for individual blocks using PICKGUI. Refer to manual for
%   operation (pickgui_man.docx).
%   
%   MERGEGUI requires that the functions TOPOCORR AND SMOOTH_LOWESS be
%   available within the user's path. If the Parallel Computing Toolbox is
%   licensed and available, then several calculations related to data
%   flattening will be parallelized. The Mapping Toolbox is necessary to
%   plot a map of the transect location.
% 
% Joe MacGregor (UTIG)
% Last updated: 03/18/13

if ~exist('topocorr', 'file')
    error('mergegui:topocorr', 'Function TOPOCORR is not available within this user''s path.')
end
if ~exist('smooth_lowess', 'file')
    error('mergegui:smoothlowess', 'Function SMOOTH_LOWESS is not available within this user''s path.')
end

%% Intialize variables

pk                          = struct;

% elevation/depth/distance defaults
[elev_min_ref, elev_max_ref]= deal(0, 1);
[elev_min, elev_max]        = deal(elev_min_ref, elev_max_ref);
[depth_min, depth_min_ref]  = deal(elev_min_ref);
[depth_max, depth_max_ref]  = deal(elev_max_ref);
[dist_min_ref, dist_max_ref]= deal(0, 1);
[dist_min, dist_max]        = deal(dist_min_ref, dist_max_ref);

% dB default
[db_min_ref, db_max_ref]    = deal(-130, 0);
[db_min, db_max]            = deal(-80, -20);

% more default values
speed_vacuum                = 299792458; % m/s
permitt_ice                 = 3.15;
speed_ice                   = speed_vacuum / sqrt(permitt_ice);
decim                       = 10; % decimate radargram for display
length_chunk                = 100; % chunk length in km
disp_type                   = 'amplitude';
cmaps                       = {'bone' 'jet'}';
var_cell                    = {'ind_keep_aresp' 'ind_keep_phase' 'ind_y_start_aresp' 'ind_y_start_phase' 'ind_match' 'ind_x_mean' 'poly_flat' 'file_in' 'ind_y_man' 'file_pk' 'file_block' 'ind_y_flat_mean'}; ...
                              % cell variables that are easy to merge
var_zero                    = {'freq' 'length_smooth' 'ind_trim_start' 'ind_x_start_aresp' 'ind_x_start_phase' 'ind_y_aresp_max' 'ind_y_phase_max'  'num_aresp' 'num_keep_aresp' ...
                               'num_keep_phase' 'num_ind_mean' 'num_man' 'num_phase' 'num_sample' 'num_trace' 'num_win' 'twtt_min_ref' 'twtt_match' 'twtt_max_ref' 'ind_trace_start' 'num_layer_block'};
num_cell                    = length(var_cell);
num_zero                    = length(var_zero);
var_layer                   = {'ind_y' 'twtt' 'twtt_ice' 'int' 'ind_y_smooth' 'ind_y_flat_smooth' 'twtt_smooth' 'twtt_ice_smooth' 'int_smooth' 'depth' 'depth_smooth' 'elev' 'elev_smooth'}; % layer variables
var_pos                     = {'dist' 'lat' 'lon' 'x' 'y' 'elev_air' 'twtt_surf' 'elev_surf' 'twtt_bed' 'elev_bed'}; % position variables
num_var_layer               = length(var_layer);
num_var_pos                 = length(var_pos);
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

% allocate a bunch of variables
[bed_avail, core_done, data_done, edit_flag, flat_done, merge_done, merge_file, surf_avail] ...
                            = deal(false);
[amp_flat, amp_mean, button, colors, curr_chunk, curr_layer, curr_trans, curr_year, depth, depth_bed, depth_bed_flat, depth_curr, depth_flat, depth_layer_flat, depth_mat, dist_chunk, dt, elev, ii, ind_corr, ...
    ind_decim, ind_int, ind_layer_flat, ind_x_pk, ind_x_ref, ind_y_pk, int_core, jj, kk, name_core, name_trans, num_chunk, num_data, num_decim, num_int, num_pk, num_sample, num_trans, num_year, p_bed, p_bedflat, ...
    p_block, p_blockflat, p_blocknum, p_blocknumflat, p_core, p_coreflat, p_corename, p_corenameflat, p_data, pk_all, p_pk, p_pkflat, p_surf, pkfig, p_refflat, rad_threshold, tmp1, tmp2, tmp3, tmp4, tmp5, twtt] ...
                            = deal(0);
[file_data, file_core, file_pk, file_pk_short, file_save, path_core, path_data, path_pk, path_save] ...
                            = deal('');
layer_str                   = {};

if license('checkout', 'distrib_computing_toolbox')
    if ~matlabpool('size')
        matlabpool open
    end
    parallel_check          = true;
else
    parallel_check          = false;
end

%% draw GUI

set(0, 'DefaultFigureWindowStyle', 'docked')
if ispc % windows switch
    mgui                    = figure('toolbar', 'figure', 'name', 'MERGEGUI', 'position', [1920 940 1 1], 'menubar', 'none', 'keypressfcn', @keypress, 'closerequestfcn', @closefig);
    ax_radar                = subplot('position', [0.065 0.06 1.42 0.81]);
    size_font               = 14;
    width_slide             = 0.01;
else
    mgui                    = figure('toolbar', 'figure', 'name', 'MERGEGUI', 'position', [1864 1100 1 1], 'menubar', 'none', 'keypressfcn', @keypress, 'closerequestfcn', @closefig);
    ax_radar                = subplot('position', [0.065 0.06 0.86 0.81]);
    size_font               = 18;
    width_slide             = 0.02;
end

hold on
colormap(bone)
caxis([db_min db_max])
axis xy tight
set(gca, 'fontsize', size_font, 'layer', 'top')
xlabel('Distance (km)')
yl                          = ylabel('Elevation (m)');
colorbar('fontsize', size_font)
box on
% pan/zoom callbacks
h_pan                       = pan;
set(h_pan, 'actionpostcallback', @panzoom)
h_zoom                      = zoom;
set(h_zoom, 'actionpostcallback', @panzoom)

% sliders
z_min_slide                 = uicontrol(mgui, 'style', 'slider', 'units', 'normalized', 'position', [0.005 0.07 width_slide 0.32], 'callback', @slide_z_min, 'min', 0, 'max', 1, 'value', elev_min_ref, ...
                                              'sliderstep', [0.01 0.1]);
z_max_slide                 = uicontrol(mgui, 'style', 'slider', 'units', 'normalized', 'position', [0.005 0.50 width_slide 0.32], 'callback', @slide_z_max, 'min', 0, 'max', 1, 'value', elev_max_ref, ...
                                              'sliderstep', [0.01 0.1]);
cb_min_slide                = uicontrol(mgui, 'style', 'slider', 'units', 'normalized', 'position', [0.96 0.07 width_slide 0.32], 'callback', @slide_db_min, 'min', -150, 'max', 0, 'value', db_min_ref, ...
                                              'sliderstep', [0.01 0.1]);
cb_max_slide                = uicontrol(mgui, 'style', 'slider', 'units', 'normalized', 'position', [0.96 0.50 width_slide 0.32], 'callback', @slide_db_max, 'min', -150, 'max', 0, 'value', db_max_ref, ...
                                              'sliderstep', [0.01 0.1]);
dist_min_slide              = uicontrol(mgui, 'style', 'slider', 'units', 'normalized', 'position', [0.12 width_slide 0.27 0.02], 'callback', @slide_dist_min, 'min', 0, 'max', 1, 'value', dist_min_ref, ...
                                              'sliderstep', [0.01 0.1]);
dist_max_slide              = uicontrol(mgui, 'style', 'slider', 'units', 'normalized', 'position', [0.64 width_slide 0.27 0.02], 'callback', @slide_dist_max, 'min', 0, 'max', 1, 'value', dist_max_ref, ...
                                              'sliderstep', [0.01 0.1]);
% slider values
z_min_edit                  = annotation('textbox', [0.005 0.39 0.04 0.03], 'string', num2str(elev_min_ref), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
z_max_edit                  = annotation('textbox', [0.005 0.82 0.04 0.03], 'string', num2str(elev_max_ref), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
cb_min_edit                 = annotation('textbox', [0.965 0.39 0.04 0.03], 'string', num2str(db_min_ref), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
cb_max_edit                 = annotation('textbox', [0.9665 0.82 0.04 0.03], 'string', num2str(db_max_ref), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
dist_min_edit               = annotation('textbox', [0.08 0.005 0.04 0.03], 'string', num2str(dist_min_ref), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
dist_max_edit               = annotation('textbox', [0.59 0.005 0.04 0.03], 'string', num2str(dist_max_ref), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');

% push buttons
uicontrol(mgui, 'style', 'pushbutton', 'string', 'Load picks', 'units', 'normalized', 'position', [0.005 0.965 0.07 0.03], 'callback', @load_pk, 'fontsize', size_font, 'foregroundcolor', 'b')
uicontrol(mgui, 'style', 'pushbutton', 'string', 'Load radar data', 'units', 'normalized', 'position', [0.30 0.965 0.09 0.03], 'callback', @load_data, 'fontsize', size_font, 'foregroundcolor', 'b')
uicontrol(mgui, 'style', 'pushbutton', 'string', 'Load core intersections', 'units', 'normalized', 'position', [0.64 0.925 0.11 0.03], 'callback', @load_core, 'fontsize', size_font, 'foregroundcolor', 'b')
uicontrol(mgui, 'style', 'pushbutton', 'string', 'Flatten amplitude', 'units', 'normalized', 'position', [0.30 0.885 0.09 0.03], 'callback', @flatten, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(mgui, 'style', 'pushbutton', 'string', 'Focus', 'units', 'normalized', 'position', [0.41 0.925 0.05 0.03], 'callback', @focus_layer, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(mgui, 'style', 'pushbutton', 'string', 'Match', 'units', 'normalized', 'position', [0.41 0.885 0.04 0.03], 'callback', @pk_match, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(mgui, 'style', 'pushbutton', 'string', 'Last', 'units', 'normalized', 'position', [0.46 0.885 0.04 0.03], 'callback', @pk_last, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(mgui, 'style', 'pushbutton', 'string', 'Next', 'units', 'normalized', 'position', [0.46 0.925 0.04 0.03], 'callback', @pk_next, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(mgui, 'style', 'pushbutton', 'string', 'Choose layer', 'units', 'normalized', 'position', [0.505 0.925 0.07 0.03], 'callback', @choose_pk, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(mgui, 'style', 'pushbutton', 'string', 'Delete layer', 'units', 'normalized', 'position', [0.505 0.885 0.07 0.03], 'callback', @del_layer, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(mgui, 'style', 'pushbutton', 'string', 'test', 'units', 'normalized', 'position', [0.59 0.965 0.04 0.03], 'callback', @misctest, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(mgui, 'style', 'pushbutton', 'string', 'Update picks', 'units', 'normalized', 'position', [0.78 0.925 0.08 0.03], 'callback', @update_pk, 'fontsize', size_font, 'foregroundcolor', 'g')
uicontrol(mgui, 'style', 'pushbutton', 'string', 'Map', 'units', 'normalized', 'position', [0.865 0.925 0.03 0.03], 'callback', @pop_map, 'fontsize', size_font, 'foregroundcolor', 'g')
uicontrol(mgui, 'style', 'pushbutton', 'string', 'Pop figure', 'units', 'normalized', 'position', [0.90 0.925 0.06 0.03], 'callback', @pop_fig, 'fontsize', size_font, 'foregroundcolor', 'g')
uicontrol(mgui, 'style', 'pushbutton', 'string', 'Save', 'units', 'normalized', 'position', [0.965 0.925 0.03 0.03], 'callback', @pk_save, 'fontsize', size_font, 'foregroundcolor', 'g')
uicontrol(mgui, 'style', 'pushbutton', 'string', 'reset x/y', 'units', 'normalized', 'position', [0.945 0.885 0.05 0.03], 'callback', @reset_xz, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(mgui, 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.005 0.03 0.03 0.03], 'callback', @reset_z_min, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(mgui, 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.005 0.46 0.03 0.03], 'callback', @reset_z_max, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(mgui, 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.40 0.005 0.03 0.03], 'callback', @reset_dist_min, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(mgui, 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.92 0.005 0.03 0.03], 'callback', @reset_dist_max, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(mgui, 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.965 0.03 0.03 0.03], 'callback', @reset_db_min, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(mgui, 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.965 0.46 0.03 0.03], 'callback', @reset_db_max, 'fontsize', size_font, 'foregroundcolor', 'r')

% fixed text annotations
a(1)                        = annotation('textbox', [0.13 0.965 0.03 0.03], 'string', 'N_{decimate}', 'fontsize', size_font, 'color', 'b', 'edgecolor', 'none');
a(2)                        = annotation('textbox', [0.21 0.965 0.03 0.03], 'string', 'Chunk', 'fontsize', size_font, 'color', 'b', 'edgecolor', 'none');
a(3)                        = annotation('textbox', [0.21 0.925 0.03 0.03], 'string', 'L_{chunk}', 'fontsize', size_font, 'color', 'b', 'edgecolor', 'none');
a(4)                        = annotation('textbox', [0.965 0.42 0.03 0.03], 'string', 'dB_{min}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(5)                        = annotation('textbox', [0.965 0.85 0.03 0.03], 'string', 'dB_{max}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(6)                        = annotation('textbox', [0.41 0.965 0.03 0.03], 'string', 'Layer', 'fontsize', size_font, 'color', 'm', 'edgecolor', 'none');
a(7)                        = annotation('textbox', [0.85 0.885 0.03 0.03], 'string', 'Grid', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(8)                        = annotation('textbox', [0.005 0.42 0.03 0.03], 'string', 'z_{min}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(9)                        = annotation('textbox', [0.005 0.85 0.03 0.03], 'string', 'z_{max}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(10)                       = annotation('textbox', [0.04 0.005 0.03 0.03], 'string', 'dist_{min}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(11)                       = annotation('textbox', [0.55 0.005 0.03 0.03], 'string', 'dist_{max}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(12)                       = annotation('textbox', [0.025 0.85 0.03 0.03], 'string', 'fix', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(13)                       = annotation('textbox', [0.005 0.005 0.03 0.03], 'string', 'fix', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(14)                       = annotation('textbox', [0.95 0.005 0.03 0.03], 'string', 'fix 1', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(15)                       = annotation('textbox', [0.98 0.005 0.03 0.03], 'string', '2', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(16)                       = annotation('textbox', [0.765 0.885 0.05 0.03], 'string', 'dd-mmm-yyyy', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(17)                       = annotation('textbox', [0.31 0.925 0.10 0.03], 'string', 'block divisions', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
if ~ispc
    set(a, 'fontweight', 'bold')
end

% variable text annotations
file_box                    = annotation('textbox', [0.005 0.925 0.20 0.03], 'string', '', 'color', 'k', 'fontsize', size_font, 'backgroundcolor', 'w', 'edgecolor', 'k', 'interpreter', 'none');
date_box                    = annotation('textbox', [0.64 0.885 0.12 0.03], 'string', '', 'color', 'k', 'fontsize', size_font, 'backgroundcolor', 'w', 'edgecolor', 'k', 'interpreter', 'none');
status_box                  = annotation('textbox', [0.64 0.965 0.35 0.03], 'string', '', 'color', 'k', 'fontsize', size_font, 'backgroundcolor', 'w', 'edgecolor', 'k', 'interpreter', 'none');
cbl                         = annotation('textbox', [0.93 0.03 0.03 0.03], 'string', '(dB)', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
if ~ispc
    set(cbl, 'fontweight', 'bold')
end

% value boxes
decim_edit                  = uicontrol(mgui, 'style', 'edit', 'string', num2str(decim), 'units', 'normalized', 'position', [0.175 0.965 0.03 0.03], 'fontsize', size_font, 'foregroundcolor', 'k', ...
                                               'backgroundcolor', 'w', 'callback', @adj_decim);
length_chunk_edit           = uicontrol(mgui, 'style', 'edit', 'string', num2str(length_chunk), 'units', 'normalized', 'position', [0.25 0.925 0.04 0.03], 'fontsize', size_font, 'foregroundcolor', 'k', ...
                                               'backgroundcolor', 'w', 'callback', @adj_length_chunk);
% menus
chunk_list                  = uicontrol(mgui, 'style', 'popupmenu', 'string', 'N/A', 'value', 1, 'units', 'normalized', 'position', [0.25 0.955 0.045 0.04], 'fontsize', size_font, 'foregroundcolor', 'k', ...
                                               'callback', @plot_chunk);
cmap_list                   = uicontrol(mgui, 'style', 'popupmenu', 'string', cmaps, 'value', 1, 'units', 'normalized', 'position', [0.89 0.865 0.05 0.05], 'callback', @change_cmap, 'fontsize', size_font);
layer_list                  = uicontrol(mgui, 'style', 'popupmenu', 'string', 'N/A', 'value', 1, 'units', 'normalized', 'position', [0.45 0.955 0.05 0.04], 'fontsize', size_font, 'foregroundcolor', 'k', ...
                                               'callback', @choose_layer);
block_list                  = uicontrol(mgui, 'style', 'popupmenu', 'string', 'N/A', 'value', 1, 'units', 'normalized', 'position', [0.16 0.865 0.135 0.05], 'fontsize', size_font, 'foregroundcolor', 'k');

% check boxes
zfix_check                  = uicontrol(mgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.04 0.85 0.01 0.03], 'fontsize', size_font, 'value', 0);
distfix_check               = uicontrol(mgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.02 0.005 0.01 0.03], 'fontsize', size_font, 'value', 0);
cbfix_check1                = uicontrol(mgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.97 0.005 0.01 0.03], 'fontsize', size_font, 'value', 0);
cbfix_check2                = uicontrol(mgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.985 0.005 0.01 0.03], 'callback', @narrow_cb, 'fontsize', size_font, 'value', 0);
grid_check                  = uicontrol(mgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.88 0.885 0.01 0.03], 'callback', @toggle_grid, 'fontsize', size_font, 'value', 0);
pk_check                    = uicontrol(mgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.08 0.965 0.01 0.03], 'callback', @show_pk, 'fontsize', size_font, 'value', 0);
date_check                  = uicontrol(mgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.835 0.885 0.01 0.03], 'callback', @adj_date, 'fontsize', size_font, 'value', 0);
data_check                  = uicontrol(mgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.395 0.965 0.01 0.03], 'callback', @show_data, 'fontsize', size_font, 'value', 0);
block_check                 = uicontrol(mgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.395 0.925 0.01 0.03], 'callback', @show_block, 'fontsize', size_font, 'value', 0);
core_check                  = uicontrol(mgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.75 0.925 0.01 0.03], 'callback', @show_core, 'fontsize', size_font, 'value', 0);

% display buttons
disp_group                  = uibuttongroup('position', [0.005 0.885 0.15 0.03], 'selectionchangefcn', @disp_radio);
uicontrol(mgui, 'style', 'text', 'parent', disp_group, 'units', 'normalized', 'position', [0 0.6 0.9 0.3], 'fontsize', size_font)
disp_check(1)               = uicontrol(mgui, 'style', 'radio', 'string', 'amplitude', 'units', 'normalized', 'position', [0.01 0.1 0.45 0.8], 'parent', disp_group, 'fontsize', size_font, 'handlevisibility', 'off');
disp_check(2)               = uicontrol(mgui, 'style', 'radio', 'string', 'flattened', 'units', 'normalized', 'position', [0.51 0.1 0.45 0.8], 'parent', disp_group, 'fontsize', size_font, 'handlevisibility', 'off');
set(disp_group, 'selectedobject', disp_check(1))

%% Clear plots

    function clear_plots(source, eventdata) %#ok<*INUSD>
        if (logical(p_bed) && ishandle(p_bed))
            delete(p_bed)
        end
        if (logical(p_bedflat) && ishandle(p_bedflat))
            delete(p_bedflat)
        end
        if (any(p_block) && any(ishandle(p_block)))
            delete(p_block(logical(p_block) & ishandle(p_block)))
        end
        if (any(p_blockflat) && any(ishandle(p_blockflat)))
            delete(p_blockflat(logical(p_blockflat) & ishandle(p_blockflat)))
        end
        if (any(p_blocknum) && any(ishandle(p_blocknum)))
            delete(p_blocknum(logical(p_blocknum) & ishandle(p_blocknum)))
        end
        if (any(p_blocknumflat) && any(ishandle(p_blocknumflat)))
            delete(p_blocknumflat(logical(p_blocknumflat) & ishandle(p_blocknumflat)))
        end
        if (any(p_core) && any(ishandle(p_core)))
            delete(p_core(logical(p_core) & ishandle(p_core)))
        end
        if (any(p_coreflat) && any(ishandle(p_coreflat)))
            delete(p_coreflat(logical(p_coreflat) & ishandle(p_coreflat)))
        end
        if (any(p_corename) && any(ishandle(p_corename)))
            delete(p_corename(logical(p_corename) & ishandle(p_corename)))
        end
        if (any(p_corenameflat) && any(ishandle(p_corenameflat)))
            delete(p_corenameflat(logical(p_corenameflat) & ishandle(p_corenameflat)))
        end
        if (logical(p_data) && ishandle(p_data))
            delete(p_data)
        end
        if (any(p_pk) && any(ishandle(p_pk)))
            delete(p_pk(logical(p_pk) & ishandle(p_pk)))
        end
        if (any(p_pkflat) && any(ishandle(p_pkflat)))
            delete(p_pkflat(logical(p_pkflat) & ishandle(p_pkflat)))
        end
        if (logical(p_surf) && ishandle(p_surf))
            delete(p_surf)
        end
        set([layer_list block_list], 'string', 'N/A', 'value', 1)
        if (get(disp_group, 'selectedobject') ~= disp_check(1))
            set(disp_group, 'selectedobject', disp_check(1))
            disp_type       = 'amplitude';
        end
    end

%% Clear data and picks

    function clear_data(source, eventdata)
        pk                  = struct;
        [bed_avail, data_done, edit_flag, flat_done, merge_done, merge_file, surf_avail] ...
                            = deal(false);
        [amp_mean, amp_flat, colors, curr_chunk, curr_layer, curr_trans, curr_year, depth, depth_bed, depth_bed_flat, depth_curr, depth_flat, depth_layer_flat, depth_mat, dist_chunk, dt, elev, ii, ind_corr, ...
         ind_decim, ind_int, ind_layer_flat, ind_x_pk, ind_x_ref, ind_y_pk, jj, kk, num_chunk, num_data, num_decim, num_int, num_pk, num_sample, p_bed, p_block, p_blockflat, p_blocknum, p_blocknumflat, p_core, ...
         p_coreflat, p_corename, p_corenameflat, p_data, p_pk, p_pkflat, p_refflat, p_surf, pk_all, pkfig, tmp1, tmp2, tmp3, tmp4, tmp5, twtt] ...
                            = deal(0);
        [file_data, file_pk, file_pk_short, file_save] ...
                            = deal('');
        set([file_box date_box], 'string', '')
    end

%% Load picks for this transect

    function load_pk(source, eventdata)
        
        if (merge_done || data_done)
            clear_plots
            clear_data
        end
        
        % Dialog box to choose picks file to load
        set(mgui, 'keypressfcn', [])
        set(status_box, 'string', 'Load pick files (F) or directory (D)?')
        waitforbuttonpress
        if strcmpi(get(mgui, 'currentcharacter'), 'F')
            [tmp1, tmp2]    = deal(file_pk, path_pk);
            if ~isempty(path_pk)
                [file_pk, path_pk] = uigetfile('*.mat', 'Load picks files (*_pk.mat):', path_pk, 'multiselect', 'on');
            elseif ~isempty(path_data)
                [file_pk, path_pk] = uigetfile('*.mat', 'Load picks files: (*_pk.mat)', path_data, 'multiselect', 'on');
            else
                [file_pk, path_pk] = uigetfile('*.mat', 'Load picks files: (*_pk.mat)', 'multiselect', 'on');
            end
            if isnumeric(file_pk)
                [file_pk, path_pk] = deal('', tmp2);
            end
        elseif strcmpi(get(mgui, 'currentcharacter'), 'D')
            [tmp1, tmp2]    = deal(file_pk, path_pk);
            if ~isempty(path_pk)
                path_pk     = uigetdir(path_pk, 'Load picks directory (containing *_pk.mat):');
            elseif ~isempty(path_data)
                path_pk     = uigetdir(path_data, 'Load picks directory (containing *_pk.mat):');
            else
                path_pk     = uigetdir('Load picks directory (containing *_pk.mat):');
            end
            if ~ischar(path_pk)
                [file_pk, path_pk] ...
                            = deal('', tmp2);
            else
                if ispc
                    path_pk = [path_pk '\'];
                else
                    path_pk = [path_pk '/'];
                end
                file_pk     = dir([path_pk '*_pk.mat']);
                file_pk     = {file_pk.name}';
            end
        else
            set(status_box, 'string', 'Choice unclear. Try again.')
            set(mgui, 'keypressfcn', @keypress)
            return
        end
        
        set(mgui, 'keypressfcn', @keypress)
        
        if ~isempty(file_pk)
            
            set(status_box, 'string', 'Starting initial loading and merging of picks...')
            pause(0.1)
            
            set(date_check, 'value', 0)
            
            if ischar(file_pk)
                file_pk     = {file_pk}; % only one pick file
            end
            
            % load picks files
            num_pk          = length(file_pk);
            pk_all          = cell(1, num_pk);
            [surf_avail, bed_avail] ...
                            = deal(false(1, num_pk));
            
            for ii = 1:num_pk %#ok<*FXUP>
                set(status_box, 'string', ['Loading ' file_pk{ii}(1:(end - 4)) ' (' num2str(ii) ' / ' num2str(num_pk) ')...'])
                pause(0.1)
                tmp1        = load([path_pk file_pk{ii}]);
                try
                    pk_all{ii} ...
                            = tmp1.pk;
                    tmp1    = 0;
                    % check to see if surface and bed picks are available
                    if isfield(pk_all{ii}, 'elev_surf')
                        surf_avail(ii) ...
                            = true;
                    end
                    if isfield(pk_all{ii}, 'elev_bed')
                        bed_avail(ii) ...
                            = true;
                    end
                catch % give up, force restart
                    set(status_box, 'string', [file_pk{ii} ' does not contain a pk structure. Try again.'])
                    [pk_all, num_pk] ...
                            = deal(0);
                    return
                end
            end
            
            if (~any(surf_avail) && any(~isnan(pk.elev_surf)))
                set(status_box, 'string', 'Surface pick is not included in the pick files, so data will not be elevation-corrected.')
                pause(0.1)
            end
            
            % extract date and simple name from pk files
            [tmp1, tmp2]    = strtok(file_pk{ii}, '_');
            file_pk_short   = [tmp1 tmp2(1:3)];
            if (~strcmp(tmp2(4), '_') && ~strcmp(tmp2(4:8), 'block'))
                file_pk_short ...
                            = [file_pk_short tmp2(4)];
            end
            set(file_box, 'string', file_pk_short)
            pk.date_num     = datenum([str2double(tmp1(1:4)) str2double(tmp1(5:6)) str2double(tmp1(7:8))]);
            pk.date_str     = datestr(pk.date_num);
            set(date_box, 'string', pk.date_str)
            
            if ((num_pk == 1) && isfield(pk_all{1}, 'merge_flag')) % data already merged
                
                set(date_check, 'value', 1)
                tmp1        = fieldnames(pk_all{1});
                for ii = 1:length(tmp1)
                    eval(['pk.' tmp1{ii} ' = pk_all{1}.' tmp1{ii} ';'])
                end
                merge_file  = true;
                path_save   = path_pk;
                if isfield(pk, 'poly_flat_merge')
                    if any(pk.poly_flat_merge(:))
                        flat_done ...
                            = true;
                    end
                end
                if isfield(pk, 'ind_layer_keep_phase')
                    pk.ind_keep_phase ...
                        = pk.ind_layer_keep_phase;
                end
                if isfield(pk, 'ind_layer_keep')
                    pk.ind_keep_phase ...
                        = pk.ind_layer_keep;
                end
                if isfield(pk, 'num_layer_keep_phase')
                    pk.num_keep_phase ...
                        = pk.num_layer_keep_phase;
                end
                if isfield(pk, 'num_layer_keep')
                    pk.num_keep_phase ...
                        = pk.num_layer_keep;
                end
                if isfield(pk, 'ind_y_start')
                    pk.ind_y_start_phase ...
                        = pk.ind_y_start;
                end
                if isfield(pk, 'num_layer_keep_aresp')
                    pk.num_keep_aresp ...
                        = pk.num_layer_keep_aresp;
                end
                
                set(status_box, 'string', 'Merged file loaded...')
                pause(0.1)
                
            else % time to merge
                
                for ii = 1:num_cell
                    eval(['pk.' var_cell{ii} ' = cell(num_pk, 1);'])
                end
                for ii = 1:num_zero
                    eval(['pk.' var_zero{ii} ' = zeros(num_pk, 1);'])
                end
                pk.ind_overlap  ...
                            = NaN(num_pk, 2);
                tmp1        = cell(1, num_pk);
                
                % loop through pk files and assign values
                for ii = 1:num_pk
                    pk.ind_overlap(ii, :) ...
                            = pk_all{ii}.ind_overlap;
                    pk.file_pk{ii} ...
                            = file_pk{ii};
                    pk.file_block{ii} ...
                            = pk_all{ii}.file_block;
                    pk.num_layer_block(ii) ...
                            = pk_all{ii}.num_layer;
                    for jj = 1:(num_cell - 3)
                        if isfield(pk_all{ii}, var_cell{jj})
                            eval(['pk.' var_cell{jj} '{ii} = pk_all{ii}.' var_cell{jj} ';'])
                        end
                    end
                    if isfield(pk_all{ii}, 'ind_layer_keep')
                        pk.ind_keep_phase{ii} ...
                            = pk_all{ii}.ind_layer_keep;
                    end
                    if isfield(pk_all{ii}, 'ind_layer_keep_phase')
                        pk.ind_keep_phase{ii} ...
                            = pk_all{ii}.ind_layer_keep_phase;
                    end
                    if isfield(pk_all{ii}, 'ind_layer_keep_aresp')
                        pk.ind_keep_aresp{ii} ...
                            = pk_all{ii}.ind_layer_keep_aresp;
                    end
                    if isfield(pk_all{ii}, 'ind_y_start')
                        pk.ind_y_start_phase{ii} ...
                            = pk_all{ii}.ind_y_start;
                    end
                    for jj = 1:(num_zero - 2)
                        if isfield(pk_all{ii}, var_zero{jj})
                            eval(['pk.' var_zero{jj} '(ii) = pk_all{ii}.' var_zero{jj} ';'])
                        end
                    end
                    if isfield(pk_all{ii}, 'num_layer_keep')
                        pk.num_keep_phase(ii) ...
                            = pk_all{ii}.num_layer_keep;
                    end
                    if isfield(pk_all{ii}, 'num_layer_keep_aresp')
                        pk.num_keep_aresp(ii) ...
                            = pk_all{ii}.num_layer_keep_aresp;
                    end
                    if isfield(pk_all{ii}, 'ind_x_start')
                        pk.ind_x_start_phase(ii) ...
                            = pk_all{ii}.ind_x_start;
                    end
                    if isfield(pk_all{ii}, 'ind_y_flat_mean')
                        pk.ind_y_flat_mean{ii} ...
                            = NaN(pk_all{ii}.num_layer, length(pk_all{ii}.layer(1).ind_y_flat_mean));
                        for jj = 1:pk_all{ii}.num_layer
                            pk.ind_y_flat_mean{ii}(jj, :) ...
                            = pk_all{ii}.layer(jj).ind_y_flat_mean;
                        end
                    end
                    tmp1{ii}= pk_all{ii}.ind_match';
                end
                
                pk.num_layer= length(unique([tmp1{:}])); % number of layers in this transect
                
                pk.num_trace_tot ...
                            = sum(pk.num_trace) - sum(pk.ind_overlap(~isnan(pk.ind_overlap(:, 1)), 1)); % number of non-overlapping traces
                if ~isnan(pk.ind_overlap(1, 1)) % true if not starting at beginning of transect (e.g., no layers in early blocks)
                    pk.num_trace_tot ...
                            = pk.num_trace_tot + pk.ind_overlap(1, 1);
                end
                
                % allocate variables for merging
                for ii = 1:num_var_layer
                    eval(['pk.' var_layer{ii} ' = NaN(pk.num_layer, pk.num_trace_tot, ''single'');'])
                end
                for ii = 1:num_var_pos
                    eval(['pk.' var_pos{ii} ' = NaN(1, pk.num_trace_tot);'])
                end
                
                [pk.ind_trace_start(1), pk.ind_layer_all, pk.ind_match_block] ...
                            = deal(1, 1:pk.num_layer, false(pk.num_layer, num_pk)); % start at 1, layer indices that are present in all blocks (to be eroded in loop below)
                
                % merge picks
                for ii = 1:num_pk
                    set(status_box, 'string', ['Merging picks (' num2str(ii) ' / ' num2str(num_pk) ')...'])
                    pause(0.1)
                    tmp1    = pk.ind_trace_start(ii):(pk.ind_trace_start(ii) + pk.num_trace(ii) - 1); % merged indices
                    
                    % merge layer data
                    for jj = 1:pk.num_layer
                        if any(pk.ind_match{ii} == jj)
                            tmp2                        = find((pk.ind_match{ii} == jj), 1);
                            pk.ind_match_block(jj, ii)  = true;
                            if (ii == 1)
                                for kk = 1:num_var_layer
                                    eval(['pk.' var_layer{kk} '(jj, tmp1) = single(pk_all{ii}.layer(tmp2).' var_layer{kk} ');'])
                                end
                            else
                                for kk = 1:num_var_layer
                                    eval(['pk.' var_layer{kk} '(jj, tmp1(1:pk.ind_overlap(ii, 1))) = single(nanmean([pk_all{ii}.layer(tmp2).' var_layer{kk} '(1:pk.ind_overlap(ii, 1)); pk.' var_layer{kk} ...
                                          '(jj, tmp1(1:pk.ind_overlap(ii, 1)))]));'])
                                    tmp3                = pk_all{ii}.layer(tmp2); % intermediate step because of some variable assignment restriction
                                    eval(['pk.' var_layer{kk} '(jj, tmp1((pk.ind_overlap(ii, 1) + 1):end)) = single(tmp3.' var_layer{kk} '((pk.ind_overlap(ii, 1) + 1):end));'])
                                end
                            end
                        else
                            pk.ind_layer_all            = setdiff(pk.ind_layer_all, jj);
                        end
                    end
                    
                    % merge position data
                    if (ii == 1)
                        if (surf_avail(ii) && bed_avail(ii))
                            for jj = 1:num_var_pos
                                eval(['pk.' var_pos{jj} '(tmp1) = pk_all{ii}.' var_pos{jj} ';'])
                            end
                        elseif (surf_avail(ii) && ~bed_avail(ii))
                            for jj = 1:(num_var_pos - 2)
                                eval(['pk.' var_pos{jj} '(tmp1) = pk_all{ii}.' var_pos{jj} ';'])
                            end
                        elseif (~surf_avail(ii) && bed_avail(ii))
                            for jj = [1:6 9 10]
                                eval(['pk.' var_pos{jj} '(tmp1) = pk_all{ii}.' var_pos{jj} ';'])
                            end
                        else
                            for jj = 1:(num_var_pos - 4)
                                eval(['pk.' var_pos{jj} '(tmp1) = pk_all{ii}.' var_pos{jj} ';'])
                            end
                        end
                    else % should not need to average overlapping position data
                        if (surf_avail(ii) && bed_avail(ii))
                            for jj = 1:num_var_pos
                                eval(['pk.' var_pos{jj} '(tmp1((pk.ind_overlap(ii, 1) + 1):end)) = pk_all{ii}.' var_pos{jj} '((pk.ind_overlap(ii, 1) + 1):end);'])
                            end
                        elseif (surf_avail(ii) && ~bed_avail(ii))
                            for jj = 1:(num_var_pos - 2)
                                eval(['pk.' var_pos{jj} '(tmp1((pk.ind_overlap(ii, 1) + 1):end)) = pk_all{ii}.' var_pos{jj} '((pk.ind_overlap(ii, 1) + 1):end);'])
                            end
                        elseif (~surf_avail(ii) && bed_avail(ii))
                            for jj = [1:6 9 10]
                                eval(['pk.' var_pos{jj} '(tmp1((pk.ind_overlap(ii, 1) + 1):end)) = pk_all{ii}.' var_pos{jj} '((pk.ind_overlap(ii, 1) + 1):end);'])
                            end
                        else
                            for jj = 1:(num_var_pos - 4)
                                eval(['pk.' var_pos{jj} '(tmp1((pk.ind_overlap(ii, 1) + 1):end)) = pk_all{ii}.' var_pos{jj} '((pk.ind_overlap(ii, 1) + 1):end);'])
                            end
                        end
                    end
                    
                    % prepare for next block and correct for overlap
                    if (ii < num_pk)
                        pk.ind_trace_start(ii + 1)      = pk.ind_trace_start(ii) + pk.num_trace(ii);
                        if ~isnan(pk.ind_overlap(ii, 2));
                            pk.ind_trace_start(ii + 1)  = pk.ind_trace_start(ii + 1) - (pk.num_trace(ii) - pk.ind_overlap(ii, 2) + 1);
                        end
                    end
                end
                
                % sort layers by decreasing elevation
                tmp1        = zeros(pk.num_layer, 1);
                for ii = 1:pk.num_layer
                    tmp1(ii)= nanmean(pk.elev_smooth(ii, :));
                end
                [~, tmp1]   = sort(tmp1);
                for ii = 1:num_var_layer
                    eval(['pk.' var_layer{ii} ' = pk.' var_layer{ii} '(flipud(tmp1), :);'])
                end
                
                pk.dist_lin = interp1([1 pk.num_trace_tot], pk.dist([1 end]), 1:pk.num_trace_tot); % new linear distance vector (should not be merged because each linear distance vector is different)
                
            end
            
            % get rid of weird empty layers
            tmp1            = [];
            for ii = 1:pk.num_layer
                if all(isnan(pk.ind_y(ii, :)))
                    tmp1    = [tmp1 ii]; %#ok<AGROW>
                end
            end
            if ~isempty(tmp1)
                for ii = 1:num_var_layer
                    eval(['pk.' var_layer{ii} ' = pk.' var_layer{ii} '(setdiff(1:pk.num_layer, tmp1), :);'])
                end
                pk.num_layer= pk.num_layer - length(tmp1);
            end
            
            [pk_all, tmp3]  = deal(0);
            
            % decimated vectors for display
            if (decim > 1)
                ind_decim   = (1 + ceil(decim / 2)):decim:(pk.num_trace_tot - ceil(decim / 2));
            else
                ind_decim   = 1:pk.num_trace_tot;
            end
            num_decim       = length(ind_decim);
            
            if (merge_file && flat_done && isfield(pk, 'ind_x_ref'))
                ind_x_ref   = interp1(ind_decim, 1:num_decim, pk.ind_x_ref, 'nearest', 'extrap');
            else
                [ind_x_ref, pk.ind_x_ref] ...
                            = deal(1);
            end
            
            % make chunks
            adj_length_chunk
            
            % display merged picks
            colors          = repmat(colors_def, ceil(pk.num_layer / size(colors_def, 1)), 1); % extend predefined color pattern
            colors          = colors(1:pk.num_layer, :);
            p_pk            = deal(zeros(1, pk.num_layer));
            layer_str       = num2cell(1:pk.num_layer);
            for ii = 1:pk.num_layer
                if ~any(~isnan(pk.elev_smooth(ii, ind_decim)))
                    p_pk(ii)=plot(0, 0, 'w.', 'markersize', 1, 'visible', 'off');
                    layer_str{ii} ...
                            = [num2str(layer_str{ii}) ' H'];
                else
                    p_pk(ii)= plot(pk.dist_lin(ind_decim(~isnan(pk.elev_smooth(ii, ind_decim)))), pk.elev_smooth(ii, ind_decim(~isnan(pk.elev_smooth(ii, ind_decim)))), ...
                                   '.', 'color', colors(ii, :), 'markersize', 12, 'visible', 'off');
                end
            end
            set(layer_list, 'string', layer_str, 'value', 1)
            
            [dist_min_ref, dist_max_ref, dist_min, dist_max] ...
                            = deal(pk.dist_lin(1), pk.dist_lin(end), pk.dist_lin(1), pk.dist_lin(end));
            if (any(surf_avail) && any(~isnan(pk.elev_surf)))
                [elev_max_ref, elev_max] ...
                            = deal(max(pk.elev_surf(~isinf(pk.elev_surf))) + (0.1 * (max(pk.elev_surf(~isinf(pk.elev_surf))) - min(pk.elev_surf(~isinf(pk.elev_surf))))));
            else
                [elev_max_ref, elev_max] ...
                            = deal(max(pk.elev_smooth(:)) + (0.1 * (max(pk.elev_smooth(:)) - min(pk.elev_smooth(:)))));
            end
            if (any(bed_avail) && any(~isnan(pk.elev_bed)))
                [elev_min_ref, elev_min] ...
                            = deal(min(pk.elev_bed(~isinf(pk.elev_bed))) - (0.1 * (max(pk.elev_bed(~isinf(pk.elev_bed))) - min(pk.elev_bed(~isinf(pk.elev_bed))))));
            else
                [elev_min_ref, elev_min] ...
                            = deal(min(pk.elev_smooth(:)) - (0.1 * (max(pk.elev_smooth(:)) - min(pk.elev_smooth(:)))));
            end
            set(z_min_slide, 'min', elev_min_ref, 'max', elev_max_ref, 'value', elev_min_ref)
            set(z_max_slide, 'min', elev_min_ref, 'max', elev_max_ref, 'value', elev_max_ref)
            set(dist_min_slide, 'min', dist_min_ref, 'max', dist_max_ref, 'value', dist_min_ref)
            set(dist_max_slide, 'min', dist_min_ref, 'max', dist_max_ref, 'value', dist_max_ref)
            set(z_min_edit, 'string', sprintf('%4.0f', elev_min_ref))
            set(z_max_edit, 'string', sprintf('%4.0f', elev_max_ref))
            set(dist_min_edit, 'string', sprintf('%3.1f', dist_min_ref))
            set(dist_max_edit, 'string', sprintf('%3.1f', dist_max_ref))
            update_dist_range
            update_z_range
            
            [p_block, p_blocknum] ...
                            = deal(zeros(1, length(pk.ind_trace_start)));
            for ii = 1:length(pk.ind_trace_start)
                p_block(ii) = plot(repmat(pk.dist_lin(pk.ind_trace_start(ii)), 1, 2), [elev_min_ref elev_max_ref], 'm--', 'linewidth', 1, 'visible', 'off');
                if ~isempty(pk.file_block{ii})
                    p_blocknum(ii) ...
                            = text((pk.dist_lin(pk.ind_trace_start(ii)) + 1), (elev_max_ref - 100), pk.file_block{ii}((end - 1):end), 'color', 'm', 'fontsize', (size_font - 2), 'visible', 'off');
                end
            end
            
            set(pk_check, 'value', 1)
            merge_done      = true;
            set(disp_group, 'selectedobject', disp_check(1))
            disp_type       = 'amplitude';
            axes(ax_radar)
            axis xy
            show_pk
            show_block
            if core_done
                load_core_breakout
            end
            
            if merge_file
                set(status_box, 'string', ['Loaded merged pick file from ' file_pk_short '.'])
            else
                set(status_box, 'string', ['Loaded and merged ' num2str(num_pk) ' pick files from ' file_pk_short '.'])
            end
            
            % list with block filenames to load
            set(block_list, 'string', pk.file_block)
            
        else
            file_pk         = tmp1;
            set(status_box, 'string', 'No picks loaded.')
        end
    end

%% Load radar data

    function load_data(source, eventdata)
        
        if ~merge_done
            set(status_box, 'string', 'Load picks before data.')
            return
        end
        
        % check if data are in expected location based on picks' path
        tmp1                = file_pk_short;
        if isnan(str2double(tmp1(end))) % check for a/b/c/etc in file_pk_short
            tmp1            = tmp1(1:(end - 1));
        end
        
        if merge_file
            if ispc
                if (~isempty(path_pk) && exist([path_pk '..\block\' tmp1], 'dir'))
                    path_data = [path_pk(1:strfind(path_pk, '\merge')) 'block\' tmp1 '\'];
                end
            else
                if (~isempty(path_pk) && exist([path_pk '../block/' tmp1], 'dir'))
                    path_data = [path_pk(1:strfind(path_pk, '/merge')) 'block/' tmp1 '/'];
                end
            end
        else
            if ispc
                if (~isempty(path_pk) && (exist([path_pk '..\..\block\' tmp1], 'dir') || exist([path_pk '..\..\..\block\' tmp1], 'dir')))
                    path_data = [path_pk(1:strfind(path_pk, '\pk')) 'block\' tmp1 '\'];
                end
            else
                if (~isempty(path_pk) && (exist([path_pk '../../block/' file_pk_short], 'dir') || exist([path_pk '../../../block/' tmp1], 'dir')))
                    path_data = [path_pk(1:strfind(path_pk, '/pk')) 'block/' tmp1 '/'];
                end
            end
        end
        
        % check filenames are available
        tmp1                = false(length(pk.file_block), 1);
        for ii = 1:length(pk.file_block)
            if ~isempty(pk.file_block{ii})
                tmp1(ii)    = true;
            end
        end
        
        if (~isempty(path_data) && all(tmp1)) % get data filenames automatically if path found
            if (merge_file && exist([path_data pk.file_block{1} '.mat'], 'file'))
                file_data   = pk.file_block;
                for ii = 1:num_pk
                    file_data{ii} ...
                            = [pk.file_block{ii} '.mat'];
                end
            elseif exist([path_data file_pk{1}(1:(end - 7)) '.mat'], 'file')
                file_data   = cell(num_pk, 1);
                for ii = 1:num_pk
                    file_data{ii} ...
                            = [file_pk{ii}(1:(end - 7)) '.mat'];
                end
            end
        else % dialog box to choose radar data file to load
            [tmp1, tmp2]    = deal(file_data, path_data);
            if ~isempty(path_pk)
                [file_data, path_data] = uigetfile('*.mat', 'Load radar data:', path_pk, 'multiselect', 'on');
            else
                [file_data, path_data] = uigetfile('*.mat', 'Load radar data:', 'multiselect', 'on');
            end
            if isnumeric(file_data)
                [file_data, path_data] = deal('', tmp2);
            end
        end
        
        if ~isempty(file_data)
            
            pause(0.1)
            
            if ischar(file_data)
                file_data   = {file_data};
            end
            
            num_data        = length(file_data);
            if ((num_data ~= num_pk) && ~merge_file)
                set(status_box, 'string', ['Number of data blocks (' num2str(num_data) ') does not match number of picks files (' num2str(num_pk) ').'])
                return
            end
            
            % get rid of all data handle if somehow present
            if data_done
                data_done   = false;
            end
            
            load_data_breakout
            
        else
            file_data       = tmp1;
            set(status_box, 'string', 'No radar data loaded.')
        end
    end

%% Loading function broken out (so that changing decimation is automatic)

    function load_data_breakout(source, eventdata)
        
        if (logical(p_data) && ishandle(p_data))
            delete(p_data)
        end
        if (any(p_block) && any(ishandle(p_block)))
            delete(p_block(logical(p_block) & ishandle(p_block)))
        end
        if (any(p_blocknum) && any(ishandle(p_blocknum)))
            delete(p_blocknum(logical(p_blocknum) & ishandle(p_blocknum)))
        end
        
        % attempt to load data
        set(status_box, 'string', 'Loading radar data...')
        pause(0.1)
        
        for ii = 1:num_data
            
            set(status_box, 'string', ['Loading ' file_data{ii}(1:(end - 4)) ' (' num2str(ii) ' / ' num2str(num_data) ')...'])
            pause(0.1)
            tmp1            = load([path_data file_data{ii}]);
            try
                tmp1        = tmp1.block;
            catch %#ok<*CTCH>
                set(status_box, 'string', [file_data{ii} ' does not contain a block structure. Try again.'])
                return
            end
            
            if (ii == 1) % start decimation
                
                amp_mean    = NaN(pk.num_sample(ii), num_decim, 'single');
                [dt, twtt, num_sample] ...
                            = deal(tmp1.dt, tmp1.twtt, tmp1.num_sample);
                tmp2        = (1 + ceil(decim / 2)):decim:(pk.num_trace(ii) - ceil(decim / 2));
                tmp3        = floor(decim / 2);
                for jj = 1:length(tmp2)
                    amp_mean(:, jj) ...
                            = nanmean(tmp1.amp(:, (tmp2(jj) - tmp3):(tmp2(jj) + tmp3)), 2);
                end
                
            else % middle/end
                
                tmp2        = repmat((pk.ind_trace_start(ii) + pk.ind_overlap(ii, 1)), 1, 2);
                tmp2(2)     = tmp2(2) + pk.num_trace(ii) - pk.ind_overlap(ii, 1) - 1;
                if (ii < num_data)
                    tmp3    = (find((tmp2(1) > ind_decim), 1, 'last') + 1):find((tmp2(2) <= ind_decim), 1);
                else
                    tmp3    = (find((tmp2(1) > ind_decim), 1, 'last') + 1):num_decim;
                end
                if (ii == 2)
                    tmp3    = [(tmp3(1) - 2) (tmp3(1) - 1) tmp3]; %#ok<AGROW> % correct weirdness for first block
                end
                tmp2        = ind_decim(tmp3) - pk.ind_trace_start(ii) + 1;
                tmp4        = tmp2 - floor(decim / 2);
                tmp5        = tmp2 + floor(decim / 2);
                tmp4(tmp4 < 1) ...
                            = 1;
                tmp5(tmp5 > pk.num_trace(ii)) ...
                            = pk.num_trace(ii);
                for jj = 1:length(tmp3)
                    amp_mean(:, tmp3(jj)) ...
                            = nanmean(tmp1.amp(:, tmp4(jj):tmp5(jj)), 2);
                end
            end
            
            tmp1            = 0;
            
        end
        
        % convert to dB and determine elevation/depth vectors
        amp_mean(isinf(amp_mean)) ...
                            = NaN;
        amp_mean            = 10 .* log10(abs(amp_mean));
        num_sample          = size(amp_mean, 1);
        depth               = (speed_ice / 2) .* twtt;
        tmp1                = NaN(size(amp_mean), 'single');
        
        if (any(surf_avail) && any(~isnan(pk.elev_surf)))
            tmp2            = interp1(twtt, 1:num_sample, pk.twtt_surf(ind_decim), 'nearest', 'extrap'); % surface traveltime indices
            if any(isnan(tmp2))
                tmp2(isnan(tmp2)) ...
                            = round(interp1(find(~isnan(tmp2)), tmp2(~isnan(tmp2)), find(isnan(tmp2)), 'linear', 'extrap'));
            end
        else
            tmp2            = ones(1, num_decim);
        end
        tmp3                = pk.elev_surf(ind_decim);
        if any(isnan(tmp3))
            tmp3(isnan(tmp3)) ...
                            = interp1(find(~isnan(tmp3)), tmp3(~isnan(tmp3)), find(isnan(tmp3)), 'linear', 'extrap');
        end
        
        for ii = 1:num_decim
            tmp1(1:(num_sample - tmp2(ii) + 1), ii) ...
                            = amp_mean(tmp2(ii):num_sample, ii); % shift data up to surface
        end
        [amp_mean, ind_corr]= topocorr(tmp1, depth, tmp3); % topographically correct data
        tmp1                = 0;
        amp_mean            = flipud(amp_mean); % flip for axes
        ind_corr            = max(ind_corr) - ind_corr + 1;
        depth               = (speed_ice / 2) .* (0:dt:((num_sample - 1) * dt))'; % simple monotonically increasing depth vector
        elev                = flipud(max(pk.elev_surf(ind_decim)) - depth); % elevation vector
        depth               = depth(1:(num_sample - max(ind_corr)));
        
        if any(surf_avail & bed_avail)
            depth_bed       = pk.elev_surf(ind_decim) - pk.elev_bed(ind_decim);
        end
        
        % assign traveltime and distance reference values/sliders based on data
        [elev_min_ref, elev_max_ref, db_min_ref, db_max_ref, elev_min, elev_max, db_min, db_max, depth_min_ref, depth_max_ref, depth_min, depth_max] ...
                            = deal(min(elev), max(elev), min(amp_mean(~isinf(amp_mean(:)) & ~isnan(amp_mean(:)))), max(amp_mean(~isinf(amp_mean(:)) & ~isnan(amp_mean(:)))), min(elev), max(elev), ...
                                   min(amp_mean(~isinf(amp_mean(:)) & ~isnan(amp_mean(:)))), max(amp_mean(~isinf(amp_mean(:)) & ~isnan(amp_mean(:)))), min(depth), max(depth), min(depth), max(depth));
        set(z_min_slide, 'min', elev_min_ref, 'max', elev_max_ref, 'value', elev_min_ref)
        set(z_max_slide, 'min', elev_min_ref, 'max', elev_max_ref, 'value', elev_max_ref)
        set(cb_min_slide, 'min', db_min_ref, 'max', db_max_ref, 'value', db_min_ref)
        set(cb_max_slide, 'min', db_min_ref, 'max', db_max_ref, 'value', db_max_ref)
        set(z_min_edit, 'string', sprintf('%4.0f', elev_min_ref))
        set(z_max_edit, 'string', sprintf('%4.0f', elev_max_ref))
        set(cb_min_edit, 'string', sprintf('%3.0f', db_min_ref))
        set(cb_max_edit, 'string', sprintf('%3.0f', db_max_ref))
        update_z_range
        
        % plot block divisions
        [p_block, p_blocknum] ...
                            = deal(zeros(1, length(pk.ind_trace_start)));
        for ii = 1:length(pk.ind_trace_start)
            p_block(ii)     = plot(repmat(pk.dist_lin(pk.ind_trace_start(ii)), 1, 2), [elev_min_ref elev_max_ref], 'm--', 'linewidth', 1, 'visible', 'off');
            if ~isempty(pk.file_block{ii})                
                p_blocknum(ii) ...
                            = text((pk.dist_lin(pk.ind_trace_start(ii)) + 1), (elev_max_ref - 100), pk.file_block{ii}((end - 1):end), 'color', 'm', 'fontsize', (size_font - 2), 'visible', 'off');
            end
        end
        
        % plot data
        data_done           = true;
        set(data_check, 'value', 1)
        plot_db
        
        if flat_done
            ind_layer_flat  = NaN(pk.num_layer, 1);
            ind_layer_flat(~isnan(pk.elev_smooth(:, pk.ind_x_ref))) ...
                            = 1;
            % fix polynomials to current decimation vector
            if (num_decim ~= size(pk.poly_flat_merge, 2))
                pk.poly_flat_merge ...
                            = interp2(pk.poly_flat_merge, linspace(1, num_decim, num_decim), (1:size(pk.poly_flat_merge, 1))');
            end
            flatten_breakout
        end
        
        set(status_box, 'string', 'Transect radar data loaded.')
        
    end

%% Flatten amplitudes using layers

    function flatten(source, eventdata)
        
        if (~merge_done || ~data_done)
            set(status_box, 'string', 'Picks must be loaded and merged and data must be loaded prior to flattening.')
            return
        end
        if (pk.num_layer < 4)
            set(status_box, 'string', 'Not enough layers to flatten (need 4+ including surface).')
            return
        end
        
        flat_done           = false;
        
        if (logical(p_bedflat) && ishandle(p_bedflat))
            delete(p_bedflat)
        end
        if (any(p_pkflat) && any(ishandle(p_pkflat)))
            delete(p_pkflat(logical(p_pkflat) & ishandle(p_pkflat)))
        end
        if (any(p_blockflat) && any(ishandle(p_blockflat)))
            delete(p_blockflat(logical(p_blockflat) & ishandle(p_blockflat)))
        end
        if (any(p_blocknumflat) && any(ishandle(p_blocknumflat)))
            delete(p_blocknumflat(logical(p_blocknumflat) & ishandle(p_blocknumflat)))
        end
        if (any(p_coreflat) && any(ishandle(p_coreflat)))
            delete(p_coreflat(logical(p_coreflat) & ishandle(p_coreflat)))
        end
        if (any(p_corenameflat) && any(ishandle(p_corenameflat)))
            delete(p_corenameflat(logical(p_corenameflat) & ishandle(p_corenameflat)))
        end
        
        tic
        
        % compile layer depths to use in flattening
        pk.poly_flat_merge  = NaN(4, num_decim, 'single');
        depth_curr          = [zeros(1, num_decim); pk.depth_smooth(:, ind_decim)];
        
        % keep the maximum number of layers for first go
        ind_x_ref           = sum(isnan(depth_curr));
        ind_x_ref           = find((ind_x_ref == min(ind_x_ref)), 1); % reference trace, i.e., the trace with the most layers
        ind_layer_flat      = depth_curr(:, ind_x_ref); % layer depths at reference trace
        tmp1                = depth_curr(~isnan(ind_layer_flat), :); % depth of all layers that are not nan at the reference trace
        tmp2                = tmp1(:, ind_x_ref); % depths of non-nan layers at reference trace
        tmp3                = find(sum(~isnan(tmp1)) > 3); % traces where it will be worth doing the polynomial
        
        set(status_box, 'string', ['Starting polynomial fits using ' num2str(length(tmp2)) ' / ' num2str(pk.num_layer) ' layers...'])
        pause(0.1)
        
        set(mgui, 'keypressfcn', [])
        
        if parallel_check
            pctRunOnAll warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale')
            pctRunOnAll warning('off', 'MATLAB:polyfit:PolyNotUnique')
        else
            warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale')
            warning('off', 'MATLAB:polyfit:PolyNotUnique')
        end
        
        % now polyfit
        if parallel_check
            tmp1            = tmp1(:, tmp3);
            tmp4            = pk.poly_flat_merge(:, tmp3);
            parfor ii = 1:length(tmp3)
                tmp4(:, ii) = polyfit(tmp2(~isnan(tmp1(:, ii))), tmp1(~isnan(tmp1(:, ii)), ii), 3)'; %#ok<PFBNS>
            end
            pk.poly_flat_merge(:, tmp3) ...
                            = tmp4;
            tmp4            = 0;
        else
            for ii = 1:length(tmp3)
                pk.poly_flat_merge(:, tmp3(ii)) ...
                            = polyfit(tmp2(~isnan(tmp1(:, tmp3(ii)))), tmp1(~isnan(tmp1(:, tmp3(ii))), tmp3(ii)), 3)';
            end
        end
        
        set(status_box, 'string', 'Smoothing initial polynomials...')
        pause(0.1)
        
        % smooth polynomials
        tmp2                = round(mean(pk.length_smooth) / nanmean(diff(pk.dist(ind_decim))));
        if parallel_check
            tmp1            = pk.poly_flat_merge;
            parfor ii = 1:4
                tmp1(ii, :) = smooth_lowess(tmp1(ii, :), tmp2);
            end
            pk.poly_flat_merge ...
                            = tmp1;
            tmp1            = 0;
        else
            for ii = 1:4
                pk.poly_flat_merge(ii, :) ...
                            = smooth_lowess(pk.poly_flat_merge(ii, :), tmp2);
            end
        end
        
        % iterate using other layers if any are available
        if any(isnan(ind_layer_flat))
            
            % determine which layers overlap with original polyfit layers and order polyfit iteration based on the length of their overlap
            tmp1            = zeros(1, length(ind_layer_flat));
            for ii = find(isnan(ind_layer_flat))'
                tmp1(ii)    = length(find(~isnan(depth_curr(ii, tmp3))));
            end
            [tmp1, tmp2]    = sort(tmp1, 'descend');
            tmp5            = tmp2(logical(tmp1));
            set(status_box, 'string', ['Iterating flattening for ' num2str(length(tmp2(logical(tmp1)))) ' overlapping layers out of ' num2str(length(find(isnan(ind_layer_flat)))) ' not fitted initially...'])
            
            pause(0.1)
            
            kk              = 0;
            
            for ii = tmp5
                
                kk                  = kk + 1;
                set(status_box, 'string', ['Adding layer #' num2str(ii) ' (' num2str(kk) ' / ' num2str(length(tmp5)) ') to flattening...'])
                pause(0.1)
                
                % calculate flattening matrix based on current polynomials
                tmp1                = find(~isnan(depth_curr(ii, :)));
                depth_mat           = single(depth(:, ones(1, length(tmp1)))); % depth matrix
                depth_flat          = ((depth_mat .^ 3) .* pk.poly_flat_merge(ones((num_sample - max(ind_corr)), 1), tmp1)) + ...
                                      ((depth_mat .^ 2) .* pk.poly_flat_merge((2 .* ones((num_sample - max(ind_corr)), 1)), tmp1)) + ...
                                      (depth_mat .* (pk.poly_flat_merge((3 .* ones((num_sample - max(ind_corr)), 1)), tmp1))) + ...
                                      pk.poly_flat_merge((4 .* ones((num_sample - max(ind_corr)), 1)), tmp1);
                depth_mat           = 0;
                depth_flat(depth_flat < 0) ...
                                    = 0; % limit too-low depths
                depth_flat(depth_flat > depth(end)) ...
                                    = depth(end); % limit too-high depths
                tmp3                = find(sum(~isnan(depth_flat)));
                
                % find best-fit value at original trace
                tmp4                = NaN(1, length(tmp1));
                tmp3                = tmp3(~isnan(depth_curr(ii, tmp1(tmp3)))); % indices where both flattening and new layer exist
                
                if isempty(tmp3)
                    set(status_box, 'string', ['Layer #' num2str(ii) ' has no overlap...'])
                    pause(0.1)
                    continue
                end
                
                % depth of current layer at reference trace for all overlapping traces
                for jj = tmp3
                    [~, tmp2]       = unique(depth_flat(:, jj), 'last');
                    tmp2            = intersect((1 + find(diff(depth_flat(:, jj)) > 0)), tmp2);
                    tmp4(jj)        = interp1(depth_flat(tmp2, jj), depth(tmp2), depth_curr(ii, tmp1(jj)), 'linear', NaN);
                end
                
                ind_layer_flat(ii)  = nanmean(tmp4); % best guess depth at reference trace 
                
                % extract best layers again, now including the new layer
                tmp4                = tmp1(tmp3);
                tmp1                = depth_curr(~isnan(ind_layer_flat), tmp4);
                tmp2                = ind_layer_flat(~isnan(ind_layer_flat));
                tmp3                = find(sum(~isnan(tmp1)) > 3);
                tmp4                = tmp4(sum(~isnan(tmp1)) > 3);
                
                % new polynomials using additional "depth" for this layer
                if parallel_check
                    tmp1            = tmp1(:, tmp3);
                    tmp3            = pk.poly_flat_merge(:, tmp3);
                    parfor jj = 1:length(tmp4)
                        tmp3(:, jj) = polyfit(tmp2(~isnan(tmp1(:, jj))), tmp1(~isnan(tmp1(:, jj)), jj), 3)'; %#ok<PFBNS>
                    end
                    pk.poly_flat_merge(:, tmp4) ...
                                    = tmp3;
                else
                    for jj = 1:length(tmp3)
                        pk.poly_flat_merge(:, tmp4(jj)) ...
                                    = polyfit(tmp2(~isnan(tmp1(:, tmp3(jj)))), tmp1(~isnan(tmp1(:, tmp3(jj))), tmp3(jj)), 3)';
                    end
                end
                
                % smooth polynomials again
                tmp2                = round(mean(pk.length_smooth) / nanmean(diff(pk.dist(ind_decim))));
                if parallel_check
                    tmp1            = pk.poly_flat_merge;
                    parfor jj = 1:4
                        tmp1(jj, :) = smooth_lowess(tmp1(jj, :), tmp2);
                    end
                    pk.poly_flat_merge ...
                                    = tmp1;
                    tmp1            = 0;
                else
                    for jj = 1:4
                        pk.poly_flat_merge(jj, :) ...
                                    = smooth_lowess(pk.poly_flat_merge(jj, :), tmp2);
                    end
                end
            end
        end
        
        ind_layer_flat      = ind_layer_flat(2:end); % drop surface
        
        if parallel_check
            pctRunOnAll warning('on', 'MATLAB:polyfit:RepeatedPointsOrRescale')
            pctRunOnAll warning('on', 'MATLAB:polyfit:PolyNotUnique')
        else
            warning('on', 'MATLAB:polyfit:RepeatedPointsOrRescale')
            warning('on', 'MATLAB:polyfit:PolyNotUnique')
        end
        
        flatten_breakout
        
        set(mgui, 'keypressfcn', @keypress)
        
    end

%% Flatten amplitudes using layers

    function flatten_breakout(source, eventdata)
        
        if (logical(p_bedflat) && ishandle(p_bedflat))
            delete(p_bedflat)
        end
        if (any(p_pkflat) && any(ishandle(p_pkflat)))
            delete(p_pkflat(logical(p_pkflat) & ishandle(p_pkflat)))
        end
        if (any(p_blockflat) && any(ishandle(p_blockflat)))
            delete(p_blockflat(logical(p_blockflat) & ishandle(p_blockflat)))
        end
        if (any(p_blocknumflat) && any(ishandle(p_blocknumflat)))
            delete(p_blocknumflat(logical(p_blocknumflat) & ishandle(p_blocknumflat)))
        end
        if (any(p_coreflat) && any(ishandle(p_coreflat)))
            delete(p_coreflat(logical(p_coreflat) & ishandle(p_coreflat)))
        end
        if (any(p_corenameflat) && any(ishandle(p_corenameflat)))
            delete(p_corenameflat(logical(p_corenameflat) & ishandle(p_corenameflat)))
        end
        
        set(status_box, 'string', 'Flattening data based on polynomials...')
        
        tic
        
        % flattening matrix based on layer depth fits
        depth_mat           = single(depth(:, ones(1, num_decim))); % depth matrix
        depth_flat          = ((depth_mat .^ 3) .* pk.poly_flat_merge(ones((num_sample - max(ind_corr)), 1), :)) + ((depth_mat .^ 2) .* pk.poly_flat_merge((2 .* ones((num_sample - max(ind_corr)), 1)), :)) + ...
                              (depth_mat .* (pk.poly_flat_merge((3 .* ones((num_sample - max(ind_corr)), 1)), :))) + pk.poly_flat_merge((4 .* ones((num_sample - max(ind_corr)), 1)), :);
        depth_mat           = 0;
        if any(isnan(depth_flat(:))) % fix nans
            for ii = find(any(isnan(depth_flat)))
                if (length(find(~isnan(depth_flat(:, ii)))) > 2)
                    depth_flat(isnan(depth_flat(:, ii)), ii) ...
                            = interp1(find(~isnan(depth_flat(:, ii))), depth_flat(~isnan(depth_flat(:, ii)), ii), find(isnan(depth_flat(:, ii))), 'linear', 'extrap');
                end
            end
        end
        depth_flat(depth_flat < 0) ...
                            = 0;
        depth_flat(depth_flat > depth(end)) ...
                            = depth(end);
        set(status_box, 'string', ['Calculated remapping in ' num2str(toc, '%.0f') ' s...'])
        pause(0.1)
        
        tic
        
        % flattened radargram based on the polyfits
        tmp1                = flipud(amp_mean);
        [amp_flat, tmp2]    = deal(NaN((num_sample - max(ind_corr)), num_decim, 'single'));
        for ii = 1:num_decim
            tmp2(:, ii)     = tmp1(ind_corr(ii):(end - (max(ind_corr) - ind_corr(ii)) - 1), ii); % grab starting at surface
        end
        tmp1                = tmp2;
        tmp2                = find(sum(~isnan(depth_flat)));
        if parallel_check
            tmp1            = tmp1(:, tmp2);
            tmp3            = depth_flat(:, tmp2);
            tmp4            = amp_flat(:, tmp2);
            pctRunOnAll warning('off', 'MATLAB:interp1:NaNinY')
            parfor ii = 1:length(tmp2)
                tmp4(:, ii) = interp1(depth, tmp1(:, ii), tmp3(:, ii));
            end
            pctRunOnAll warning('on', 'MATLAB:interp1:NaNinY')
            amp_flat(:, tmp2) ...
                            = tmp4;
        else
            warning('off', 'MATLAB:interp1:NaNinY')
            for ii = 1:length(tmp2)
                amp_flat(:, tmp2(ii)) ...
                            = interp1(depth, tmp1(:, tmp2(ii)), depth_flat(:, tmp2(ii)), 'linear');
            end
            warning('on', 'MATLAB:interp1:NaNinY')
        end
        tmp1                = 0;
        
        set(status_box, 'string', ['Flattened amplitude in ' num2str(toc, '%.0f') ' s...'])
        pause(0.1)
        
        % flatten bed pick
        if (any(bed_avail) && any(~isnan(pk.elev_bed)))
            depth_bed_flat  = NaN(1, num_decim);
            for ii = 1:length(tmp2)
                [~, tmp1]   = unique(depth_flat(:, tmp2(ii)), 'last');
                tmp1        = intersect((1 + find(diff(depth_flat(:, tmp2(ii))) > 0)), tmp1);
                if (length(tmp1) > 1)
                    depth_bed_flat(tmp2(ii)) ...
                            = interp1(depth_flat(tmp1, tmp2(ii)), depth(tmp1), depth_bed(tmp2(ii)), 'nearest', 'extrap');
                end
            end
            depth_bed_flat((depth_bed_flat < 0) | (depth_bed_flat > depth(end))) ...
                            = NaN;
        end
        
        % flatten merged layer depths
        warning('off', 'MATLAB:interp1:NaNinY')
        depth_layer_flat    = NaN(pk.num_layer, num_decim);
        for ii = 1:length(tmp2)
            [~, tmp1]       = unique(depth_flat(:, tmp2(ii)), 'last');
            tmp1            = intersect((1 + find(diff(depth_flat(:, tmp2(ii))) > 0)), tmp1);
            if any(~isnan(pk.depth_smooth(:, ind_decim(tmp2(ii)))))
                depth_layer_flat(~isnan(pk.depth_smooth(:, ind_decim(tmp2(ii)))), tmp2(ii)) ...
                            = interp1(depth_flat(tmp1, tmp2(ii)), depth(tmp1), pk.depth_smooth(~isnan(pk.depth_smooth(:, ind_decim(tmp2(ii)))), ind_decim(tmp2(ii))), 'nearest', 'extrap');
            end
        end
        warning('on', 'MATLAB:interp1:NaNinY')
        
        % plot flat layers
        axes(ax_radar)
        p_pkflat            = zeros(1, pk.num_layer);
        for ii = 1:pk.num_layer
            if ~isempty(find(~isnan(depth_layer_flat(ii, :)), 1))
                p_pkflat(ii)    = plot(pk.dist_lin(ind_decim(~isnan(depth_layer_flat(ii, :)))), depth_layer_flat(ii, ~isnan(depth_layer_flat(ii, :))), '.', 'markersize', 12, 'color', colors(ii, :), 'visible', 'off');
            else
                p_pkflat(ii)    = plot(0, 0, '.', 'markersize', 12, 'color', colors(ii, :), 'visible', 'off');
            end
            if isnan(ind_layer_flat(ii))
                set(p_pkflat(ii), 'markersize', 6)
            end
        end
        
        [p_blockflat, p_blocknumflat] ...
                            = deal(zeros(1, length(pk.ind_trace_start)));
        for ii = 1:length(pk.ind_trace_start)
            p_blockflat(ii) = plot(repmat(pk.dist_lin(pk.ind_trace_start(ii)), 1, 2), [depth_min_ref depth_max_ref], 'm--', 'linewidth', 1, 'visible', 'off');
            if ~isempty(pk.file_block{ii})
                p_blocknumflat(ii) ...
                            = text((pk.dist_lin(pk.ind_trace_start(ii)) + 1), (depth_min_ref + 100), pk.file_block{ii}((end - 1):end), 'color', 'm', 'fontsize', (size_font - 2), 'visible', 'off');
            end
        end
        
        if core_done
            for ii = 1:num_int
                p_coreflat(ii) ...
                            = plot(repmat(pk.dist_lin(ind_int(ii)), 1, 2), [depth_min_ref depth_max_ref], 'm', 'linewidth', 2, 'visible', 'off');
                p_corenameflat(ii) ...
                            = text((pk.dist_lin(ind_int(ii)) + 1), (depth_min_ref + 250), name_core{int_core{curr_year}{curr_trans}(ii, 3)}, 'color', 'm', 'fontsize', size_font, 'visible', 'off');
            end
        end
        
        flat_done           = true;
        edit_flag           = false;
        set(disp_group, 'selectedobject', disp_check(2))
        disp_type           = 'flattened';
        choose_layer
        [depth_min, depth_max] ...
                            = deal(depth_min_ref, depth_max_ref);
        plot_flat
        set(status_box, 'string', 'Merged radargram flattened.')
        
    end

%% Choose the current layer

    function choose_layer(source, eventdata)
        curr_layer          = get(layer_list, 'value');
        if merge_done
            if (any(p_pk) && any(ishandle(p_pk)))
                set(p_pk(logical(p_pk) & ishandle(p_pk)), 'markersize', 12)
            end
            if (logical(p_pk(curr_layer)) && ishandle(p_pk(curr_layer)))
                set(p_pk(curr_layer), 'markersize', 24)
            end
            if flat_done
                if (any(p_pkflat) && any(ishandle(p_pkflat)))
                    set(p_pkflat(logical(p_pkflat) & ishandle(p_pkflat)), 'markersize', 12)
                end
                for ii = find(isnan(ind_layer_flat))'
                    if (logical(p_pkflat(ii)) && ishandle(p_pkflat(ii)))
                        set(p_pkflat(ii), 'markersize', 6)
                    end
                end
                if (logical(p_pkflat(curr_layer)) && ishandle(p_pkflat(curr_layer)))
                    set(p_pkflat(curr_layer), 'markersize', 24)
                end
            end
        end
    end

%% Choose a layer manually

    function choose_pk(source, eventdata)
        if ~merge_done
            set(status_box, 'string', 'No layers to focus on.')
           return
        end
        set(status_box, 'string', 'Choose a layer to highlight...')
        pause(0.1)
        [ind_x_pk, ind_y_pk]= ginput(1);
        switch disp_type
            case 'amplitude'
                [tmp1, tmp2]= unique(pk.elev_smooth(:, interp1(pk.dist_lin(ind_decim), ind_decim, ind_x_pk, 'nearest', 'extrap')));
            case 'flattened'
                [tmp1, tmp2]= unique(depth_layer_flat(:, interp1(ind_decim, 1:num_decim, interp1(pk.dist_lin(ind_decim), ind_decim, ind_x_pk, 'nearest', 'extrap'), 'nearest', 'extrap')));
        end
        tmp2                = tmp2(~isnan(tmp1));
        tmp1                = tmp1(~isnan(tmp1));
        if (length(tmp1) > 1)
            curr_layer      = interp1(tmp1, tmp2(~isnan(tmp1)), ind_y_pk, 'nearest', 'extrap');
        else
            curr_layer      = tmp2;
        end
        set(layer_list, 'value', curr_layer)
        choose_layer
        set(status_box, 'string', ['Layer #' num2str(curr_layer) ' chosen.'])
    end

%% Focus x/y axes on current layer

    function focus_layer(source, eventdata)
        if ~any(~isnan(pk.elev_smooth(curr_layer, ind_decim)))
            set(status_box, 'string', 'Cannot focus on a layer that is hidden.')
            return
        end
        axes(ax_radar)
        xlim([pk.dist_lin(find(~isnan(pk.elev_smooth(curr_layer, :)), 1)) pk.dist_lin(find(~isnan(pk.elev_smooth(curr_layer, :)), 1, 'last'))])
        tmp1                = get(ax_radar, 'xlim');
        [tmp1(1), tmp1(2)]  = deal((tmp1(1) - diff(tmp1)), (tmp1(2) + diff(tmp1)));
        if (tmp1(1) < dist_min_ref)
            tmp1(1)         = dist_min_ref;
        end
        if (tmp1(2) > dist_max_ref)
            tmp1(2)         = dist_max_ref;
        end
        xlim(tmp1)
        [dist_min, dist_max] = deal(tmp1(1), tmp1(2));
        if (dist_min < get(dist_min_slide, 'min'))
            set(dist_min_slide, 'value', get(dist_min_slide, 'min'))
        elseif (dist_min > get(dist_min_slide, 'max'))
            set(dist_min_slide, 'value', get(dist_min_slide, 'max'))
        else
            set(dist_min_slide, 'value', dist_min)
        end
        if (dist_max < get(dist_max_slide, 'min'))
            set(dist_max_slide, 'value', get(dist_max_slide, 'min'))
        elseif (dist_max > get(dist_max_slide, 'max'))
            set(dist_max_slide, 'value', get(dist_max_slide, 'max'))
        else
            set(dist_max_slide, 'value', dist_max)
        end
        set(dist_min_edit, 'string', sprintf('%3.1f', dist_min))
        set(dist_max_edit, 'string', sprintf('%3.1f', dist_max))
        switch disp_type
            case 'amplitude'
                ylim([min(pk.elev_smooth(curr_layer, ind_decim)) max(pk.elev_smooth(curr_layer, ind_decim))])
                tmp1        = get(ax_radar, 'ylim');
                [tmp1(1), tmp1(2)] ...
                            = deal((tmp1(1) - (2 * diff(tmp1))), (tmp1(2) + (2 * diff(tmp1))));
                if (tmp1(1) < elev_min_ref)
                    tmp1(1) = elev_min_ref;
                end
                if (tmp1(2) > elev_max_ref)
                    tmp1(2) = elev_max_ref;
                end
                ylim(tmp1)
                [elev_min, elev_max] ...
                            = deal(tmp1(1), tmp1(2));
                if flat_done
                    depth_min       = depth_min - (elev_max - tmp1(2));
                    depth_max       = depth_max - (elev_min - tmp1(1));
                    if (depth_min < depth_min_ref)
                        depth_min   = depth_min_ref;
                    end
                    if (depth_max > depth_max_ref)
                        depth_max   = depth_max_ref;
                    end
                end
                if (elev_min < get(z_min_slide, 'min'))
                    set(z_min_slide, 'value', get(z_min_slide, 'min'))
                elseif (elev_min > get(z_min_slide, 'max'))
                    set(z_min_slide, 'value', get(z_min_slide, 'max'))
                else
                    set(z_min_slide, 'value', elev_min)
                end
                if (elev_max < get(z_max_slide, 'min'))
                    set(z_max_slide, 'value', get(z_max_slide, 'min'))
                elseif (elev_max > get(z_max_slide, 'max'))
                    set(z_max_slide, 'value', get(z_max_slide, 'max'))
                else
                    set(z_max_slide, 'value', elev_max)
                end
                set(z_min_edit, 'string', sprintf('%4.0f', elev_min))
                set(z_max_edit, 'string', sprintf('%4.0f', elev_max))
            case 'flattened'
                ylim([min(depth_layer_flat(curr_layer, ind_decim)) max(depth_layer_flat(curr_layer, ind_decim))])
                tmp1        = get(ax_radar, 'ylim');
                [tmp1(1), tmp1(2)] ...
                            = deal((tmp1(1) - (2 * diff(tmp1))), (tmp1(2) + (2 * diff(tmp1))));
                if (tmp1(1) < depth_min_ref)
                    tmp1(1) = depth_min_ref;
                end
                if (tmp1(2) > depth_max_ref)
                    tmp1(2) = depth_max_ref;
                end
                ylim(tmp1)
                [depth_min, depth_max] ...
                            = deal(tmp1(1), tmp1(2));
                elev_min    = elev_min - (depth_max - tmp1(2));
                elev_max    = elev_max - (depth_min - tmp1(1));
                if (elev_min < elev_min_ref)
                    elev_min= elev_min_ref;
                end
                if (elev_max > elev_max_ref)
                    elev_max= elev_max_ref;
                end
                if ((depth_max_ref - (depth_max - depth_min_ref)) < get(z_min_slide, 'min'))
                    set(z_min_slide, 'value', get(z_min_slide, 'min'))
                elseif ((depth_max_ref - (depth_max - depth_min_ref)) > get(z_min_slide, 'max'))
                    set(z_min_slide, 'value', get(z_min_slide, 'max'))
                else
                    set(z_min_slide, 'value', (depth_max_ref - (depth_max - depth_min_ref)))
                end
                if ((depth_max_ref - (depth_min - depth_min_ref)) < get(z_max_slide, 'min'))
                    set(z_max_slide, 'value', get(z_max_slide, 'min'))
                elseif ((depth_max_ref - (depth_min - depth_min_ref)) > get(z_max_slide, 'max'))
                    set(z_max_slide, 'value', get(z_max_slide, 'max'))
                else
                    set(z_max_slide, 'value', (depth_max_ref - (depth_min - depth_min_ref)))
                end
                set(z_min_edit, 'string', sprintf('%4.0f', depth_max))
                set(z_max_edit, 'string', sprintf('%4.0f', depth_min))
        end
        narrow_cb
        set(status_box, 'string', ['Focused on Layer #' num2str(curr_layer) '.'])
    end

%% Delete layer

    function del_layer(source, eventdata)
        if ~merge_done
            set(status_box, 'string', 'No layers to delete yet. Load picks first.')
            return
        end
        set(status_box, 'string', 'Delete current layer? Y: yes; otherwise: no.')
        waitforbuttonpress
        if ~strcmpi(get(mgui, 'currentcharacter'), 'Y')
            set(status_box, 'string', 'Layer deletion cancelled by user.')
            return
        end
        tmp1                = setdiff(1:pk.num_layer, curr_layer);
        if (logical(p_pk(curr_layer)) && ishandle(p_pk(curr_layer)))
            delete(p_pk(curr_layer))
        end
        for ii = 1:num_var_layer
            eval(['pk.' var_layer{ii} ' = pk. ' var_layer{ii} '(tmp1, :);'])
        end
        [pk.num_layer, p_pk]= deal((pk.num_layer - 1), p_pk(tmp1));
        layer_str           = num2cell(1:pk.num_layer);
        for ii = 1:pk.num_layer
            if ~any(~isnan(pk.elev_smooth(ii, ind_decim)))
                layer_str{ii} ...
                    = [num2str(layer_str{ii}) ' H'];
            end
        end
        if flat_done
            if (logical(p_pkflat(curr_layer)) && ishandle(p_pkflat(curr_layer)))
                delete(p_pkflat(curr_layer))
            end
            p_pkflat        = p_pkflat(tmp1);
            depth_layer_flat= depth_layer_flat(tmp1, :);
            ind_layer_flat  = ind_layer_flat(tmp1);
            edit_flag       = true;
        end
        curr_layer          = curr_layer - 1;
        if (isempty(curr_layer) || any(curr_layer < 1) || any(curr_layer > pk.num_layer))
            curr_layer      = 1;
        end
        if pk.num_layer
            set(layer_list, 'string', layer_str, 'value', curr_layer)
        else
            set(layer_list, 'string', 'N/A', 'value', 1)
        end
        set(status_box, 'string', ['Layer #' num2str(curr_layer) ' deleted.'])
        choose_layer
    end

%% Switch to previous layer in list

    function pk_last(source, eventdata)
        if (merge_done && (curr_layer > 1))
            curr_layer      = curr_layer - 1;
            set(layer_list, 'value', curr_layer)
            if (any(p_pk) && any(ishandle(p_pk)))
                set(p_pk(logical(p_pk) & ishandle(p_pk)), 'markersize', 12)
            end
            if (logical(p_pk(curr_layer)) && ishandle(p_pk(curr_layer)))
                set(p_pk(curr_layer), 'markersize', 24)
            end
            if flat_done
                if (any(p_pkflat) && any(ishandle(p_pkflat)))
                    set(p_pkflat(logical(p_pkflat) & ishandle(p_pkflat)), 'markersize', 12)
                end
                if (logical(p_pkflat(curr_layer)) && ishandle(p_pkflat(curr_layer)))
                    set(p_pkflat(curr_layer), 'markersize', 24)
                end
                for ii = find(isnan(ind_layer_flat))'
                    if (logical(p_pkflat(ii)) && ishandle(p_pkflat(ii)))
                        set(p_pkflat(ii), 'markersize', 6)
                    end
                end
            end
        end
    end

%% Switch to next layer in the list

    function pk_next(source, eventdata)
        if (merge_done && (curr_layer < pk.num_layer))
            curr_layer      = curr_layer + 1;
            set(layer_list, 'value', curr_layer)
            if (any(p_pk) && any(ishandle(p_pk)))
                set(p_pk(logical(p_pk) & ishandle(p_pk)), 'markersize', 12)
            end
            if (logical(p_pk(curr_layer)) && ishandle(p_pk(curr_layer)))
                set(p_pk(curr_layer), 'markersize', 24)
            end
            if flat_done
                if (any(p_pkflat) && any(ishandle(p_pkflat)))
                    set(p_pkflat(logical(p_pkflat) & ishandle(p_pkflat)), 'markersize', 12)
                end
                if (logical(p_pkflat(curr_layer)) && ishandle(p_pkflat(curr_layer)))
                    set(p_pkflat(curr_layer), 'markersize', 24)
                end
                for ii = find(isnan(ind_layer_flat))'
                    if (logical(p_pkflat(ii)) && ishandle(p_pkflat(ii)))
                        set(p_pkflat(ii), 'markersize', 6)
                    end
                end
            end
        end
    end

%% Match two layers

    function pk_match(source, eventdata)
        
        set(status_box, 'string', ['Pick layer to merge with layer #' num2str(curr_layer) ' (Q: cancel)...'])
        pause(0.1)
        
        axes(ax_radar)
        
        % get pick and convert to indices
        [ind_x_pk, ind_y_pk, button] ...
                            = ginput(1);
        
        if strcmpi(char(button), 'Q')
            set(status_box, 'string', 'Layer merging cancelled by user.')
            return
        end
        
        % get current layer positions at ind_x_pk, depending on what we're working with
        switch disp_type
            case 'amplitude'
                [tmp1, tmp2]= unique(pk.elev_smooth(:, interp1(pk.dist_lin(ind_decim), ind_decim, ind_x_pk, 'nearest', 'extrap')));
            case 'flattened'
                [tmp1, tmp2]= unique(depth_layer_flat(:, interp1(ind_decim, 1:num_decim, interp1(pk.dist_lin(ind_decim), ind_decim, ind_x_pk, 'nearest', 'extrap'), 'nearest', 'extrap')));
        end
        tmp2                = tmp2(~isnan(tmp1));
        tmp1                = tmp1(~isnan(tmp1));
        if (length(tmp1) > 1)
            tmp1            = interp1(tmp1, tmp2(~isnan(tmp1)), ind_y_pk, 'nearest', 'extrap');
        else
            tmp1            = tmp2;
        end
        if (tmp1 == curr_layer)
            set(status_box, 'string', 'Aborted merging because picked layer for merging is the same as current layer.')
            return
        end
        if isempty(tmp1)
            set(status_box, 'string', 'No layer picked to be merged. Pick more precisely.')
            return
        end
        
        if (logical(p_pk(tmp1)) && ishandle(p_pk(tmp1)))
            set(p_pk(tmp1), 'color', colors(curr_layer, :))
        end
        if (flat_done && logical(p_pkflat(tmp1)) && ishandle(p_pkflat(tmp1)))
            set(p_pkflat(tmp1), 'color', colors(curr_layer, :))
        end
        pause(0.1)
        
        set(status_box, 'string', 'Merging correct? (Y: yes; otherwise: cancel)...')
        
        waitforbuttonpress
        if ~strcmpi(get(mgui, 'currentcharacter'), 'Y')
            set(p_pk(tmp1), 'color', colors(tmp1, :))
            if flat_done
                set(p_pkflat(tmp1), 'markersize', 6, 'color', colors(tmp1, :))
            end
            set(status_box, 'string', 'Layer merging cancelled by user.')
            return
        end
        
        % assign layer variables in tmp1 that aren't nan to the same layer variable in curr_layer
        for ii = 1:num_var_layer
            eval(['pk.' var_layer{ii} '(curr_layer, ~isnan(pk.' var_layer{ii} '(tmp1, :))) = pk.' var_layer{ii} '(tmp1, ~isnan(pk. ' var_layer{ii} '(tmp1, :)));']);
        end
        if flat_done
            depth_layer_flat(curr_layer, ~isnan(depth_layer_flat(tmp1, :))) ...
                            = depth_layer_flat(tmp1, ~isnan(depth_layer_flat(tmp1, :)));
        end
        
        % redo plots
        if (any(p_pk([curr_layer tmp1])) && any(ishandle(p_pk([curr_layer tmp1]))))
            delete(p_pk([curr_layer tmp1]))
            p_pk(curr_layer)= plot(pk.dist_lin(~isnan(pk.elev_smooth(curr_layer, :))), pk.elev_smooth(curr_layer, ~isnan(pk.elev_smooth(curr_layer, :))), ...
                                   '.', 'color', colors(curr_layer, :), 'markersize', 24, 'visible', 'off');
            if flat_done
                if (any(p_pkflat([curr_layer tmp1])) && any(ishandle(p_pkflat([curr_layer tmp1]))))
                    delete(p_pkflat([curr_layer tmp1]))
                end
                if any(~isnan(depth_layer_flat(curr_layer, :)))
                    p_pkflat(curr_layer)= plot(pk.dist_lin(ind_decim(~isnan(depth_layer_flat(curr_layer, :)))), depth_layer_flat(curr_layer, ~isnan(depth_layer_flat(curr_layer, :))), '.', ...
                                           'color', colors(curr_layer, :), 'markersize', 24, 'visible', 'off');
                else
                    p_pkflat(curr_layer)= plot(0, 0, '.', 'markersize', 12, 'color', colors(ii, :), 'visible', 'off');
                end
                if isnan(ind_layer_flat(curr_layer))
                    set(p_pkflat(curr_layer), 'markersize', 6)
                end
                edit_flag   = true;
            end
        end
        for ii = 1:num_var_layer
            eval(['pk.' var_layer{ii} ' = pk.' var_layer{ii} '(setdiff(1:pk.num_layer, tmp1), :);']) % remove old layer
        end
        [p_pk, colors]      = deal(p_pk(setdiff(1:pk.num_layer, tmp1)), colors(setdiff(1:pk.num_layer, tmp1), :));
        if flat_done
            [p_pkflat, depth_layer_flat, ind_layer_flat] ...
                            = deal(p_pkflat(setdiff(1:pk.num_layer, tmp1)), depth_layer_flat(setdiff(1:pk.num_layer, tmp1), :), ind_layer_flat(setdiff(1:pk.num_layer, tmp1))');
        end
        if isrow(ind_layer_flat)
            ind_layer_flat  = ind_layer_flat';
        end
        pk.num_layer        = pk.num_layer - 1;
        show_pk
        
        layer_str           = num2cell(1:pk.num_layer);
        for ii = 1:pk.num_layer
            if ~any(~isnan(pk.elev_smooth(ii, ind_decim)))
                layer_str{ii} ...
                    = [num2str(layer_str{ii}) ' H'];
            end
        end
        set(layer_list, 'string', layer_str, 'value', 1)
        set(status_box, 'string', ['Layers #' num2str(curr_layer) ' and #' num2str(tmp1) ' merged.'])
        if (curr_layer > tmp1)
            curr_layer      = curr_layer - 1;
        end
        if (isempty(curr_layer) || any(curr_layer < 1) || any(curr_layer > pk.num_layer))
            curr_layer      = 1;
        end
        set(layer_list, 'value', curr_layer)
        choose_layer
        
    end

%% Update original pick files with merged picks

    function update_pk(source, eventdata)
        
        if ~merge_done
            set(status_box, 'string', 'Layers not merged yet, or need to be re-merged.')
            return
        end
        
        set(status_box, 'string', 'Update original pick files with merged picks? Y: yes; otherwise: no...')
        waitforbuttonpress
        if ~strcmpi(get(mgui, 'currentcharacter'), 'Y')
            set(status_box, 'string', 'Cancelled updating of original pick files.')
            return
        end
        
        if merge_file
            [tmp1, tmp2]    = deal(file_pk, path_pk);
            set(status_box, 'string', 'Locate original picks...')
            if ~isempty(path_pk)
                [file_pk, path_pk] = uigetfile('*.mat', 'Locate original picks...', path_pk, 'multiselect', 'on');
            elseif ~isempty(path_data)
                [file_pk, path_pk] = uigetfile('*.mat', 'Locate original picks...', path_data, 'multiselect', 'on');
            else
                [file_pk, path_pk] = uigetfile('*.mat', 'Locate original picks...', 'multiselect', 'on');
            end
            if isnumeric(file_pk)
                [file_pk, path_pk] = deal('', tmp2);
            end
            if ~isempty(file_pk)
                if ischar(file_pk)
                    file_pk = {file_pk}; % only one pick file
                end
                num_pk      = length(file_pk);
                if (num_pk ~= length(pk.file_pk))
                    set(status_box, 'string', 'Number of pick files does not match data from merged picks.')
                    return
                end
            else
                file_pk     = tmp1;
                set(status_box, 'string', 'Original picks not located. Merging cancelled.')
                return
            end
        end
        
        tmp1                = pk;
        
        try
            for ii = 1:num_pk
                set(status_box, 'string', ['Updating ' file_pk{ii} ' (' num2str(ii) ' / ' num2str(num_pk) ')...'])
                pause(0.1)
                pk          = load([path_pk file_pk{ii}]);
                pk          = pk.pk;
                if isfield(pk, 'ind_layer_keep')
                    pk.ind_keep_phase ...
                            = pk.ind_layer_keep;
                    pk      = rmfield(pk, 'ind_layer_keep');
                end
                if isfield(pk, 'ind_layer_keep_phase')
                    pk.ind_keep_phase ...
                            = pk.ind_layer_keep_phase;
                    pk      = rmfield(pk, 'ind_layer_keep_phase');
                end
                if isfield(pk, 'ind_layer_keep_aresp')
                    pk.ind_keep_aresp ...
                            = pk.ind_layer_keep_aresp;
                    pk      = rmfield(pk, 'ind_layer_keep_aresp');
                end
                if isfield(pk, 'num_layer_keep')
                    pk.num_keep_phase ...
                            = pk.num_layer_keep;
                    pk      = rmfield(pk, 'num_layer_keep');
                end
                if isfield(pk, 'num_layer_keep_aresp')
                    pk.num_keep_aresp ...
                            = pk.num_layer_keep_aresp;
                    pk      = rmfield(pk, 'num_layer_keep_aresp');
                end
                if isfield(pk, 'ind_x_start')
                    pk.ind_x_start_phase ...
                            = pk.ind_x_start;
                    pk      = rmfield(pk, 'ind_x_start');
                end
                if isfield(pk, 'ind_y_start')
                    pk.ind_y_start_phase ...
                            = pk.ind_y_start;
                    pk      = rmfield(pk, 'ind_y_start');
                end
                pk.layer    = struct;
                tmp3        = tmp1.ind_trace_start(ii):(tmp1.ind_trace_start(ii) + pk.num_trace - 1); % indices of current pick file in the merged picks
                pk.ind_match= find(sum(~isnan(tmp1.elev(:, tmp3)), 2)); % layers that exist in merged picks within range of current pick file
                pk.num_layer= length(pk.ind_match); % number of layers within range of current pick file
                for jj = 1:pk.num_layer
                    for kk = 1:num_var_layer
                        eval(['pk.layer(jj).' var_layer{kk} ' = tmp1.' var_layer{kk} '(pk.ind_match(jj), tmp3);'])
                    end
                    pk.layer(jj).type ...
                            = 'max';
                end
                pk          = orderfields(pk);
                save([path_pk file_pk{ii}], '-v7.3', 'pk')
            end
        catch
            pk              = tmp1;
            set(status_box, 'string', ['Something crashed. Updating cancelled. Got through pick file #' num2str(ii - 1) '.'])
            return
        end
        
        pk                  = tmp1;
        
        set(status_box, 'string', 'Original pick files updated.')
    end

%% Save merged picks

    function pk_save(source, eventdata)

        if ~merge_done
            set(status_box, 'string', 'Layers not merged yet, or need to be re-merged.')
            return
        end
        if ~get(date_check, 'value')
            set(status_box, 'string', 'Extracted date has not been confirmed (check date box).')
            return
        end
        if edit_flag
            set(status_box, 'string', 'Layers merged or deleted since flattening. Re-flatten prior to saving.')
            return
        end
        
        reset_xz
        pause(0.1)
        
        % didn't get to flatten :(
        if ~flat_done
            set(status_box, 'string', 'Merged layers not flattened...continue saving? Y: yes; otherwise: no.')
            waitforbuttonpress
            if ~strcmpi(get(mgui, 'currentcharacter'), 'Y')
                set(status_box, 'string', 'Saving cancelled by user.')
                return
            end
            pk.poly_flat_merge ...
                            = [];
        end
        
        pk.merge_flag       = true; % will help when merged picks are loaded later
        pk.ind_x_ref        = ind_decim(ind_x_ref);
        pk                  = orderfields(pk);
        
        if isempty(path_save)
            if ispc
                if (exist([path_pk '..\..\merge'], 'dir') || exist([path_pk '..\..\..\merge'], 'dir'))
                    path_save = [path_pk(1:strfind(path_pk, '\pk')) 'merge\']; 
                end
            else
                if (exist([path_pk '../../merge'], 'dir') || exist([path_pk '../../../merge'], 'dir'))
                    path_save = [path_pk(1:strfind(path_pk, '/pk')) 'merge/']; 
                end
            end
        end

        if (~isempty(path_save) && ~isempty(file_save))
            [file_save, path_save] = uiputfile('*.mat', 'Save picks:', [path_save file_save]);
        elseif ~isempty(path_save)
            [file_save, path_save] = uiputfile('*.mat', 'Save picks:', [path_save file_pk_short '_pk_merge.mat']);
        elseif ~isempty(path_pk)
            [file_save, path_save] = uiputfile('*.mat', 'Save picks:', [path_pk file_pk_short '_pk_merge.mat']);
        elseif ~isempty(path_data)
            [file_save, path_save] = uiputfile('*.mat', 'Save picks:', [path_data file_pk_short '_pk_merge.mat']);
        else
            [file_save, path_save] = uiputfile('*.mat', 'Save picks:', [file_pk_short '_pk_merge.mat']);
        end
        
        if ~ischar(file_save)
            [file_save, path_save] ...
                            = deal('');
            set(status_box, 'string', 'Saving cancelled.')
        else
            
            set(status_box, 'string', 'Saving merged picks...')
            pause(0.1)
            
            save([path_save file_save], '-v7.3', 'pk')
            
            % make a simple figure that also gets saved
            set(0, 'DefaultFigureWindowStyle', 'default')
            pkfig           = figure('position', [10 10 1600 1000]);
            imagesc(pk.dist_lin(ind_decim), elev, amp_mean, [db_min db_max])
            hold on
            axis xy tight
            colormap(bone)
            for ii = 1:pk.num_layer
                plot(pk.dist_lin(ind_decim), pk.elev_smooth(ii, ind_decim), '.', 'color', colors(ii, :), 'markersize', 12)
            end
            for ii = 1:num_int
                plot(pk.dist_lin(ind_int(ii)), get(gca, 'ylim'), 'm', 'linewidth', 2)
                text((pk.dist_lin(ind_int(ii)) + 1), (elev_max_ref - 250), name_core{int_core{curr_year}{curr_trans}(ii, 3)}, 'color', 'm', 'fontsize', size_font)
            end
            set(gca, 'fontsize', 18)
            xlabel('Distance (km)')
            ylabel('Elevation (m)')
            title(file_pk_short, 'fontweight', 'bold', 'interpreter', 'none')
            grid on
            box on
            print(pkfig, '-dpng', [path_save file_save(1:(end - 4)) '.png'])
            
            set(status_box, 'string', ['Merged picks saved as ' file_save(1:(end - 4)) ' in ' path_save '.'])
            
        end
    end

%% Load core intersection data

    function load_core(source, eventdata)
        
        if core_done
            set(status_box, 'string', 'Core intersections already loaded.')
            return
        end
        
        if merge_file
            if ispc
                if (~isempty(path_pk) && exist([path_pk '..\..\mat'], 'dir'))
                    path_core = [path_pk '..\..\mat\'];
                end
            else
                if (~isempty(path_pk) && exist([path_pk '../../mat/'], 'dir'))
                    path_core = [path_pk '../../mat/'];
                end
            end
        else
            if ispc
                if (~isempty(path_pk) && exist([path_pk '..\..\..\mat\'], 'dir'))
                    path_core = [path_pk '..\..\..\mat\'];
                elseif (~isempty(path_pk) && exist([path_pk '..\..\..\..\mat\'], 'dir'))
                    path_core = [path_pk '..\..\..\..\mat\'];
                end
            else
                if (~isempty(path_pk) && exist([path_pk '../../../mat/'], 'dir'))
                    path_core = [path_pk '../../../mat/'];
                elseif (~isempty(path_pk) && exist([path_pk '../../../../mat/'], 'dir'))
                    path_core = [path_pk '../../../../mat/'];
                end
            end
        end
        
        if (~isempty(path_core) && exist([path_core 'core_int.mat'], 'file'))
            file_core       = 'core_int.mat';
        else
            % Dialog box to choose picks file to load
            if ~isempty(path_core)
                [file_core, path_core] = uigetfile('*.mat', 'Load core intersections (core_int.mat):', path_core);
            elseif ~isempty(path_pk)
                [file_core, path_core] = uigetfile('*.mat', 'Load core intersections (core_int.mat):', path_pk);
            elseif ~isempty(path_data)
                [file_core, path_core] = uigetfile('*.mat', 'Load core intersections (core_int.mat):', path_data);
            else
                [file_core, path_core] = uigetfile('*.mat', 'Load core intersections (core_int.mat):');
            end
            if isnumeric(file_core)
                [file_core, path_core] = deal('');
            end
        end
        
        if ~isempty(file_core)
            
            set(status_box, 'string', 'Loading core intersections...')
            pause(0.1)
            
            % load core intersection file
            tmp1        = load([path_core file_core]);
            try
                [int_core, name_core, rad_threshold, name_trans, num_year, num_trans] ...
                        = deal(tmp1.int_core, tmp1.name_core, tmp1.rad_threshold, tmp1.name_trans, tmp1.num_year, tmp1.num_trans);
            catch % give up, force restart
                set(status_box, 'string', [file_core ' does not contain the expected variables. Try again.'])
                return
            end
            
            load_core_breakout
            
        else
            set(status_box, 'string', 'No core intersections loaded.')
        end
    end

%% Load core breakout

    function load_core_breakout(source, eventdata)
        
        if merge_done
            
            if (any(p_core) && any(ishandle(p_core)))
                delete(p_core(logical(p_core) & ishandle(p_core)))
            end
            if (any(p_corename) && any(ishandle(p_corename)))
                delete(p_corename(logical(p_corename) & ishandle(p_corename)))
            end
            if (any(p_coreflat) && any(ishandle(p_coreflat)))
                delete(p_coreflat(logical(p_coreflat) & ishandle(p_coreflat)))
            end
            if (any(p_corenameflat) && any(ishandle(p_corenameflat)))
                delete(p_corenameflat(logical(p_corenameflat) & ishandle(p_corenameflat)))
            end
            
            for ii = 1:num_year
                for jj = 1:num_trans(ii)
                    if strcmp(file_pk_short, name_trans{ii}{jj})
                        break
                    end
                end
                if strcmp(file_pk_short, name_trans{ii}{jj})
                    break
                end
            end
            
            curr_year       = ii;
            curr_trans      = jj;
            
            if ~isempty(int_core{curr_year}{curr_trans})
                
                num_int     = size(int_core{curr_year}{curr_trans}, 1);
                ind_int     = zeros(1, num_int);
                for ii = 1:num_int
                    [~, ind_int(ii)] ...
                            = min(sqrt(((pk.x - int_core{curr_year}{curr_trans}(ii, 4)) .^ 2) + ((pk.y - int_core{curr_year}{curr_trans}(ii, 5)) .^ 2)));
                end
                
                [p_core, p_corename] ...
                            = deal(zeros(1, num_int));
                for ii = 1:num_int
                    p_core(ii) ...
                            = plot(repmat(pk.dist_lin(ind_int(ii)), 1, 2), [elev_min_ref elev_max_ref], 'm', 'linewidth', 2, 'visible', 'off');
                    p_corename(ii) ...
                            = text((pk.dist_lin(ind_int(ii)) + 1), (elev_max_ref - 250), name_core{int_core{curr_year}{curr_trans}(ii, 3)}, 'color', 'm', 'fontsize', size_font, 'visible', 'off');
                end
                
                if (data_done && flat_done)
                    for ii = 1:num_int
                        p_coreflat(ii) ...
                            = plot(repmat(pk.dist_lin(ind_int(ii)), 1, 2), [depth_min_ref depth_max_ref], 'm', 'linewidth', 2, 'visible', 'off');
                        p_corenameflat(ii) ...
                            = text((pk.dist_lin(ind_int(ii)) + 1), (depth_min_ref + 250), name_core{int_core{curr_year}{curr_trans}(ii, 3)}, 'color', 'm', 'fontsize', size_font, 'visible', 'off');
                    end
                end
                
                set(status_box, 'string', ['Core intersections loaded. ' num2str(num_int) ' for this transect within ' num2str(rad_threshold) ' km.'])
                
            else
                set(status_box, 'string', 'Core intersections loaded but none for this transect.')
            end
            
        else
            set(status_box, 'string', 'Core intersections loaded.')
        end
        
        core_done       = true;
        set(core_check, 'value', 1)
        show_core
        
    end

%% Update minimum elevation/depth value

    function slide_z_min(source, eventdata)
        switch disp_type
            case 'amplitude'
                if (get(z_min_slide, 'value') < elev_max)
                    if flat_done
                        tmp2        = [elev_min elev_max];
                    end
                    if get(zfix_check, 'value')
                        tmp1        = elev_max - elev_min;
                    end
                    elev_min        = get(z_min_slide, 'value');
                    if get(zfix_check, 'value')
                        elev_max    = elev_min + tmp1;
                        if (elev_max > elev_max_ref)
                            elev_max = elev_max_ref;
                            elev_min = elev_max - tmp1;
                            if (elev_min < elev_min_ref)
                                elev_min = elev_min_ref;
                            end
                            if (elev_min < get(z_min_slide, 'min'))
                                set(z_min_slide, 'value', get(z_min_slide, 'min'))
                            else
                                set(z_min_slide, 'value', elev_min)
                            end
                        end
                        set(z_max_edit, 'string', sprintf('%4.0f', elev_max))
                        if (elev_max > get(z_max_slide, 'max'))
                            set(z_max_slide, 'value', get(z_max_slide, 'max'))
                        else
                            set(z_max_slide, 'value', elev_max)
                        end
                    end
                    if flat_done
                        if get(zfix_check, 'value')
                            depth_min       = depth_min - (elev_max - tmp2(2));
                            depth_max       = depth_max - (elev_min - tmp2(1));
                            if (depth_min < depth_min_ref)
                                depth_min   = depth_min_ref;
                            end
                            if (depth_max > depth_max_ref)
                                depth_max   = depth_max_ref;
                            end
                        elseif ((depth_max - (elev_min - tmp2(1))) > depth_min)
                            depth_max       = depth_max - (elev_min - tmp2(1));
                        end
                    end
                    set(z_min_edit, 'string', sprintf('%4.0f', elev_min))
                    update_z_range
                else
                    if (elev_min < get(z_min_slide, 'min'))
                        set(z_min_slide, 'value', get(z_min_slide, 'min'))
                    else
                        set(z_min_slide, 'value', elev_min)
                    end
                end
            case 'flattened'
                if ((depth_max_ref - (get(z_min_slide, 'value') - depth_min_ref)) > depth_min)
                    tmp2            = [depth_min depth_max];
                    if get(zfix_check, 'value')
                        tmp1        = depth_max - depth_min;
                    end
                    depth_max       = depth_max_ref - (get(z_min_slide, 'value') - depth_min_ref);
                    if get(zfix_check, 'value')
                        depth_min    = depth_max - tmp1;
                        if (depth_min < depth_min_ref)
                            depth_min = depth_min_ref;
                            depth_max = depth_min + tmp1;
                            if (depth_max > depth_max_ref)
                                depth_max = depth_max_ref;
                            end
                            if ((depth_max_ref - (depth_max - depth_min_ref)) < get(z_min_slide, 'min'))
                                set(z_min_slide, 'value', get(z_min_slide, 'min'))
                            elseif ((depth_max_ref - (depth_max - depth_min_ref)) > get(z_min_slide, 'max'))
                                set(z_min_slide, 'value', get(z_min_slide, 'max'))
                            else
                                set(z_min_slide, 'value', (depth_max_ref - (depth_max - depth_min_ref)))
                            end
                        end
                        set(z_min_edit, 'string', sprintf('%4.0f', depth_max))
                        if ((depth_max_ref - (depth_min - depth_min_ref)) < get(z_max_slide, 'min'))
                            set(z_max_slide, 'value', get(z_max_slide, 'min'))
                        elseif ((depth_max_ref - (depth_min - depth_min_ref)) > get(z_max_slide, 'max'))
                            set(z_max_slide, 'value', get(z_max_slide, 'max'))
                        else
                            set(z_max_slide, 'value', (depth_max_ref - (depth_min - depth_min_ref)))
                        end
                    end
                    if get(zfix_check, 'value')
                        elev_min        = elev_min - (depth_max - tmp2(2));
                        elev_max        = elev_max - (depth_min - tmp2(1));
                        if (elev_min < elev_min_ref)
                            elev_min    = elev_min_ref;
                        end
                        if (elev_max > elev_max_ref)
                            elev_max    = elev_max_ref;
                        end
                    elseif ((elev_max - (depth_min - tmp2(1))) > elev_min)
                        elev_max        = elev_max - (depth_min - tmp2(1));
                    end
                    set(z_min_edit, 'string', sprintf('%4.0f', depth_max))
                    update_z_range
                else
                    if ((depth_max_ref - (depth_max - depth_min_ref)) < get(z_min_slide, 'min'))
                        set(z_min_slide, 'value', get(z_min_slide, 'min'))
                    elseif ((depth_max_ref - (depth_max - depth_min_ref)) > get(z_min_slide, 'max'))
                        set(z_min_slide, 'value', get(z_min_slide, 'max'))
                    else
                        set(z_min_slide, 'value', (depth_max_ref - (depth_max - depth_min_ref)))
                    end
                end
        end
        set(z_min_slide, 'enable', 'off')
        drawnow
        set(z_min_slide, 'enable', 'on')
    end

%% Update maximum elevation/depth value

    function slide_z_max(source, eventdata)
        switch disp_type
            case 'amplitude'
                if (get(z_max_slide, 'value') > elev_min)
                    if flat_done
                        tmp2        = [elev_min elev_max];
                    end
                    if get(zfix_check, 'value')
                        tmp1        = elev_max - elev_min;
                    end
                    elev_max        = get(z_max_slide, 'value');
                    if get(zfix_check, 'value')
                        elev_min    = elev_max - tmp1;
                        if (elev_min < elev_min_ref)
                            elev_min = elev_min_ref;
                            elev_max = elev_min + tmp1;
                            if (elev_max > elev_max_ref)
                                elev_max = elev_max_ref;
                            end
                            if (elev_max > get(z_max_slide, 'max'))
                                set(z_max_slide, 'value', get(z_max_slide, 'max'))
                            else
                                set(z_max_slide, 'value', elev_max)
                            end
                        end
                        set(z_min_edit, 'string', sprintf('%4.0f', elev_min))
                        if (elev_min < get(z_min_slide, 'min'))
                            set(z_min_slide, 'value', get(z_min_slide, 'min'))
                        else
                            set(z_min_slide, 'value', elev_min)
                        end
                    end
                    if flat_done
                        if get(zfix_check, 'value')
                            depth_min   = depth_min - (elev_max - tmp2(2));
                            depth_max   = depth_max - (elev_min - tmp2(1));
                            if (depth_min < depth_min_ref)
                                depth_min = depth_min_ref;
                            end
                            if (depth_max > depth_max_ref)
                                depth_max = depth_max_ref;
                            end
                        elseif ((depth_min - (elev_max - tmp2(2))) < depth_max)
                            depth_min   = depth_min - (elev_max - tmp2(2));
                        end
                    end
                    set(z_max_edit, 'string', sprintf('%4.0f', elev_max))
                    update_z_range
                else
                    if (elev_max > get(z_max_slide, 'max'))
                        set(z_max_slide, 'value', get(z_max_slide, 'max'))
                    else
                        set(z_max_slide, 'value', elev_max)
                    end
                end
            case 'flattened'
                if ((depth_max_ref - (get(z_max_slide, 'value') - depth_min_ref)) < depth_max)
                    tmp2            = [depth_min depth_max];
                    if get(zfix_check, 'value')
                        tmp1        = depth_max - depth_min;
                    end
                    depth_min       = depth_max_ref - (get(z_max_slide, 'value') - depth_min_ref);
                    if get(zfix_check, 'value')
                        depth_max   = depth_min + tmp1;
                        if (depth_max > depth_max_ref)
                            depth_max = depth_max_ref;
                            depth_min = depth_max - tmp1;
                            if (depth_min < depth_min_ref)
                                depth_min = depth_min_ref;
                            end
                            if ((depth_max_ref - (depth_min - depth_min_ref)) < get(z_max_slide, 'min'))
                                set(z_max_slide, 'value', get(z_max_slide, 'min'))
                            elseif ((depth_max_ref - (depth_min - depth_min_ref)) > get(z_max_slide, 'max'))
                                set(z_max_slide, 'value', get(z_max_slide, 'max'))
                            else
                                set(z_max_slide, 'value', (depth_max_ref - (depth_min - depth_min_ref)))
                            end
                        end
                        set(z_min_edit, 'string', sprintf('%4.0f', depth_max))
                        if ((depth_max_ref - (depth_max - depth_min_ref)) < get(z_min_slide, 'min'))
                            set(z_min_slide, 'value', get(z_min_slide, 'min'))
                        elseif ((depth_max_ref - (depth_max - depth_min_ref)) > get(z_min_slide, 'max'))
                            set(z_min_slide, 'value', get(z_min_slide, 'max'))
                        else
                            set(z_min_slide, 'value', (depth_max_ref - (depth_max - depth_min_ref)))
                        end
                    end
                    if get(zfix_check, 'value')
                        elev_min        = elev_min - (depth_max - tmp2(2));
                        elev_max        = elev_max - (depth_min - tmp2(1));
                        if (elev_min < elev_min_ref)
                            elev_min    = elev_min_ref;
                        end
                        if (elev_max > elev_max_ref)
                            elev_max    = elev_max_ref;
                        end
                    elseif ((elev_min - (depth_max - tmp2(2))) < elev_max)
                        elev_min        = elev_min - (depth_max - tmp2(2));
                    end
                    set(z_max_edit, 'string', sprintf('%4.0f', depth_min))
                    update_z_range
                else
                    if ((depth_max_ref - (depth_min - depth_min_ref)) < get(z_max_slide, 'min'))
                        set(z_max_slide, 'value', get(z_max_slide, 'min'))
                    elseif ((depth_max_ref - (depth_min - depth_min_ref)) > get(z_max_slide, 'max'))
                        set(z_max_slide, 'value', get(z_max_slide, 'max'))
                    else
                        set(z_max_slide, 'value', (depth_max_ref - (depth_min - depth_min_ref)))
                    end
                end
        end
        set(z_max_slide, 'enable', 'off')
        drawnow
        set(z_max_slide, 'enable', 'on')
    end

%% Reset minimum elevation/depth value

    function reset_z_min(source, eventdata)
        elev_min            = elev_min_ref;
        depth_max           = depth_max_ref;
        switch disp_type
            case 'amplitude'
                if (elev_min_ref < get(z_min_slide, 'min'))
                    set(z_min_slide, 'value', get(z_min_slide, 'min'))
                else
                    set(z_min_slide, 'value', elev_min_ref)
                end
                set(z_min_edit, 'string', sprintf('%4.0f', elev_min_ref))
            case 'flattened'
                if (depth_min_ref < get(z_min_slide, 'min'))
                    set(z_min_slide, 'value', get(z_min_slide, 'min'))
                elseif (depth_min_ref > get(z_min_slide, 'max'))
                    set(z_min_slide, 'value', get(z_min_slide, 'max'))
                else
                    set(z_min_slide, 'value', depth_min_ref)
                end
                set(z_min_edit, 'string', sprintf('%4.0f', depth_max_ref))
        end
        update_z_range
    end

%% Reset maximum elevation/depth value

    function reset_z_max(source, eventdata)
        elev_max            = elev_max_ref;
        depth_min           = depth_min_ref;
        switch disp_type
            case 'amplitude'
                if (elev_max_ref > get(z_max_slide, 'max'))
                    set(z_max_slide, 'value', get(z_max_slide, 'max'))
                else
                    set(z_max_slide, 'value', elev_max_ref)
                end
                set(z_max_edit, 'string', sprintf('%4.0f', elev_max_ref))
            case 'flattened'
                if (depth_max_ref < get(z_max_slide, 'min'))
                    set(z_max_slide, 'value', get(z_max_slide, 'min'))
                elseif (depth_max_ref > get(z_max_slide, 'max'))
                    set(z_max_slide, 'value', get(z_max_slide, 'max'))
                else
                    set(z_max_slide, 'value', depth_max_ref)
                end
                set(z_max_edit, 'string', sprintf('%4.0f', depth_min_ref))
        end
        update_z_range
    end

%% Update elevation/depth-axis range

    function update_z_range(source, eventdata)
        axes(ax_radar)
        switch disp_type
            case 'amplitude'
                ylim([elev_min elev_max])
            case 'flattened'
                ylim([depth_min depth_max])
        end
        narrow_cb
    end

%% Update minimum dB

    function slide_db_min(source, eventdata)
        if (get(cb_min_slide, 'value') < db_max)
            if get(cbfix_check1, 'value')
                tmp1 = db_max - db_min;
            end
            db_min          = get(cb_min_slide, 'value');
            if get(cbfix_check1, 'value')
                db_max = db_min + tmp1;
                if (db_max > db_max_ref)
                    db_max  = db_max_ref;
                    db_min  = db_max - tmp1;
                    if (db_min < db_min_ref)
                        db_min = db_min_ref;
                    end
                    if (db_min < get(cb_min_slide, 'min'))
                        set(cb_min_slide, 'value', get(cb_min_slide, 'min'))
                    else
                        set(cb_min_slide, 'value', db_min)
                    end
                end
                set(cb_max_edit, 'string', sprintf('%3.0f', db_max))
                if (db_max > get(cb_max_slide, 'max'))
                    set(cb_max_slide, 'value', get(cb_max_slide, 'max'))
                else
                    set(cb_max_slide, 'value', db_max)
                end
            end
            set(cb_min_edit, 'string', sprintf('%3.0f', db_min))
            update_db_range
        else
            if (db_min < get(cb_min_slide, 'min'))
                set(cb_min_slide, 'value', get(cb_min_slide, 'min'))
            else
                set(cb_min_slide, 'value', db_min)
            end
        end
        set(cb_min_slide, 'enable', 'off')
        drawnow
        set(cb_min_slide, 'enable', 'on')
    end

%% Update maximum dB

    function slide_db_max(source, eventdata)
        if (get(cb_max_slide, 'value') > db_min)
            if get(cbfix_check1, 'value')
                tmp1        = db_max - db_min;
            end
            db_max = get(cb_max_slide, 'value');
            if get(cbfix_check1, 'value')
                db_min = db_max - tmp1;
                if (db_min < db_min_ref)
                    db_min  = db_min_ref;
                    db_max  = db_min + tmp1;
                    if (db_max > db_max_ref)
                        db_max = db_max_ref;
                    end
                    if (db_max > get(cb_max_slide, 'max'))
                        set(cb_max_slide, 'value', get(cb_max_slide, 'max'))
                    else
                        set(cb_max_slide, 'value', db_max)
                    end
                end
                set(cb_min_edit, 'string', sprintf('%3.0f', db_min))
                if (db_min < get(cb_min_slide, 'min'))
                    set(cb_min_slide, 'value', get(cb_min_slide, 'min'))
                else
                    set(cb_min_slide, 'value', db_min)
                end
            end
            set(cb_max_edit, 'string', sprintf('%3.0f', db_max))
            update_db_range
        else
            if (db_max > get(cb_max_slide, 'max'))
                set(cb_max_slide, 'value', get(cb_max_slide, 'max'))
            else
                set(cb_max_slide, 'value', db_max)
            end
        end
        set(cb_max_slide, 'enable', 'off')
        drawnow
        set(cb_max_slide, 'enable', 'on')
    end

%% Reset minimum dB

    function reset_db_min(source, eventdata)
        if (db_min_ref < get(cb_min_slide, 'min'))
            set(cb_min_slide, 'value', get(cb_min_slide, 'min'))
        else
            set(cb_min_slide, 'value', db_min_ref)
        end
        set(cb_min_edit, 'string', num2str(db_min_ref))
        db_min      = db_min_ref;
        update_db_range
    end

%% Reset maximum dB

    function reset_db_max(source, eventdata)
        if (db_max_ref > get(cb_max_slide, 'max'))
            set(cb_max_slide, 'value', get(cb_max_slide, 'max'))
        else
            set(cb_max_slide, 'value', db_max_ref)
        end
        set(cb_max_edit, 'string', num2str(db_max_ref))
        db_max      = db_max_ref;
        update_db_range
    end

%% Update dB range

    function update_db_range(source, eventdata)
        axes(ax_radar)
        caxis([db_min db_max])
    end

%% Update minimum distance

    function slide_dist_min(source, eventdata)
        if (get(dist_min_slide, 'value') < dist_max)
            if get(distfix_check, 'value')
                tmp1        = dist_max - dist_min;
            end
            dist_min        = get(dist_min_slide, 'value');
            if get(distfix_check, 'value')
                dist_max    = dist_min + tmp1;
                if (dist_max > dist_max_ref)
                    dist_max= dist_max_ref;
                    dist_min= dist_max - tmp1;
                    if (dist_min < dist_min_ref)
                        dist_min ...
                            = dist_min_ref;
                    end
                    if (dist_min < get(dist_min_slide, 'min'))
                        set(dist_min_slide, 'value', get(dist_min_slide, 'min'))
                    else
                        set(dist_min_slide, 'value', dist_min)
                    end
                end
                if (dist_max > get(dist_max_slide, 'max'))
                    set(dist_max_slide, 'value', get(dist_max_slide, 'max'))
                else
                    set(dist_max_slide, 'value', dist_max)
                end
                set(dist_max_slide, 'value', dist_max)
                set(dist_max_edit, 'string', sprintf('%3.1f', dist_max))
            end
            set(dist_min_edit, 'string', sprintf('%3.1f', dist_min))
            update_dist_range
        else
            if (dist_min < get(dist_min_slide, 'min'))
                set(dist_min_slide, 'value', get(dist_min_slide, 'min'))
            else
                set(dist_min_slide, 'value', dist_min)
            end
        end
        set(dist_min_slide, 'enable', 'off')
        drawnow
        set(dist_min_slide, 'enable', 'on')
    end

%% Update maximum distance

    function slide_dist_max(source, eventdata)
        if (get(dist_max_slide, 'value') > dist_min)
            if get(distfix_check, 'value')
                tmp1        = dist_max - dist_min;
            end
            dist_max        = get(dist_max_slide, 'value');
            if get(distfix_check, 'value')
                dist_min    = dist_max - tmp1;
                if (dist_min < dist_min_ref)
                    dist_min= dist_min_ref;
                    dist_max= dist_min + tmp1;
                    if (dist_max > dist_max_ref)
                        dist_max ...
                            = dist_max_ref;
                    end
                    if (dist_max > get(dist_max_slide, 'max'))
                        set(dist_max_slide, 'value', get(dist_max_slide, 'max'))
                    else
                        set(dist_max_slide, 'value', dist_max)
                    end
                end
                if (dist_min < get(dist_min_slide, 'min'))
                    set(dist_min_slide, 'value', get(dist_min_slide, 'min'))
                else
                    set(dist_min_slide, 'value', dist_min)
                end
                set(dist_min_edit, 'string', sprintf('%3.1f', dist_min))
            end
            set(dist_max_edit, 'string', sprintf('%3.1f', dist_max))
            update_dist_range
        else
            if (dist_max > get(dist_max_slide, 'max'))
                set(dist_max_slide, 'value', get(dist_max_slide, 'max'))
            else
                set(dist_max_slide, 'value', dist_max)
            end
        end
        set(dist_max_slide, 'enable', 'off')
        drawnow
        set(dist_max_slide, 'enable', 'on')
    end

%% Reset minimum distance

    function reset_dist_min(source, eventdata)
        if (dist_min_ref < get(dist_min_slide, 'min'))
            set(dist_min_slide, 'value', get(dist_min_slide, 'min'))
        else
            set(dist_min_slide, 'value', dist_min_ref)
        end
        set(dist_min_edit, 'string', sprintf('%3.1f', dist_min_ref))
        dist_min            = dist_min_ref;
        update_dist_range
    end

%% Reset maximum distance

    function reset_dist_max(source, eventdata)
        if (dist_max_ref > get(dist_max_slide, 'max'))
            set(dist_max_slide, 'value', get(dist_max_slide, 'max'))
        else
            set(dist_max_slide, 'value', dist_max_ref)
        end
        set(dist_max_edit, 'string', sprintf('%3.1f', dist_max_ref))
        dist_max            = dist_max_ref;
        update_dist_range
    end

%% Update distance range

    function update_dist_range(source, eventdata)
        axes(ax_radar)
        xlim([dist_min dist_max])
        narrow_cb
    end

%% Reset both distance (x) and elevation (y)

    function reset_xz(source, eventdata)
        reset_z_min
        reset_z_max
        reset_dist_min
        reset_dist_max
    end

%% Adjust slider limits after panning or zooming

    function panzoom(source, eventdata)
        tmp1                = get(ax_radar, 'xlim');
        if (tmp1(1) < dist_min_ref)
            reset_dist_min
        else
            if (tmp1(1) < get(dist_min_slide, 'min'))
                set(dist_min_slide, 'value', get(dist_min_slide, 'min'))
            else
                set(dist_min_slide, 'value', tmp1(1))
            end
            set(dist_min_edit, 'string', sprintf('%3.1f', tmp1(1)))
            dist_min        = tmp1(1);
        end
        if (tmp1(2) > dist_max_ref)
            reset_dist_max
        else
            if (tmp1(2) > get(dist_max_slide, 'max'))
                set(dist_max_slide, 'value', get(dist_max_slide, 'max'))
            else
                set(dist_max_slide, 'value', tmp1(2))
            end
            set(dist_max_edit, 'string', sprintf('%3.1f', tmp1(2)))
            dist_max        = tmp1(2);
        end
        tmp1                = get(ax_radar, 'ylim');
        switch disp_type
            case 'amplitude'
                if flat_done
                    tmp2    = [elev_min elev_max];
                end
                if (tmp1(1) < elev_min_ref)
                    reset_z_min
                else
                    if (tmp1(1) < get(z_min_slide, 'min'))
                        set(z_min_slide, 'value', get(z_min_slide, 'min'))
                    else
                        set(z_min_slide, 'value', tmp1(1))
                    end
                    set(z_min_edit, 'string', sprintf('%4.0f', tmp1(1)))
                    elev_min= tmp1(1);
                end
                if (tmp1(2) > elev_max_ref)
                    reset_z_max
                else
                    if (tmp1(2) > get(z_max_slide, 'max'))
                        set(z_max_slide, 'value', get(z_max_slide, 'max'))
                    else
                        set(z_max_slide, 'value', tmp1(2))
                    end
                    set(z_max_edit, 'string', sprintf('%4.0f', tmp1(2)))
                    elev_max= tmp1(2);
                end
                if flat_done
                    depth_min       = depth_min - (elev_max - tmp2(2));
                    depth_max       = depth_max - (elev_min - tmp2(1));
                    if (depth_min < depth_min_ref)
                        depth_min   = depth_min_ref;
                    end
                    if (depth_max > depth_max_ref)
                        depth_max   = depth_max_ref;
                    end
                end
            case 'flattened'
                tmp2        = [depth_min depth_max];
                if (tmp1(1) < depth_min_ref)
                    reset_z_min
                else
                    if (tmp1(1) < get(z_min_slide, 'min'))
                        set(z_min_slide, 'value', get(z_min_slide, 'min'))
                    elseif (tmp1(1) > get(z_min_slide, 'max'))
                        set(z_min_slide, 'value', get(z_min_slide, 'max'))
                    else
                        set(z_min_slide, 'value', tmp1(1))
                    end
                    set(z_min_edit, 'string', sprintf('%4.0f', tmp1(1)))
                    depth_min ...
                            = tmp1(1);
                end
                if (tmp1(2) > depth_max_ref)
                    reset_z_max
                else
                    if (tmp1(2) < get(z_max_slide, 'min'))
                        set(z_max_slide, 'value', get(z_max_slide, 'min'))
                    elseif (tmp1(2) > get(z_max_slide, 'max'))
                        set(z_max_slide, 'value', get(z_max_slide, 'max'))
                    else
                        set(z_max_slide, 'value', tmp1(2))
                    end
                    set(z_max_edit, 'string', sprintf('%4.0f', tmp1(2)))
                    depth_max ...
                            = tmp1(2);
                end
                elev_min    = elev_min - (depth_max - tmp2(2));
                elev_max    = elev_max - (depth_min - tmp2(1));
                if (elev_min < elev_min_ref)
                    elev_min= elev_min_ref;
                end
                if (elev_max > elev_max_ref)
                    elev_max= elev_max_ref;
                end
        end
        narrow_cb
    end

%% Plot data in dB (amplitude)

    function plot_db(source, eventdata)
        if ~data_done
            set(status_box, 'string', 'Data not loaded yet.')
            return
        end
        if (logical(p_data) && ishandle(p_data))
            delete(p_data)
        end
        if (logical(p_bedflat) && ishandle(p_bedflat))
            set(p_bedflat, 'visible', 'off')
        end
        if (logical(p_refflat) && ishandle(p_refflat))
            set(p_refflat, 'visible', 'off')
        end
        if (logical(p_surf) && ishandle(p_surf))
            delete(p_surf)
        end
        if (logical(p_bed) && ishandle(p_bed))
            delete(p_bed)
        end
        axes(ax_radar)
        axis xy
        set(z_min_slide, 'min', elev_min_ref, 'max', elev_max_ref, 'value', elev_min)
        set(z_max_slide, 'min', elev_min_ref, 'max', elev_max_ref, 'value', elev_max)
        set(z_min_edit, 'string', sprintf('%4.0f', elev_min))
        set(z_max_edit, 'string', sprintf('%4.0f', elev_max))
        ylim([elev_min elev_max])
        if data_done
            p_data          = imagesc(pk.dist_lin(ind_decim), elev, amp_mean, [db_min db_max]);
        end
        if (any(surf_avail) && any(~isnan(pk.elev_surf)))
            p_surf          = plot(pk.dist_lin(ind_decim), pk.elev_surf(ind_decim), 'g--', 'linewidth', 2);
            if any(isnan(pk.elev_surf(ind_decim)))
                set(p_surf, 'marker', '.', 'linestyle', 'none', 'markersize', 12)
            end
        end
        if (any(bed_avail) && any(~isnan(pk.elev_bed)))
            p_bed           = plot(pk.dist_lin(ind_decim), pk.elev_bed(ind_decim), 'g--', 'linewidth', 2);
            if any(isnan(depth_bed))
                set(p_bed, 'marker', '.', 'linestyle', 'none', 'markersize', 12)
            end
        end
        if any(surf_avail)
            set(yl, 'string', 'Elevation (m)')
        else
            set(yl, 'string', 'Depth (m)')
        end
        narrow_cb
        show_data
        show_core
        show_block
        show_pk
    end

%% Plot layer-flattened radargram

    function plot_flat(source, eventdata)
        if ~flat_done
            set(disp_group, 'selectedobject', disp_check(1))
            disp_type       = 'amplitude';
            plot_db
            return
        end
        if (logical(p_data) && ishandle(p_data)) % get rid of old plotted data
            delete(p_data)
        end
        if (logical(p_bedflat) && ishandle(p_bedflat))
            delete(p_bedflat)
        end
        if (logical(p_refflat) && ishandle(p_refflat))
            delete(p_refflat)
        end
        if (logical(p_bed) && ishandle(p_bed))
            set(p_bed, 'visible', 'off')
        end
        if (logical(p_surf) && ishandle(p_surf))
            set(p_surf, 'visible', 'off')
        end
        axes(ax_radar)
        if (get(cmap_list, 'value') ~= 1)
            set(cmap_list, 'value', 1)
            change_cmap
        end
        axis ij
        set(z_min_slide, 'min', depth_min_ref, 'max', depth_max_ref, 'value', (depth_max_ref - (depth_max - depth_min_ref)))
        set(z_max_slide, 'min', depth_min_ref, 'max', depth_max_ref, 'value', (depth_max_ref - (depth_min - depth_min_ref)))
        set(z_min_edit, 'string', sprintf('%4.0f', depth_max))
        set(z_max_edit, 'string', sprintf('%4.0f', depth_min))
        ylim([depth_min depth_max])
        p_data              = imagesc(pk.dist_lin(ind_decim), depth, amp_flat, [db_min db_max]);
        if (any(bed_avail) && any(~isnan(depth_bed_flat)) && any(depth_bed_flat))
            if any(isnan(depth_bed_flat))
                p_bedflat   = plot(pk.dist_lin(ind_decim(~isnan(depth_bed_flat))), depth_bed_flat(~isnan(depth_bed_flat)), 'g.', 'markersize', 12);
            else
                p_bedflat   = plot(pk.dist_lin(ind_decim), depth_bed_flat, 'g--', 'linewidth', 2);
            end
        end
        p_refflat           = plot(repmat(pk.dist_lin(ind_decim(ind_x_ref)), 1, 2), [depth_min_ref depth_max_ref], 'w--', 'linewidth', 2);
        disp_type           = 'flattened';
        set(yl, 'string', 'Depth (m)')
        narrow_cb
        show_core
        show_block
        show_pk
    end

%% Show radar data

    function show_data(source, eventdata)
        if data_done
            if (get(data_check, 'value') && logical(p_data) && ishandle(p_data))
                set(p_data, 'visible', 'on');
            elseif (logical(p_data) && ishandle(p_data))
                set(p_data, 'visible', 'off')
            end
        elseif get(data_check, 'value')
            set(data_check, 'value', 0)
        end
    end

%% Show picked layers (smooth only)

    function show_pk(source, eventdata)
        if merge_done
            if get(pk_check, 'value')
                switch disp_type
                    case 'amplitude'
                        if (any(p_pk) && any(ishandle(p_pk)))
                            set(p_pk(logical(p_pk) & ishandle(p_pk)), 'visible', 'on')
                            uistack(p_pk(logical(p_pk) & ishandle(p_pk)), 'top')
                        end
                        if (any(p_pkflat) && any(ishandle(p_pkflat)))
                            set(p_pkflat(logical(p_pkflat) & ishandle(p_pkflat)), 'visible', 'off')
                        end
                    case 'flattened'
                        if (any(p_pkflat) && any(ishandle(p_pkflat)))
                            set(p_pkflat(logical(p_pkflat) & ishandle(p_pkflat)), 'visible', 'on')
                            uistack(p_pkflat(logical(p_pkflat) & ishandle(p_pkflat)), 'top')
                        end
                        if (any(p_pk) && any(ishandle(p_pk)))
                            set(p_pk(logical(p_pk) & ishandle(p_pk)), 'visible', 'off')
                        end
                end
            else
                if (any(p_pk) && any(ishandle(p_pk)))
                    set(p_pk(logical(p_pk) & ishandle(p_pk)), 'visible', 'off')
                end
                if (any(p_pkflat) && any(ishandle(p_pkflat)))
                    set(p_pkflat(logical(p_pkflat) & ishandle(p_pkflat)), 'visible', 'off')
                end
            end
        elseif get(pk_check, 'value')
            set(pk_check, 'value', 0)
        end
    end

%% Show block divisions

    function show_block(source, eventdata)
        if merge_done
            if get(block_check, 'value')
                switch disp_type
                    case 'amplitude'
                        if (any(p_block) && any(ishandle(p_block)))
                            set(p_block(logical(p_block) & ishandle(p_block)), 'visible', 'on')
                            uistack(p_block(logical(p_block) & ishandle(p_block)), 'top')
                        end
                        if (any(p_blocknum) && any(ishandle(p_blocknum)))
                            set(p_blocknum(logical(p_blocknum) & ishandle(p_blocknum)), 'visible', 'on')
                            uistack(p_blocknum(logical(p_blocknum) & ishandle(p_blocknum)), 'top')
                        end
                        if (any(p_blockflat) && any(ishandle(p_blockflat)))
                            set(p_blockflat(logical(p_blockflat) & ishandle(p_blockflat)), 'visible', 'off')
                        end
                        if (any(p_blocknumflat) && any(ishandle(p_blocknumflat)))
                            set(p_blocknumflat(logical(p_blocknumflat) & ishandle(p_blocknumflat)), 'visible', 'off')
                        end
                    case 'flattened'
                        if (any(p_blockflat) && any(ishandle(p_blockflat)))
                            set(p_blockflat(logical(p_blockflat) & ishandle(p_blockflat)), 'visible', 'on')
                            uistack(p_blockflat(logical(p_blockflat) & ishandle(p_blockflat)), 'top')
                        end
                        if (any(p_blocknumflat) && any(ishandle(p_blocknumflat)))
                            set(p_blocknumflat(logical(p_blocknumflat) & ishandle(p_blocknumflat)), 'visible', 'on')
                            uistack(p_blocknumflat(logical(p_blocknumflat) & ishandle(p_blocknumflat)), 'top')
                        end
                        if (any(p_block) && any(ishandle(p_block)))
                            set(p_block(logical(p_block) & ishandle(p_block)), 'visible', 'off')
                        end
                        if (any(p_blocknum) && any(ishandle(p_blocknum)))
                            set(p_blocknum(logical(p_blocknum) & ishandle(p_blocknum)), 'visible', 'off')
                        end
                end
            else
                if (any(p_block) && any(ishandle(p_block)))
                    set(p_block(logical(p_block) & ishandle(p_block)), 'visible', 'off')
                end
                if (any(p_blocknum) && any(ishandle(p_blocknum)))
                    set(p_blocknum(logical(p_blocknum) & ishandle(p_blocknum)), 'visible', 'off')
                end
                if (any(p_blockflat) && any(ishandle(p_blockflat)))
                    set(p_blockflat(logical(p_blockflat) & ishandle(p_blockflat)), 'visible', 'off')
                end
                if (any(p_blocknumflat) && any(ishandle(p_blocknumflat)))
                    set(p_blocknumflat(logical(p_blocknumflat) & ishandle(p_blocknumflat)), 'visible', 'off')
                end
            end
        elseif get(block_check, 'value')
            set(block_check, 'value', 0)
        end
    end

%% Show core intersections

    function show_core(source, eventdata)
        if core_done
            if get(core_check, 'value')
                switch disp_type
                    case 'amplitude'
                        if (any(p_core) && any(ishandle(p_core)))
                            set(p_core(logical(p_core) & ishandle(p_core)), 'visible', 'on')
                            uistack(p_core(logical(p_core) & ishandle(p_core)), 'top')
                        end
                        if (any(p_corename) && any(ishandle(p_corename)))
                            set(p_corename(logical(p_corename) & ishandle(p_corename)), 'visible', 'on')
                            uistack(p_corename(logical(p_corename) & ishandle(p_corename)), 'top')
                        end
                        if (any(p_coreflat) && any(ishandle(p_coreflat)))
                            set(p_coreflat(logical(p_coreflat) & ishandle(p_coreflat)), 'visible', 'off')
                        end
                        if (any(p_corenameflat) && any(ishandle(p_corenameflat)))
                            set(p_corenameflat(logical(p_corenameflat) & ishandle(p_corenameflat)), 'visible', 'off')
                        end
                    case 'flattened'
                        if (any(p_coreflat) && any(ishandle(p_coreflat)))
                            set(p_coreflat(logical(p_coreflat) & ishandle(p_coreflat)), 'visible', 'on')
                            uistack(p_coreflat(logical(p_coreflat) & ishandle(p_coreflat)), 'top')
                        end
                        if (any(p_corenameflat) && any(ishandle(p_corenameflat)))
                            set(p_corenameflat(logical(p_corenameflat) & ishandle(p_corenameflat)), 'visible', 'on')
                            uistack(p_corenameflat(logical(p_corenameflat) & ishandle(p_corenameflat)), 'top')
                        end
                        if (any(p_core) && any(ishandle(p_core)))
                            set(p_core(logical(p_core) & ishandle(p_core)), 'visible', 'off')
                        end
                        if (any(p_corename) && any(ishandle(p_corename)))
                            set(p_corename(logical(p_corename) & ishandle(p_corename)), 'visible', 'off')
                        end
                end
            else
                if (any(p_core) && any(ishandle(p_core)))
                    set(p_core(logical(p_core) & ishandle(p_core)), 'visible', 'off')
                end
                if (any(p_corename) && any(ishandle(p_corename)))
                    set(p_corename(logical(p_corename) & ishandle(p_corename)), 'visible', 'off')
                end
                if (any(p_coreflat) && any(ishandle(p_coreflat)))
                    set(p_coreflat(logical(p_coreflat) & ishandle(p_coreflat)), 'visible', 'off')
                end
                if (any(p_corenameflat) && any(ishandle(p_corenameflat)))
                    set(p_corenameflat(logical(p_corenameflat) & ishandle(p_corenameflat)), 'visible', 'off')
                end
            end
        elseif get(core_check, 'value')
            set(core_check, 'value', 0)
        end
    end

%% Switch to a chunk

    function plot_chunk(source, eventdata)
        curr_chunk          = get(chunk_list, 'value');
        axes(ax_radar) %#ok<*MAXES>
        if (curr_chunk <= num_chunk)
            xlim(dist_chunk([curr_chunk (curr_chunk + 1)]))
        else
            xlim(pk.dist_lin([1 end]))
        end
        tmp1                = get(ax_radar, 'xlim');
        [dist_min, dist_max] = deal(tmp1(1), tmp1(2));
        set(dist_min_slide, 'value', dist_min)
        set(dist_max_slide, 'value', dist_max)
        set(dist_min_edit, 'string', sprintf('%3.1f', dist_min))
        set(dist_max_edit, 'string', sprintf('%3.1f', dist_max))
    end

%% Adjust date

    function adj_date(source, eventdata)
        tmp1                = pk.date_str;
        pk.date_str         = get(date_box, 'string');
        try
            pk.date_num     = datenum(pk.date_str);
            set(status_box, 'string', ['Transect date set to ' pk.date_str '.'])
        catch
            pk.date_str     = tmp1;
            set(date_box, 'string', pk.date_str)
            set(status_box, 'string', 'Transect date must be in dd-mmm-yyyy format, e.g., 12-Oct-1984.')
        end
    end

%% Adjust number of indices to display

    function adj_decim(source, eventdata)
        decim               = abs(round(str2double(get(decim_edit, 'string'))));
        if merge_done
            if (decim > 1)
                ind_decim   = (1 + ceil(decim / 2)):decim:(pk.num_trace_tot - ceil(decim / 2));
            else
                ind_decim   = 1:pk.num_trace_tot;
            end
            num_decim       = length(ind_decim);
            if (any(p_pk) && any(ishandle(p_pk)))
                delete(p_pk(logical(p_pk) & ishandle(p_pk)))
            end
            layer_str       = num2cell(1:pk.num_layer);
            for ii = 1:pk.num_layer
                if ~any(~isnan(pk.elev_smooth(ii, ind_decim)))
                    p_pk(ii)=plot(0, 0, 'w.', 'markersize', 1);
                    layer_str{ii} ...
                            = [num2str(layer_str{ii}) ' H'];
                else
                    p_pk(ii)= plot(pk.dist_lin(ind_decim(~isnan(pk.elev_smooth(ii, ind_decim)))), pk.elev_smooth(ii, ind_decim(~isnan(pk.elev_smooth(ii, ind_decim)))), '.', 'color', colors(ii, :), ...
                                   'markersize', 12, 'visible', 'off');
                end
            end
            if curr_layer
                set(layer_list, 'string', layer_str, 'value', curr_layer)
            else
                set(layer_list, 'string', layer_str, 'value', 1)
            end
            show_pk
        end
        if data_done
            load_data_breakout
        elseif flat_done
            if (num_decim ~= size(pk.poly_flat_merge, 2))
                pk.poly_flat_merge ...
                            = interp2(pk.poly_flat_merge, linspace(1, num_decim, num_decim), (1:size(pk.poly_flat_merge, 1))');
            end
        end
        set(decim_edit, 'string', num2str(decim))
        set(status_box, 'string', ['Displayed data decimated to 1/' num2str(decim) ' samples.'])
    end

%% Adjust length of each chunk

    function adj_length_chunk(source, eventdata)
        length_chunk        = abs(round(str2double(get(length_chunk_edit, 'string'))));
        % make the horizontal chunks of data to pick (too big is unwieldy)
        num_chunk           = floor((pk.dist_lin(end) - pk.dist_lin(1)) ./ length_chunk);
        dist_chunk          = pk.dist_lin(1):length_chunk:pk.dist_lin(end);
        if (dist_chunk(end) ~= pk.dist_lin(end))
            dist_chunk(end) = pk.dist_lin(end);
        end
        set(chunk_list, 'string', [num2cell(1:num_chunk) 'full'], 'value', (num_chunk + 1))
        set(status_box, 'string', ['Display chunk length adjusted to ' num2str(length_chunk) ' km.'])
        set(length_chunk_edit, 'string', num2str(length_chunk))
    end

%% Change colormap

    function change_cmap(source, eventdata)
        axes(ax_radar)
        colormap(cmaps{get(cmap_list, 'value')})
    end

%% Switch display type

    function disp_radio(~, eventdata)
        disp_type           = get(eventdata.NewValue, 'string');
        switch disp_type
            case 'amplitude'
                plot_db
            case 'flattened'
                plot_flat
        end
    end

%% Toggle gridlines

    function toggle_grid(source, eventdata)
        if get(grid_check, 'value')
            axes(ax_radar)
            grid on
        else
            grid off
        end
    end

%% Narrow color axis to +/- 2 standard deviations of current mean value

    function narrow_cb(source, eventdata)
        if (get(cbfix_check2, 'value') && data_done)
            axes(ax_radar)
            tmp4            = zeros(2);
            tmp4(1, :)      = interp1(pk.dist_lin(ind_decim), 1:num_decim, [dist_min dist_max], 'nearest', 'extrap');
            switch disp_type
                case 'amplitude'
                    tmp4(2, :) = interp1(elev, 1:num_sample, [elev_min elev_max], 'nearest', 'extrap');
                    tmp4(2, :) = flipud(tmp4(2, :));
                    tmp4    = amp_mean(tmp4(2, 1):tmp4(2, 2), tmp4(1, 1):tmp4(1, 2));
                case 'flattened'
                    tmp4(2, :) = interp1(depth, 1:(num_sample - max(ind_corr)), [depth_min depth_max], 'nearest', 'extrap');
                    tmp4    = amp_flat(tmp4(2, 1):tmp4(2, 2), tmp4(1, 1):tmp4(1, 2));
            end
            [tmp5(1), tmp5(2)] ...
                            = deal(nanmean(tmp4(:)), nanstd(tmp4(:)));
            if any(isnan(tmp5))
                return
            end
            tmp4            = zeros(1, 2);
            if ((tmp5(1) - (2 * tmp5(2))) < db_min_ref)
                tmp4(1)     = db_min_ref;
            else
                tmp4(1)     = tmp5(1) - (2 * tmp5(2));
            end
            if ((tmp5(1) + (2 * tmp5(2))) > db_max_ref)
                tmp4(2)     = db_max_ref;
            else
                tmp4(2)     = tmp5(1) + (2 * tmp5(2));
            end
            [db_min, db_max]= deal(tmp4(1), tmp4(2));
            if (db_min < get(cb_min_slide, 'min'))
                set(cb_min_slide, 'value', get(cb_min_slide, 'min'))
            else
                set(cb_min_slide, 'value', db_min)
            end
            if (db_max > get(cb_max_slide, 'max'))
                set(cb_max_slide, 'value', get(cb_max_slide, 'max'))
            else
                set(cb_max_slide, 'value', db_max)
            end
            set(cb_min_edit, 'string', sprintf('%3.0f', db_min))
            set(cb_max_edit, 'string', sprintf('%3.0f', db_max))
            caxis([db_min db_max])
        end
    end

%% Pop out map

    function pop_map(source, eventdata)
        if ~merge_done
            set(status_box, 'string', 'Cannot make map until picks are loaded.')
            return
        end
        if ~license('checkout', 'map_toolbox')
            set(status_box, 'string', 'Cannot make map without Mapping Toolbox.')
            return
        end
        set(0, 'DefaultFigureWindowStyle', 'default')
        tmp1                = figure;
        worldmap('greenland'); % displays empty map
        tmp2                = shaperead('landareas', 'UseGeoCoords', true, 'Selector', {@(name) strcmp(name, 'Greenland'), 'Name'}); % Greenland outline
        tmp3                = defaultm('ups'); %  initializes a map projection structure
        [tmp3.geoid, tmp3.mapparallels, tmp3.falsenorthing, tmp3.falseeasting, tmp3.origin] ...
                            = deal(almanac('earth', 'wgs84', 'meters'), 70, 0, 0, [90 -45 0]); % assigns ps_struct to variables
        tmp4                = cell(1, 2);
        [tmp4{1}, tmp4{2}]  = mfwdtran(tmp3, tmp2.Lat, tmp2.Lon); % convert lat/lon to polar stereographic x/y
        [tmp4{1}, tmp4{2}]  = deal((1e-3 .* tmp4{1}), (1e-3 .* tmp4{2})); % m to km
        close(tmp1)
        figure('position', [100 100 600 800]);
        hold on
        plot(tmp4{1}, tmp4{2}, 'k', 'linewidth', 1)
        plot(pk.x, pk.y, 'k', 'linewidth', 2)
        plot(pk.x(1), pk.y(1), 'ko', 'markersize', 12, 'markerfacecolor', 'g')
        plot(pk.x(end), pk.y(end), 'ko', 'markersize', 12, 'markerfacecolor', 'r')
        tmp1                = (100 * ceil(pk.dist_lin(1) / 100)):100:(100 * floor(pk.dist_lin(end) / 100));
        plot(pk.x(interp1(pk.dist_lin, 1:pk.num_trace_tot, tmp1, 'nearest')), pk.y(interp1(pk.dist_lin, 1:pk.num_trace_tot, tmp1, 'nearest')), 'ko', 'markersize', 6, 'markerfacecolor', 'b')
        set(gca, 'fontsize', 20)
        xlabel('Polar stereographic X (km)')
        ylabel('Polar stereographic Y (km)')
        axis equal tight
        grid on
        box on
        set(0, 'DefaultFigureWindowStyle', 'docked')
    end

%% Pop out figure

    function pop_fig(source, eventdata)
        if (~merge_done && ~data_done)
            set(status_box, 'string', 'Picks or data must be loaded prior to popping out a figure.')
            return
        end
        set(0, 'DefaultFigureWindowStyle', 'default')
        figure('position', [10 10 1200 800]);
        axis tight
        hold on
        colormap(cmaps{get(cmap_list, 'value')})
        caxis([db_min db_max])
        switch disp_type
            case 'amplitude'
                axis xy
                axis([dist_min dist_max elev_min elev_max])
                imagesc(pk.dist_lin(ind_decim), elev, amp_mean, [db_min db_max])
                if (any(surf_avail) && any(~isnan(pk.elev_surf)))
                    tmp1    = plot(pk.dist_lin(ind_decim), pk.elev_surf(ind_decim), 'g--', 'linewidth', 2);
                    if any(isnan(pk.elev_surf(ind_decim)))
                        set(tmp1, 'marker', '.', 'linestyle', 'none', 'markersize', 12)
                    end
                end
                if (any(bed_avail) && any(~isnan(pk.elev_bed)))
                    tmp1    = plot(pk.dist_lin(ind_decim), pk.elev_bed(ind_decim), 'g--', 'linewidth', 2);
                    if any(isnan(depth_bed))
                        set(tmp1, 'marker', '.', 'linestyle', 'none', 'markersize', 12)
                    end
                end
                if get(block_check, 'value')
                    for ii = 1:length(pk.ind_trace_start)
                        plot(repmat(pk.dist_lin(pk.ind_trace_start(ii)), 1, 2), [elev_min elev_max], 'm--', 'linewidth', 1)
                        if ~isempty(pk.file_block{ii})
                            text((pk.dist_lin(pk.ind_trace_start(ii)) + 1), (elev_max - (0.1 * (elev_max - elev_min))), pk.file_block{ii}((end - 1):end), 'color', 'm', 'fontsize', (size_font - 2))
                        end
                    end
                end
                if get(core_check, 'value')
                    for ii = 1:num_int
                        plot(repmat(pk.dist_lin(ind_int(ii)), 1, 2), [elev_min elev_max], 'm', 'linewidth', 2)
                        text((pk.dist_lin(ind_int(ii)) + 1), (elev_max - (0.2 * (elev_max - elev_min))), name_core{int_core{curr_year}{curr_trans}(ii, 3)}, 'color', 'm', 'fontsize', size_font)
                    end
                end
                if get(pk_check, 'value')
                    for ii = 1:pk.num_layer
                        if any(~isnan(pk.elev_smooth(ii, ind_decim)))
                            tmp1 = plot(pk.dist_lin(ind_decim(~isnan(pk.elev_smooth(ii, ind_decim)))), pk.elev_smooth(ii, ind_decim(~isnan(pk.elev_smooth(ii, ind_decim)))), '.', 'color', colors(ii, :), ...
                                        'markersize', 12);
                            if any(ii == curr_layer)
                                set(tmp1, 'markersize', 24)
                            end
                        end
                    end
                end
                ylabel('Elevation (m)', 'fontsize', 20)
                text((dist_max + (0.015 * (dist_max - dist_min))), (elev_max + (0.04 * (elev_max - elev_min))), '(dB)', 'color', 'k', 'fontsize', 20)
            case 'flattened'
                axis ij
                axis([dist_min dist_max depth_min depth_max])                
                imagesc(pk.dist_lin(ind_decim), depth, amp_flat, [db_min db_max])
                if (any(bed_avail) && any(~isnan(pk.elev_bed)))
                    if any(isnan(depth_bed_flat))
                        plot(pk.dist_lin(ind_decim(~isnan(depth_bed_flat))), depth_bed_flat(~isnan(depth_bed_flat)), 'g.', 'markersize', 12)
                    else
                        plot(pk.dist_lin(ind_decim), depth_bed_flat, 'g--', 'linewidth', 2)
                    end
                end
                if get(block_check, 'value')
                    for ii = 1:length(pk.ind_trace_start)
                        plot(repmat(pk.dist_lin(pk.ind_trace_start(ii)), 1, 2), [depth_min depth_max], 'm--', 'linewidth', 1)
                        if ~isempty(pk.file_block{ii})
                            text((pk.dist_lin(pk.ind_trace_start(ii)) + 1), (depth_min + (0.1 * (depth_max - depth_min))), pk.file_block{ii}((end - 1):end), 'color', 'm', 'fontsize', (size_font - 2))
                        end
                    end
                end
                if get(core_check, 'value')
                    for ii = 1:num_int
                        plot(repmat(pk.dist_lin(ind_int(ii)), 1, 2), [depth_min depth_max], 'm', 'linewidth', 2)
                        text((pk.dist_lin(ind_int(ii)) + 1), (depth_min + (0.2 * (depth_max - depth_min))), name_core{int_core{curr_year}{curr_trans}(ii, 3)}, 'color', 'm', 'fontsize', size_font)
                    end
                end
                if get(pk_check, 'value')
                    for ii = 1:pk.num_layer
                        if ~isempty(find(~isnan(depth_layer_flat(ii, :)), 1))
                            tmp1    = plot(pk.dist_lin(ind_decim(~isnan(depth_layer_flat(ii, :)))), depth_layer_flat(ii, ~isnan(depth_layer_flat(ii, :))), '.', 'markersize', 12, 'color', colors(ii, :));
                        end
                        if isnan(ind_layer_flat(ii))
                            set(tmp1, 'markersize', 6)
                        end
                    end
                end
                ylabel('Depth (m)', 'fontsize', 20)
                text((dist_max + (0.015 * (dist_max - dist_min))), (depth_min - (0.04 * (depth_max - depth_min))), '(dB)', 'color', 'k', 'fontsize', 20)
        end
        set(gca, 'fontsize', 20, 'layer', 'top')
        xlabel('Distance (km)')
        title(file_pk_short, 'fontweight', 'bold', 'interpreter', 'none')
        colorbar('fontsize', 20)
        if get(grid_check, 'value')
            grid on
        end
        box on
        set(0, 'DefaultFigureWindowStyle', 'docked')
    end

%% Keyboard shortcuts for various functions

    function keypress(~, eventdata)

        switch eventdata.Key
            case '1'
                if get(pk_check, 'value')
                    set(pk_check, 'value', 0)
                else
                    set(pk_check, 'value', 1)
                end
                show_pk
            case '2'
                if get(data_check, 'value')
                    set(data_check, 'value', 0)
                else
                    set(data_check, 'value', 1)
                end
                show_data
            case '3'
                if get(block_check, 'value')
                    set(block_check, 'value', 0)
                else
                    set(block_check, 'value', 1)
                end
                show_block
            case '4'
                if get(core_check, 'value')
                    set(core_check, 'value', 0)
                else
                    set(core_check, 'value', 1)
                end
                show_core
            case 'a'
                pop_map
            case 'b'
                if get(cbfix_check2, 'value')
                    set(cbfix_check2, 'value', 0)
                else
                    set(cbfix_check2, 'value', 1)
                end
                narrow_cb
            case 'c'
                if ~core_done
                    load_core
                else
                    set(status_box, 'string', 'Core intersection data already loaded.')
                end
            case 'd'
                del_layer
            case 'e'
                reset_xz
            case 'f'
                flatten
            case 'g'
                if get(grid_check, 'value')
                    set(grid_check, 'value', 0)
                else
                    set(grid_check, 'value', 1)
                end
                toggle_grid
            case 'h'
                if get(date_check, 'value')
                    set(date_check, 'value', 0)
                else
                    set(date_check, 'value', 1)
                end
            case 'l'
                load_data
            case 'm'
                pk_match
            case 'n'
                pk_next
            case 'o'
                focus_layer
            case 'p'
                load_pk
            case 'q'
                pop_fig
            case 'v'
                pk_save
            case 'w'
                if (get(cmap_list, 'value') == 1)
                    set(cmap_list, 'value', 2)
                else
                    set(cmap_list, 'value', 1)
                end
                change_cmap
            case 'x'
                if get(distfix_check, 'value')
                    set(distfix_check, 'value', 0)
                else
                    set(distfix_check, 'value', 1)
                end
            case 'y'
                if get(zfix_check, 'value')
                    set(zfix_check, 'value', 0)
                else
                    set(zfix_check, 'value', 1)
                end
            case 'z'
                choose_pk
            case 'downarrow'
                switch disp_type
                    case 'amplitude'
                        tmp1        = elev_max - elev_min;
                        tmp2        = elev_min - (0.25 * tmp1);
                        if flat_done
                            tmp3    = [elev_min elev_max];
                        end
                        if (tmp2 < elev_min_ref)
                            elev_min= elev_min_ref;
                        else
                            elev_min= tmp2;
                        end
                        elev_max    = elev_min + tmp1;
                        if (elev_max > elev_max_ref)
                            elev_max= elev_max_ref;
                        end
                        if flat_done
                            depth_min       = depth_min - (elev_max - tmp3(2));
                            depth_max       = depth_max - (elev_min - tmp3(1));
                            if (depth_min < depth_min_ref)
                                depth_min   = depth_min_ref;
                            end
                            if (depth_max > depth_max_ref)
                                depth_max   = depth_max_ref;
                            end
                        end
                        set(z_min_edit, 'string', sprintf('%4.0f', elev_min))
                        set(z_max_edit, 'string', sprintf('%4.0f', elev_max))
                        if (elev_min < get(z_min_slide, 'min'))
                            set(z_min_slide, 'value', get(z_min_slide, 'min'))
                        else
                            set(z_min_slide, 'value', elev_min)
                        end
                        if (elev_max > get(z_max_slide, 'max'))
                            set(z_max_slide, 'value', get(z_max_slide, 'max'))
                        else
                            set(z_max_slide, 'value', elev_max)
                        end
                    case 'flattened'
                        tmp1        = depth_max - depth_min;
                        tmp2        = depth_max + (0.25 * tmp1);
                        tmp3        = [depth_min depth_max];
                        if (tmp2 > depth_max_ref)
                            depth_max= depth_max_ref;
                        else
                            depth_max= tmp2;
                        end
                        depth_min    = depth_max - tmp1;
                        if (depth_min < depth_min_ref)
                            depth_min=depth_min_ref;
                        end
                        elev_min       = elev_min - (depth_max - tmp3(2));
                        elev_max       = elev_max - (depth_min - tmp3(1));
                        if (elev_min < elev_min_ref)
                            elev_min   = elev_min_ref;
                        end
                        if (elev_max > elev_max_ref)
                            elev_max   = elev_max_ref;
                        end
                        set(z_min_edit, 'string', sprintf('%4.0f', depth_max))
                        set(z_max_edit, 'string', sprintf('%4.0f', depth_min))
                        if ((depth_max_ref - (depth_max - depth_min_ref)) < get(z_min_slide, 'min'))
                            set(z_min_slide, 'value', get(z_min_slide, 'min'))
                        elseif ((depth_max_ref - (depth_max - depth_min_ref)) > get(z_min_slide, 'max'))
                            set(z_min_slide, 'value', get(z_min_slide, 'max'))
                        else
                            set(z_min_slide, 'value', (depth_max_ref - (depth_max - depth_min_ref)))
                        end
                        if ((depth_max_ref - (depth_min - depth_min_ref)) < get(z_max_slide, 'min'))
                            set(z_max_slide, 'value', get(z_max_slide, 'min'))
                        elseif ((depth_max_ref - (depth_min - depth_min_ref)) > get(z_max_slide, 'max'))
                            set(z_max_slide, 'value', get(z_max_slide, 'max'))
                        else
                            set(z_max_slide, 'value', (depth_max_ref - (depth_min - depth_min_ref)))
                        end
                end
                update_z_range
            case 'leftarrow'
                tmp1        = dist_max - dist_min;
                tmp2        = dist_min - (0.25 * tmp1);
                if (tmp2 < dist_min_ref)
                    dist_min= dist_min_ref;
                else
                    dist_min= tmp2;
                end
                dist_max    = dist_min + tmp1;
                if (dist_max > dist_max_ref)
                    dist_max= dist_max_ref;
                end
                set(dist_min_edit, 'string', sprintf('%3.1f', dist_min))
                set(dist_max_edit, 'string', sprintf('%3.1f', dist_max))
                if (dist_min < get(dist_min_slide, 'min'))
                    set(dist_min_slide, 'value', get(dist_min_slide, 'min'))
                else
                    set(dist_min_slide, 'value', dist_min)
                end
                if (dist_max > get(dist_max_slide, 'max'))
                    set(dist_max_slide, 'value', get(dist_max_slide, 'max'))
                else
                    set(dist_max_slide, 'value', dist_max)
                end
                update_dist_range
            case 'rightarrow'
                tmp1        = dist_max - dist_min;
                tmp2        = dist_max + (0.25 * tmp1);
                if (tmp2 > dist_max_ref);
                    dist_max= dist_max_ref;
                else
                    dist_max= tmp2;
                end
                dist_min    = dist_max - tmp1;
                if (dist_min < dist_min_ref)
                    dist_min= dist_min_ref;
                end
                set(dist_min_edit, 'string', sprintf('%3.1f', dist_min))
                set(dist_max_edit, 'string', sprintf('%3.1f', dist_max))
                if (dist_min < get(dist_min_slide, 'min'))
                    set(dist_min_slide, 'value', get(dist_min_slide, 'min'))
                else
                    set(dist_min_slide, 'value', dist_min)
                end
                if (dist_max > get(dist_max_slide, 'max'))
                    set(dist_max_slide, 'value', get(dist_max_slide, 'max'))
                else
                    set(dist_max_slide, 'value', dist_max)
                end
                update_dist_range
            case 'uparrow'
                switch disp_type
                    case 'amplitude'
                        tmp1        = elev_max - elev_min;
                        tmp2        = elev_max + (0.25 * tmp1);
                        if flat_done
                            tmp3    = [elev_min elev_max];
                        end
                        if (tmp2 > elev_max_ref)
                            elev_max= elev_max_ref;
                        else
                            elev_max= tmp2;
                        end
                        elev_min    = elev_max - tmp1;
                        if (elev_min < elev_min_ref)
                            elev_min= elev_min_ref;
                        end
                        if flat_done
                            depth_min       = depth_min - (elev_max - tmp3(2));
                            depth_max       = depth_max - (elev_min - tmp3(1));
                            if (depth_min < depth_min_ref)
                                depth_min   = depth_min_ref;
                            end
                            if (depth_max > depth_max_ref)
                                depth_max   = depth_max_ref;
                            end
                        end
                        set(z_min_edit, 'string', sprintf('%4.0f', elev_min))
                        set(z_max_edit, 'string', sprintf('%4.0f', elev_max))
                        if (elev_min < get(z_min_slide, 'min'))
                            set(z_min_slide, 'value', get(z_min_slide, 'min'))
                        else
                            set(z_min_slide, 'value', elev_min)
                        end
                        if (elev_max > get(z_max_slide, 'max'))
                            set(z_max_slide, 'value', get(z_max_slide, 'max'))
                        else
                            set(z_max_slide, 'value', elev_max)
                        end
                    case 'flattened'
                        tmp1        = depth_max - depth_min;
                        tmp2        = depth_min - (0.25 * tmp1);
                        tmp3        = [depth_min depth_max];
                        if (tmp2 < depth_min_ref)
                            depth_min= depth_min_ref;
                        else
                            depth_min= tmp2;
                        end
                        depth_max   = depth_min + tmp1;
                        depth_min   = depth_max - tmp1;
                        elev_min    = elev_min - (depth_max - tmp3(2));
                        elev_max    = elev_max - (depth_min - tmp3(1));
                        if (elev_min < elev_min_ref)
                            elev_min= elev_min_ref;
                        end
                        if (elev_max > elev_max_ref)
                            elev_max= elev_max_ref;
                        end
                        set(z_min_edit, 'string', sprintf('%4.0f', depth_max))
                        set(z_max_edit, 'string', sprintf('%4.0f', depth_min))
                        if ((depth_max_ref - (depth_max - depth_min_ref)) < get(z_min_slide, 'min'))
                            set(z_min_slide, 'value', get(z_min_slide, 'min'))
                        elseif ((depth_max_ref - (depth_max - depth_min_ref)) > get(z_min_slide, 'max'))
                            set(z_min_slide, 'value', get(z_min_slide, 'max'))
                        else
                            set(z_min_slide, 'value', (depth_max_ref - (depth_max - depth_min_ref)))
                        end
                        if ((depth_max_ref - (depth_min - depth_min_ref)) < get(z_max_slide, 'min'))
                            set(z_max_slide, 'value', get(z_max_slide, 'min'))
                        elseif ((depth_max_ref - (depth_min - depth_min_ref)) > get(z_max_slide, 'max'))
                            set(z_max_slide, 'value', get(z_max_slide, 'max'))
                        else
                            set(z_max_slide, 'value', (depth_max_ref - (depth_min - depth_min_ref)))
                        end
                end
                update_z_range
            case 'space'
                switch get(disp_group, 'selectedobject')
                    case disp_check(1)
                        if flat_done
                            set(disp_group, 'selectedobject', disp_check(2))
                            disp_type = 'flattened';
                            plot_flat
                        end
                    case disp_check(2)
                        set(disp_group, 'selectedobject', disp_check(1))
                        disp_type = 'amplitude';
                        plot_db
                end
        end
    end

%% Close mergegui

    function closefig(source, eventdata)
        if parallel_check
            matlabpool close
        end
        close(mgui)
    end

%% Test something

    function misctest(source, eventdata)
        
    end

%%
end