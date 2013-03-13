function fencegui
% FENCEGUI Interactive comparison of layer picks between intersecting radar transects.
%   
%   FENCEGUI loads a pair of GUIs (3D and 2D) for interactive comparison of
%   a maseter transect's radar layers with those from transects that
%   intersect it. These layers must have been traced previously across
%   individual blocks using PICKGUI and then merged using MERGEGUI. Refer
%   to manual for operation (pickgui_man.docx).
%   
%   FENCEGUI requires that the functions INTERSECTI and TOPOCORR be
%   available within the user's path. If the Parallel Computing Toolbox is
%   licensed and available, then several calculations related to data
%   flattening will be parallelized.
% 
% Joe MacGregor (UTIG)
% Last updated: 03/13/13

if ~exist('intersecti', 'file')
    error('fencegui:intersecti', 'Function INTERSECTI is not available within this user''s path.')
end
if ~exist('topocorr', 'file')
    error('fencegui:topocorr', 'Function TOPOCORR is not available within this user''s path.')
end

%% Intialize variables

% elevation/depth defaults
[elev_min_ref, depth_min, depth_min_ref] ...
                            = deal(0);
[elev_max_ref, depth_max, depth_max_ref] ...
                            = deal(1);
elev_min                    = zeros(1, 2);
elev_max                    = ones(1, 2);

% x/y defaults
[x_min_ref, x_max_ref]      = deal(0, 1);
[y_min_ref, y_max_ref]      = deal(0, 1);
[dist_min_ref, dist_max_ref]= deal(zeros(1, 2), ones(1, 2));
[x_min, x_max]              = deal(x_min_ref, x_max_ref);
[y_min, y_max]              = deal(y_min_ref, y_max_ref);
[dist_min, dist_max]        = deal(dist_min_ref, dist_max_ref);

% dB default
[db_min_ref, db_max_ref]    = deal(repmat(-130, 1, 2), zeros(1, 2));
[db_min, db_max]            = deal(repmat(-80, 1, 2), repmat(-20, 1, 2));

% some default values
speed_vacuum                = 299792458; % m/s
permitt_ice                 = 3.15;
speed_ice                   = speed_vacuum / sqrt(permitt_ice);
decim                       = [25 25; 10 10]; % decimate radargram for display
cmaps                       = {'bone' 'jet'}';
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
[core_done, int_done]       = deal(false);
[bed_avail, data_done, flat_done, pk_done, surf_avail] ...
                            = deal(false(1, 2));
[amp_mean, amp_flat, colors, depth, depth_bed, depth_bed_flat, depth_flat, depth_layer_flat, depth_mat, elev, file_data, file_pk, file_pk_short, ind_corr, ind_int_core, layer_str, path_data, path_pk, pk, p_pkflat, ...
    twtt]                   = deal(cell(1, 2));
p_int1                      = cell(2, 3);
[ind_decim, p_core, p_corename, p_coreflat, p_corenameflat, p_int2, p_pk] ...
                            = deal(cell(2));
[curr_layer, curr_trans, curr_year, dt, num_data, num_int_core, num_sample, p_bedflat, p_data] ...
                            = deal(zeros(1, 2));
[decim_edit, layer_list, num_decim, p_bed, pk_check, p_surf] ...
                            = deal(zeros(2));
[curr_az2, curr_el2, curr_ind_int, h_link1, h_link2, h_link3, h_link4, ii, ind_x_pk, ind_y_pk, int_all, int_core, jj, kk, name_core, name_trans, num_int, num_trans, num_year, rad_threshold, tmp1, tmp2, tmp3, tmp4, ...
 tmp5, tmp6, tmp7, tmp8, x_core, y_core] ...
                            = deal(0);
[curr_ax, curr_gui, curr_int, curr_rad] ...
                            = deal(1);
curr_rad_alt                = 2;
disp_type                   = {'amp.' 'amp.' 'amp.'};
curr_dim                    = '3D';
[file_core, file_ref, path_core, path_ref] ...
                            = deal('');

if license('test', 'distrib_computing_toolbox')
    if matlabpool('size')
        matlabpool close
    end
    parallel_check          = true;
    matlabpool open
else
    parallel_check          = false;
end

%% draw first GUI

set(0, 'DefaultFigureWindowStyle', 'docked')
if ispc % windows switch
    fgui(1)                 = figure('toolbar', 'figure', 'name', 'FENCEGUI 3D', 'position', [1920 940 1 1], 'menubar', 'none', 'keypressfcn', @keypress1);
    size_font               = 14;
    width_slide             = 0.01;
else
    fgui(1)                 = figure('toolbar', 'figure', 'name', 'FENCEGUI 3D', 'position', [1864 1100 1 1], 'menubar', 'none', 'keypressfcn', @keypress1);
    size_font               = 18;
    width_slide             = 0.02;
end

ax(1)                       = subplot('position', [0.08 0.10 0.89 0.81]);
hold on
axis([x_min_ref x_max_ref y_min_ref y_max_ref elev_min_ref elev_max_ref])
box off
grid on
[curr_az3, curr_el3]        = view(3);
set(gca, 'fontsize', size_font, 'layer', 'top')
xlabel('X (km)')
ylabel('Y (km)')
zlabel('Elevation (m)')

% sliders
x_min_slide                 = uicontrol(fgui(1), 'style', 'slider', 'units', 'normalized', 'position', [0.11 0.04 0.24 width_slide], 'callback', @slide_x_min, 'min', 0, 'max', 1, 'value', x_min_ref, ...
                                               'sliderstep', [0.01 0.1]);
x_max_slide                 = uicontrol(fgui(1), 'style', 'slider', 'units', 'normalized', 'position', [0.11 0.01 0.24 width_slide], 'callback', @slide_x_max, 'min', 0, 'max', 1, 'value', x_max_ref, ...
                                               'sliderstep', [0.01 0.1]);
y_min_slide                 = uicontrol(fgui(1), 'style', 'slider', 'units', 'normalized', 'position', [0.67 0.04 0.24 width_slide], 'callback', @slide_y_min, 'min', 0, 'max', 1, 'value', x_min_ref, ...
                                               'sliderstep', [0.01 0.1]);
y_max_slide                 = uicontrol(fgui(1), 'style', 'slider', 'units', 'normalized', 'position', [0.67 0.01 0.24 width_slide], 'callback', @slide_y_max, 'min', 0, 'max', 1, 'value', x_max_ref, ...
                                               'sliderstep', [0.01 0.1]);
z_min_slide(1)              = uicontrol(fgui(1), 'style', 'slider', 'units', 'normalized', 'position', [0.005 0.07 width_slide 0.32], 'callback', @slide_z_min1, 'min', 0, 'max', 1, 'value', elev_min_ref, ...
                                               'sliderstep', [0.01 0.1]);
z_max_slide(1)              = uicontrol(fgui(1), 'style', 'slider', 'units', 'normalized', 'position', [0.005 0.52 width_slide 0.32], 'callback', @slide_z_max1, 'min', 0, 'max', 1, 'value', elev_max_ref, ...
                                               'sliderstep', [0.01 0.1]);

% slider values
x_min_edit                  = annotation('textbox', [0.07 0.045 0.04 0.03], 'string', num2str(x_min_ref), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
x_max_edit                  = annotation('textbox', [0.07 0.005 0.04 0.03], 'string', num2str(x_max_ref), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
y_min_edit                  = annotation('textbox', [0.96 0.045 0.04 0.03], 'string', num2str(y_min_ref), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
y_max_edit                  = annotation('textbox', [0.96 0.005 0.04 0.03], 'string', num2str(y_max_ref), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
z_min_edit(1)               = annotation('textbox', [0.005 0.39 0.04 0.03], 'string', num2str(elev_min_ref), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
z_max_edit(1)               = annotation('textbox', [0.005 0.84 0.04 0.03], 'string', num2str(elev_max_ref), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');

% push buttons
uicontrol(fgui(1), 'style', 'pushbutton', 'string', 'Load master picks', 'units', 'normalized', 'position', [0.005 0.965 0.10 0.03], 'callback', @load_pk1, 'fontsize', size_font, 'foregroundcolor', 'b')
uicontrol(fgui(1), 'style', 'pushbutton', 'string', 'Next', 'units', 'normalized', 'position', [0.245 0.925 0.03 0.03], 'callback', @pk_next1, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(fgui(1), 'style', 'pushbutton', 'string', 'Next', 'units', 'normalized', 'position', [0.50 0.925 0.03 0.03], 'callback', @pk_next2, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(fgui(1), 'style', 'pushbutton', 'string', 'Last', 'units', 'normalized', 'position', [0.215 0.925 0.03 0.03], 'callback', @pk_last1, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(fgui(1), 'style', 'pushbutton', 'string', 'Last', 'units', 'normalized', 'position', [0.47 0.925 0.03 0.03], 'callback', @pk_last2, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(fgui(1), 'style', 'pushbutton', 'string', 'test', 'units', 'normalized', 'position', [0.865 0.965 0.03 0.03], 'callback', @misctest, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(fgui(1), 'style', 'pushbutton', 'string', 'Save', 'units', 'normalized', 'position', [0.965 0.965 0.03 0.03], 'callback', @pk_save, 'fontsize', size_font, 'foregroundcolor', 'g')
uicontrol(fgui(1), 'style', 'pushbutton', 'string', 'transects', 'units', 'normalized', 'position', [0.62 0.965 0.05 0.03], 'callback', @load_int, 'fontsize', size_font, 'foregroundcolor', 'b')
uicontrol(fgui(1), 'style', 'pushbutton', 'string', 'cores', 'units', 'normalized', 'position', [0.685 0.965 0.03 0.03], 'callback', @load_core, 'fontsize', size_font, 'foregroundcolor', 'b')
uicontrol(fgui(1), 'style', 'pushbutton', 'string', 'reset x/y/z', 'units', 'normalized', 'position', [0.94 0.925 0.055 0.03], 'callback', @reset_xyz, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(fgui(1), 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.36 0.04 0.03 0.03], 'callback', @reset_x_min, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(fgui(1), 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.36 0.005 0.03 0.03], 'callback', @reset_x_max, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(fgui(1), 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.635 0.04 0.03 0.03], 'callback', @reset_y_min, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(fgui(1), 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.635 0.005 0.03 0.03], 'callback', @reset_y_max, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(fgui(1), 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.005 0.005 0.03 0.03], 'callback', @reset_z_min1, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(fgui(1), 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.005 0.48 0.03 0.03], 'callback', @reset_z_max1, 'fontsize', size_font, 'foregroundcolor', 'r')

% fixed text annotations
a(1)                        = annotation('textbox', [0.12 0.965 0.025 0.03], 'string', 'N_{decim}', 'fontsize', size_font, 'color', 'b', 'edgecolor', 'none');
a(2)                        = annotation('textbox', [0.37 0.965 0.025 0.03], 'string', 'N_{decim}', 'fontsize', size_font, 'color', 'b', 'edgecolor', 'none');
a(3)                        = annotation('textbox', [0.195 0.965 0.03 0.03], 'string', 'Layer', 'fontsize', size_font, 'color', 'm', 'edgecolor', 'none');
a(4)                        = annotation('textbox', [0.445 0.965 0.03 0.03], 'string', 'Layer', 'fontsize', size_font, 'color', 'm', 'edgecolor', 'none');
a(5)                        = annotation('textbox', [0.90 0.925 0.03 0.03], 'string', 'Grid', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(6)                        = annotation('textbox', [0.04 0.04 0.03 0.03], 'string', 'x_{min}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(7)                        = annotation('textbox', [0.04 0.005 0.03 0.03], 'string', 'x_{max}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(8)                        = annotation('textbox', [0.93 0.04 0.03 0.03], 'string', 'y_{min}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(9)                        = annotation('textbox', [0.93 0.005 0.03 0.03], 'string', 'y_{max}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(10)                       = annotation('textbox', [0.005 0.42 0.03 0.03], 'string', 'z_{min}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(11)                       = annotation('textbox', [0.005 0.87 0.03 0.03], 'string', 'z_{max}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(12)                       = annotation('textbox', [0.005 0.89 0.03 0.03], 'string', 'fix', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(13)                       = annotation('textbox', [0.395 0.005 0.03 0.03], 'string', 'fix', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(14)                       = annotation('textbox', [0.595 0.005 0.03 0.03], 'string', 'fix', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(15)                       = annotation('textbox', [0.28 0.965 0.10 0.03], 'string', 'Intersecting picks', 'fontsize', size_font, 'color', 'b', 'edgecolor', 'none');
a(16)                       = annotation('textbox', [0.535 0.965 0.10 0.03], 'string', 'Load intersections', 'fontsize', size_font, 'color', 'b', 'edgecolor', 'none');
if ~ispc
    set(a, 'fontweight', 'bold')
end

% variable text annotations
file_box(1)                 = annotation('textbox', [0.005 0.925 0.16 0.03], 'string', '', 'color', 'k', 'fontsize', size_font, 'backgroundcolor', 'w', 'edgecolor', 'k', 'interpreter', 'none');
status_box(1)               = annotation('textbox', [0.535 0.925 0.36 0.03], 'string', '', 'color', 'k', 'fontsize', size_font, 'backgroundcolor', 'w', 'edgecolor', 'k', 'interpreter', 'none');

dim_group                   = uibuttongroup('position', [0.90 0.965 0.06 0.03], 'selectionchangefcn', @choose_dim);
uicontrol(fgui(1), 'style', 'text', 'parent', dim_group, 'units', 'normalized', 'position', [0 0.6 0.9 0.3], 'fontsize', size_font)
dim_check(1)                = uicontrol(fgui(1), 'style', 'radio', 'string', '2D', 'units', 'normalized', 'position', [0.01 0.1 0.45 0.8], 'parent', dim_group, 'fontsize', size_font, 'handlevisibility', 'off');
dim_check(2)                = uicontrol(fgui(1), 'style', 'radio', 'string', '3D', 'units', 'normalized', 'position', [0.51 0.1 0.45 0.8], 'parent', dim_group, 'fontsize', size_font, 'handlevisibility', 'off');
set(dim_group, 'selectedobject', dim_check(2))

% value boxes
decim_edit(1, 1)            = uicontrol(fgui(1), 'style', 'edit', 'string', num2str(decim(1, 1)), 'units', 'normalized', 'position', [0.16 0.965 0.03 0.03], 'fontsize', size_font, 'foregroundcolor', 'k', ...
                                               'backgroundcolor', 'w', 'callback', @adj_decim1);
decim_edit(1, 2)            = uicontrol(fgui(1), 'style', 'edit', 'string', num2str(decim(1, 2)), 'units', 'normalized', 'position', [0.405 0.965 0.03 0.03], 'fontsize', size_font, 'foregroundcolor', 'k', ...
                                               'backgroundcolor', 'w', 'callback', @adj_decim2);
% menus
layer_list(1, 1)            = uicontrol(fgui(1), 'style', 'popupmenu', 'string', 'N/A', 'value', 1, 'units', 'normalized', 'position', [0.225 0.955 0.05 0.04], 'fontsize', size_font, 'foregroundcolor', 'k', ...
                                               'callback', @choose_layer1);
layer_list(1, 2)            = uicontrol(fgui(1), 'style', 'popupmenu', 'string', 'N/A', 'value', 1, 'units', 'normalized', 'position', [0.48 0.955 0.05 0.04], 'fontsize', size_font, 'foregroundcolor', 'k', ...
                                               'callback', @choose_layer2);
int_list                    = uicontrol(fgui(1), 'style', 'popupmenu', 'string', 'N/A', 'value', 1, 'units', 'normalized', 'position', [0.28 0.915 0.16 0.04], 'fontsize', size_font, 'foregroundcolor', 'k', ...
                                               'callback', @load_pk2);

% check boxes
xfix_check                  = uicontrol(fgui(1), 'style', 'checkbox', 'units', 'normalized', 'position', [0.415 0.005 0.01 0.03], 'fontsize', size_font, 'value', 0);
yfix_check                  = uicontrol(fgui(1), 'style', 'checkbox', 'units', 'normalized', 'position', [0.615 0.005 0.01 0.03], 'fontsize', size_font, 'value', 0);
zfix_check(1)               = uicontrol(fgui(1), 'style', 'checkbox', 'units', 'normalized', 'position', [0.02 0.89 0.01 0.03], 'fontsize', size_font, 'value', 0);
grid_check(1)               = uicontrol(fgui(1), 'style', 'checkbox', 'units', 'normalized', 'position', [0.925 0.925 0.01 0.03], 'callback', @toggle_grid1, 'fontsize', size_font, 'value', 1);
pk_check(1, 1)              = uicontrol(fgui(1), 'style', 'checkbox', 'units', 'normalized', 'position', [0.105 0.965 0.01 0.03], 'callback', @show_pk1, 'fontsize', size_font, 'value', 0);
pk_check(1, 2)              = uicontrol(fgui(1), 'style', 'checkbox', 'units', 'normalized', 'position', [0.36 0.965 0.01 0.03], 'callback', @show_pk2, 'fontsize', size_font, 'value', 0);
int_check(1)                = uicontrol(fgui(1), 'style', 'checkbox', 'units', 'normalized', 'position', [0.67 0.965 0.01 0.03], 'callback', @show_int1, 'fontsize', size_font, 'value', 0);
core_check(1)               = uicontrol(fgui(1), 'style', 'checkbox', 'units', 'normalized', 'position', [0.72 0.965 0.01 0.03], 'callback', @show_core1, 'fontsize', size_font, 'value', 0);

%% draw second GUI

if ispc % windows switch
    fgui(2)                 = figure('toolbar', 'figure', 'name', 'FENCEGUI 2D', 'position', [1920 940 1 1], 'menubar', 'none', 'keypressfcn', @keypress2);
    ax(2)                   = subplot('position', [0.065 0.06 0.41 0.81]);
    ax(3)                   = subplot('position', [0.55 0.06 0.41 0.81]);
else
    fgui(2)                 = figure('toolbar', 'figure', 'name', 'FENCEGUI 2D', 'position', [1864 1100 1 1], 'menubar', 'none', 'keypressfcn', @keypress2);
    ax(2)                   = subplot('position', [0.065 0.06 0.41 0.81]);
    ax(3)                   = subplot('position', [0.55 0.06 0.41 0.81]);
end

set(ax(2:3), 'fontsize', size_font, 'layer', 'top')

axes(ax(2))
hold on
colormap(bone)
caxis([db_min(1) db_max(1)])
axis([dist_min_ref(1) dist_max_ref(1) elev_min_ref elev_max_ref])
box on
grid off
ylabel('(m)');
colorbar('fontsize', size_font)
% pan/zoom callbacks
h_pan(1)                    = pan;
set(h_pan(1), 'actionpostcallback', @panzoom1)
h_zoom(1)                   = zoom;
set(h_zoom(1), 'actionpostcallback', @panzoom1)

axes(ax(3))
hold on
colormap(bone)
caxis([db_min(2) db_max(2)])
axis([dist_min_ref(2) dist_max_ref(2) elev_min_ref elev_max_ref])
box on
grid off
colorbar('fontsize', size_font)
% pan/zoom callbacks
h_pan(2)                    = pan;
set(h_pan(2), 'actionpostcallback', @panzoom2)
h_zoom(2)                   = zoom;
set(h_zoom(2), 'actionpostcallback', @panzoom2)

% sliders
cb_min_slide(1)             = uicontrol(fgui(2), 'style', 'slider', 'units', 'normalized', 'position', [0.475 0.07 width_slide 0.32], 'callback', @slide_db_min1, 'min', db_min_ref(1), 'max', db_max_ref(1), ...
                                                 'value', db_min(1), 'sliderstep', [0.01 0.1]);
cb_max_slide(1)             = uicontrol(fgui(2), 'style', 'slider', 'units', 'normalized', 'position', [0.475 0.50 width_slide 0.32], 'callback', @slide_db_max1, 'min', db_min_ref(1), 'max', db_max_ref(1), ...
                                                 'value', db_max(1), 'sliderstep', [0.01 0.1]);
cb_min_slide(2)             = uicontrol(fgui(2), 'style', 'slider', 'units', 'normalized', 'position', [0.97 0.07 width_slide 0.32], 'callback', @slide_db_min2, 'min', db_min_ref(2), 'max', db_max_ref(2), ...
                                                 'value', db_min(2), 'sliderstep', [0.01 0.1]);
cb_max_slide(2)             = uicontrol(fgui(2), 'style', 'slider', 'units', 'normalized', 'position', [0.97 0.50 width_slide 0.32], 'callback', @slide_db_max2, 'min', db_min_ref(2), 'max', db_max_ref(2), ...
                                                 'value', db_max(2), 'sliderstep', [0.01 0.1]);
dist_min_slide(1)           = uicontrol(fgui(2), 'style', 'slider', 'units', 'normalized', 'position', [0.10 0.005 0.12 0.02], 'callback', @slide_dist_min1, 'min', dist_min_ref(1), 'max', dist_max_ref(1), ...
                                                 'value', dist_min_ref(1), 'sliderstep', [0.01 0.1]);
dist_max_slide(1)           = uicontrol(fgui(2), 'style', 'slider', 'units', 'normalized', 'position', [0.33 0.005 0.12 0.02], 'callback', @slide_dist_max1, 'min', dist_min_ref(1), 'max', dist_max_ref(1), ...
                                                 'value', dist_max_ref(1), 'sliderstep', [0.01 0.1]);
dist_min_slide(2)           = uicontrol(fgui(2), 'style', 'slider', 'units', 'normalized', 'position', [0.58 0.005 0.12 0.02], 'callback', @slide_dist_min2, 'min', dist_min_ref(2), 'max', dist_max_ref(2), ...
                                                 'value', dist_min_ref(2), 'sliderstep', [0.01 0.1]);
dist_max_slide(2)           = uicontrol(fgui(2), 'style', 'slider', 'units', 'normalized', 'position', [0.81 0.005 0.12 0.02], 'callback', @slide_dist_max2, 'min', dist_min_ref(2), 'max', dist_max_ref(2), ...
                                                 'value', dist_max_ref(2), 'sliderstep', [0.01 0.1]);
z_min_slide(2)              = uicontrol(fgui(2), 'style', 'slider', 'units', 'normalized', 'position', [0.005 0.07 width_slide 0.32], 'callback', @slide_z_min2, 'min', elev_min_ref, 'max', elev_max_ref, ...
                                                 'value', elev_min_ref, 'sliderstep', [0.01 0.1]);
z_max_slide(2)              = uicontrol(fgui(2), 'style', 'slider', 'units', 'normalized', 'position', [0.005 0.50 width_slide 0.32], 'callback', @slide_z_max2, 'min', elev_min_ref, 'max', elev_max_ref, ...
                                                 'value', elev_max_ref, 'sliderstep', [0.01 0.1]);
z_min_slide(3)              = uicontrol(fgui(2), 'style', 'slider', 'units', 'normalized', 'position', [0.505 0.07 width_slide 0.32], 'callback', @slide_z_min3, 'min', elev_min_ref, 'max', elev_max_ref, ...
                                                 'value', elev_min_ref, 'sliderstep', [0.01 0.1]);
z_max_slide(3)              = uicontrol(fgui(2), 'style', 'slider', 'units', 'normalized', 'position', [0.505 0.50 width_slide 0.32], 'callback', @slide_z_max3, 'min', elev_min_ref, 'max', elev_max_ref, ...
                                                 'value', elev_max_ref, 'sliderstep', [0.01 0.1]);

% slider values
cb_min_edit(1)              = annotation('textbox', [0.48 0.39 0.04 0.03], 'string', num2str(db_min(1)), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
cb_min_edit(2)              = annotation('textbox', [0.965 0.39 0.04 0.03], 'string', num2str(db_min(2)), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
cb_max_edit(1)              = annotation('textbox', [0.48 0.82 0.04 0.03], 'string', num2str(db_max(1)), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
cb_max_edit(2)              = annotation('textbox', [0.965 0.82 0.04 0.03], 'string', num2str(db_max(2)), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
dist_min_edit(1)            = annotation('textbox', [0.07 0.005 0.04 0.03], 'string', num2str(dist_min_ref(1)), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
dist_min_edit(2)            = annotation('textbox', [0.55 0.005 0.04 0.03], 'string', num2str(dist_min_ref(2)), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
dist_max_edit(1)            = annotation('textbox', [0.295 0.005 0.04 0.03], 'string', num2str(dist_max_ref(1)), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
dist_max_edit(2)            = annotation('textbox', [0.775 0.005 0.04 0.03], 'string', num2str(dist_max_ref(2)), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
z_min_edit(2)               = annotation('textbox', [0.005 0.39 0.04 0.03], 'string', num2str(elev_min_ref), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
z_min_edit(3)               = annotation('textbox', [0.51 0.39 0.04 0.03], 'string', num2str(elev_min_ref), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
z_max_edit(2)               = annotation('textbox', [0.005 0.82 0.04 0.03], 'string', num2str(elev_max_ref), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
z_max_edit(3)               = annotation('textbox', [0.51 0.82 0.04 0.03], 'string', num2str(elev_max_ref), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');

% push buttons
uicontrol(fgui(2), 'style', 'pushbutton', 'string', 'Load master data', 'units', 'normalized', 'position', [0.005 0.925 0.085 0.03], 'callback', @load_data1, 'fontsize', size_font, 'foregroundcolor', 'b')
uicontrol(fgui(2), 'style', 'pushbutton', 'string', 'Load intersecting data', 'units', 'normalized', 'position', [0.345 0.925 0.11 0.03], 'callback', @load_data2, 'fontsize', size_font, 'foregroundcolor', 'b')
uicontrol(fgui(2), 'style', 'pushbutton', 'string', 'Match', 'units', 'normalized', 'position', [0.69 0.925 0.04 0.03], 'callback', @pk_match, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(fgui(2), 'style', 'pushbutton', 'string', 'Unmatch', 'units', 'normalized', 'position', [0.73 0.925 0.05 0.03], 'callback', @pk_unmatch, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(fgui(2), 'style', 'pushbutton', 'string', 'Choose', 'units', 'normalized', 'position', [0.17 0.925 0.045 0.03], 'callback', @choose_pk1, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(fgui(2), 'style', 'pushbutton', 'string', 'Choose', 'units', 'normalized', 'position', [0.535 0.925 0.045 0.03], 'callback', @choose_pk2, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(fgui(2), 'style', 'pushbutton', 'string', 'Next', 'units', 'normalized', 'position', [0.245 0.925 0.03 0.03], 'callback', @pk_next3, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(fgui(2), 'style', 'pushbutton', 'string', 'Next', 'units', 'normalized', 'position', [0.61 0.925 0.03 0.03], 'callback', @pk_next4, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(fgui(2), 'style', 'pushbutton', 'string', 'Last', 'units', 'normalized', 'position', [0.215 0.925 0.03 0.03], 'callback', @pk_last3, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(fgui(2), 'style', 'pushbutton', 'string', 'Last', 'units', 'normalized', 'position', [0.58 0.925 0.03 0.03], 'callback', @pk_last4, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(fgui(2), 'style', 'pushbutton', 'string', 'test', 'units', 'normalized', 'position', [0.91 0.925 0.03 0.03], 'callback', @misctest, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(fgui(2), 'style', 'pushbutton', 'string', 'reset x/z', 'units', 'normalized', 'position', [0.41 0.885 0.05 0.03], 'callback', @reset_xz1, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(fgui(2), 'style', 'pushbutton', 'string', 'reset x/z', 'units', 'normalized', 'position', [0.895 0.885 0.05 0.03], 'callback', @reset_xz2, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(fgui(2), 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.005 0.03 0.03 0.03], 'callback', @reset_z_min2, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(fgui(2), 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.005 0.46 0.03 0.03], 'callback', @reset_z_max2, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(fgui(2), 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.225 0.005 0.03 0.03], 'callback', @reset_dist_min1, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(fgui(2), 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.455 0.005 0.03 0.03], 'callback', @reset_dist_max1, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(fgui(2), 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.48 0.03 0.03 0.03], 'callback', @reset_db_min1, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(fgui(2), 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.48 0.46 0.03 0.03], 'callback', @reset_db_max1, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(fgui(2), 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.51 0.03 0.03 0.03], 'callback', @reset_z_min3, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(fgui(2), 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.51 0.46 0.03 0.03], 'callback', @reset_z_max3, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(fgui(2), 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.705 0.005 0.03 0.03], 'callback', @reset_dist_min2, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(fgui(2), 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.935 0.005 0.03 0.03], 'callback', @reset_dist_max2, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(fgui(2), 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.965 0.03 0.03 0.03], 'callback', @reset_db_min2, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(fgui(2), 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.965 0.46 0.03 0.03], 'callback', @reset_db_max2, 'fontsize', size_font, 'foregroundcolor', 'r')

% fixed text annotations
b(1)                        = annotation('textbox', [0.10 0.925 0.03 0.03], 'string', 'N_{decim}', 'fontsize', size_font, 'color', 'b', 'edgecolor', 'none');
b(2)                        = annotation('textbox', [0.465 0.925 0.03 0.03], 'string', 'N_{decim}', 'fontsize', size_font, 'color', 'b', 'edgecolor', 'none');
b(3)                        = annotation('textbox', [0.49 0.42 0.03 0.03], 'string', 'dB/z_{min}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
b(4)                        = annotation('textbox', [0.49 0.85 0.03 0.03], 'string', 'dB/z_{max}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
b(5)                        = annotation('textbox', [0.965 0.42 0.03 0.03], 'string', 'dB_{min}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
b(6)                        = annotation('textbox', [0.965 0.85 0.03 0.03], 'string', 'dB_{max}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
b(7)                        = annotation('textbox', [0.18 0.965 0.03 0.03], 'string', 'Layer', 'fontsize', size_font, 'color', 'm', 'edgecolor', 'none');
b(8)                        = annotation('textbox', [0.515 0.965 0.03 0.03], 'string', 'Layer', 'fontsize', size_font, 'color', 'm', 'edgecolor', 'none');
b(9)                        = annotation('textbox', [0.37 0.88 0.03 0.03], 'string', 'Grid', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
b(10)                       = annotation('textbox', [0.855 0.88 0.03 0.03], 'string', 'Grid', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
b(11)                       = annotation('textbox', [0.035 0.005 0.03 0.03], 'string', 'dist_{min}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
b(12)                       = annotation('textbox', [0.515 0.005 0.03 0.03], 'string', 'dist_{min}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
b(13)                       = annotation('textbox', [0.26 0.005 0.03 0.03], 'string', 'dist_{max}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
b(14)                       = annotation('textbox', [0.74 0.005 0.03 0.03], 'string', 'dist_{max}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
b(15)                       = annotation('textbox', [0.005 0.42 0.03 0.03], 'string', 'z_{min}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
b(16)                       = annotation('textbox', [0.005 0.85 0.03 0.03], 'string', 'z_{max}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
b(17)                       = annotation('textbox', [0.005 0.88 0.03 0.03], 'string', 'fix', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
b(18)                       = annotation('textbox', [0.525 0.88 0.03 0.03], 'string', 'fix', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
b(19)                       = annotation('textbox', [0.485 0.005 0.03 0.03], 'string', 'fix', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
b(20)                       = annotation('textbox', [0.965 0.005 0.03 0.03], 'string', 'fix', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
b(21)                       = annotation('textbox', [0.95 0.88 0.03 0.03], 'string', 'fix 1', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
b(22)                       = annotation('textbox', [0.98 0.88 0.03 0.03], 'string', '2', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
b(23)                       = annotation('textbox', [0.465 0.88 0.03 0.03], 'string', 'fix 1', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
b(24)                       = annotation('textbox', [0.4925 0.88 0.03 0.03], 'string', '2', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
b(25)                       = annotation('textbox', [0.04 0.88 0.08 0.03], 'string', 'Intersection #', 'fontsize', size_font, 'color', 'b', 'edgecolor', 'none');
b(26)                       = annotation('textbox', [0.15 0.88 0.08 0.03], 'string', 'Core', 'fontsize', size_font, 'color', 'b', 'edgecolor', 'none');
b(27)                       = annotation('textbox', [0.56 0.88 0.08 0.03], 'string', 'Intersections', 'fontsize', size_font, 'color', 'b', 'edgecolor', 'none');
b(28)                       = annotation('textbox', [0.635 0.88 0.08 0.03], 'string', 'Core', 'fontsize', size_font, 'color', 'b', 'edgecolor', 'none');
b(29)                       = annotation('textbox', [0.785 0.925 0.08 0.03], 'string', 'Nearest', 'fontsize', size_font, 'color', 'm', 'edgecolor', 'none');
if ~ispc
    set(b, 'fontweight', 'bold')
end

% variable text annotations
file_box(2)                 = annotation('textbox', [0.005 0.965 0.16 0.03], 'string', '', 'color', 'k', 'fontsize', size_font, 'backgroundcolor', 'w', 'edgecolor', 'k', 'interpreter', 'none');
file_box(3)                 = annotation('textbox', [0.345 0.965 0.16 0.03], 'string', '', 'color', 'k', 'fontsize', size_font, 'backgroundcolor', 'w', 'edgecolor', 'k', 'interpreter', 'none');
status_box(2)               = annotation('textbox', [0.69 0.965 0.30 0.03], 'string', '', 'color', 'k', 'fontsize', size_font, 'backgroundcolor', 'w', 'edgecolor', 'k', 'interpreter', 'none');

disp_group(1)               = uibuttongroup('position', [0.26 0.965 0.08 0.03], 'selectionchangefcn', @disp_radio1);
uicontrol(fgui(2), 'style', 'text', 'parent', disp_group(1), 'units', 'normalized', 'position', [0 0.6 0.9 0.3], 'fontsize', size_font)
disp_check(1, 1)            = uicontrol(fgui(2), 'style', 'radio', 'string', 'amp.', 'units', 'normalized', 'position', [0.01 0.1 0.45 0.8], 'parent', disp_group(1), 'fontsize', size_font, 'handlevisibility', 'off');
disp_check(1, 2)            = uicontrol(fgui(2), 'style', 'radio', 'string', 'flat', 'units', 'normalized', 'position', [0.51 0.1 0.45 0.8], 'parent', disp_group(1), 'fontsize', size_font, ...
                                                 'handlevisibility', 'off', 'visible', 'off');
set(disp_group(1), 'selectedobject', disp_check(1, 1))

disp_group(2)               = uibuttongroup('position', [0.60 0.965 0.08 0.03], 'selectionchangefcn', @disp_radio2);
uicontrol(fgui(2), 'style', 'text', 'parent', disp_group(2), 'units', 'normalized', 'position', [0 0.6 0.9 0.3], 'fontsize', size_font)
disp_check(2, 1)            = uicontrol(fgui(2), 'style', 'radio', 'string', 'amp.', 'units', 'normalized', 'position', [0.01 0.1 0.45 0.8], 'parent', disp_group(2), 'fontsize', size_font, 'handlevisibility', 'off');
disp_check(2, 2)            = uicontrol(fgui(2), 'style', 'radio', 'string', 'flat', 'units', 'normalized', 'position', [0.51 0.1 0.45 0.8], 'parent', disp_group(2), 'fontsize', size_font, ...
                                                 'handlevisibility', 'off', 'visible', 'off');
set(disp_group(2), 'selectedobject', disp_check(2, 1))

% value boxes
decim_edit(2, 1)            = uicontrol(fgui(2), 'style', 'edit', 'string', num2str(decim(2, 1)), 'units', 'normalized', 'position', [0.135 0.925 0.03 0.03], 'fontsize', size_font, 'foregroundcolor', 'k', ...
                                                 'backgroundcolor', 'w', 'callback', @adj_decim3);
decim_edit(2, 2)            = uicontrol(fgui(2), 'style', 'edit', 'string', num2str(decim(2, 2)), 'units', 'normalized', 'position', [0.50 0.925 0.03 0.03], 'fontsize', size_font, 'foregroundcolor', 'k', ...
                                                 'backgroundcolor', 'w', 'callback', @adj_decim4);
% menus
cmap_list                   = uicontrol(fgui(2), 'style', 'popupmenu', 'string', cmaps, 'value', 1, 'units', 'normalized', 'position', [0.945 0.925 0.05 0.03], 'callback', @change_cmap, 'fontsize', size_font);
layer_list(2, 1)            = uicontrol(fgui(2), 'style', 'popupmenu', 'string', 'N/A', 'value', 1, 'units', 'normalized', 'position', [0.21 0.955 0.05 0.04], 'fontsize', size_font, 'foregroundcolor', 'k', ...
                                                 'callback', @choose_layer3);
layer_list(2, 2)            = uicontrol(fgui(2), 'style', 'popupmenu', 'string', 'N/A', 'value', 1, 'units', 'normalized', 'position', [0.545 0.955 0.05 0.04], 'fontsize', size_font, 'foregroundcolor', 'k', ...
                                                 'callback', @choose_layer4);
intnum_list                 = uicontrol(fgui(2), 'style', 'popupmenu', 'string', 'N/A', 'value', 1, 'units', 'normalized', 'position', [0.10 0.86 0.04 0.05], 'fontsize', size_font, 'foregroundcolor', 'k', ...
                                                 'callback', @change_int);
data_list(1)                = uicontrol(fgui(2), 'style', 'popupmenu', 'string', 'N/A', 'value', 1, 'units', 'normalized', 'position', [0.19 0.86 0.175 0.05], 'fontsize', size_font, 'foregroundcolor', 'k');
data_list(2)                = uicontrol(fgui(2), 'style', 'popupmenu', 'string', 'N/A', 'value', 1, 'units', 'normalized', 'position', [0.68 0.86 0.175 0.05], 'fontsize', size_font, 'foregroundcolor', 'k');

% check boxes
distfix_check(1)            = uicontrol(fgui(2), 'style', 'checkbox', 'units', 'normalized', 'position', [0.50 0.005 0.01 0.03], 'fontsize', size_font, 'value', 0);
distfix_check(2)            = uicontrol(fgui(2), 'style', 'checkbox', 'units', 'normalized', 'position', [0.98 0.005 0.01 0.03], 'fontsize', size_font, 'value', 0);
zfix_check(2)               = uicontrol(fgui(2), 'style', 'checkbox', 'units', 'normalized', 'position', [0.02 0.88 0.01 0.03], 'fontsize', size_font, 'value', 0);
zfix_check(3)               = uicontrol(fgui(2), 'style', 'checkbox', 'units', 'normalized', 'position', [0.54 0.88 0.01 0.03], 'fontsize', size_font, 'value', 0);
cbfix_check1(1)             = uicontrol(fgui(2), 'style', 'checkbox', 'units', 'normalized', 'position', [0.485 0.88 0.01 0.03], 'fontsize', size_font, 'value', 0);
cbfix_check1(2)             = uicontrol(fgui(2), 'style', 'checkbox', 'units', 'normalized', 'position', [0.97 0.88 0.01 0.03], 'fontsize', size_font, 'value', 0);
cbfix_check2(1)             = uicontrol(fgui(2), 'style', 'checkbox', 'units', 'normalized', 'position', [0.50 0.88 0.01 0.03], 'callback', @narrow_cb1, 'fontsize', size_font, 'value', 0);
cbfix_check2(2)             = uicontrol(fgui(2), 'style', 'checkbox', 'units', 'normalized', 'position', [0.985 0.88 0.01 0.03], 'callback', @narrow_cb2, 'fontsize', size_font, 'value', 0);
grid_check(2)               = uicontrol(fgui(2), 'style', 'checkbox', 'units', 'normalized', 'position', [0.39 0.88 0.01 0.03], 'callback', @toggle_grid2, 'fontsize', size_font, 'value', 0);
grid_check(3)               = uicontrol(fgui(2), 'style', 'checkbox', 'units', 'normalized', 'position', [0.88 0.88 0.01 0.03], 'callback', @toggle_grid3, 'fontsize', size_font, 'value', 0);
pk_check(2, 1)              = uicontrol(fgui(2), 'style', 'checkbox', 'units', 'normalized', 'position', [0.17 0.965 0.01 0.03], 'callback', @show_pk3, 'fontsize', size_font, 'value', 0);
pk_check(2, 2)              = uicontrol(fgui(2), 'style', 'checkbox', 'units', 'normalized', 'position', [0.505 0.965 0.01 0.03], 'callback', @show_pk4, 'fontsize', size_font, 'value', 0);
data_check(1)               = uicontrol(fgui(2), 'style', 'checkbox', 'units', 'normalized', 'position', [0.09 0.925 0.01 0.03], 'callback', @show_data1, 'fontsize', size_font, 'value', 0);
data_check(2)               = uicontrol(fgui(2), 'style', 'checkbox', 'units', 'normalized', 'position', [0.455 0.925 0.01 0.03], 'callback', @show_data2, 'fontsize', size_font, 'value', 0);
int_check(2)                = uicontrol(fgui(2), 'style', 'checkbox', 'units', 'normalized', 'position', [0.14 0.88 0.01 0.03], 'callback', @show_int2, 'fontsize', size_font, 'value', 0);
int_check(3)                = uicontrol(fgui(2), 'style', 'checkbox', 'units', 'normalized', 'position', [0.625 0.88 0.01 0.03], 'callback', @show_int3, 'fontsize', size_font, 'value', 0);
core_check(2)               = uicontrol(fgui(2), 'style', 'checkbox', 'units', 'normalized', 'position', [0.175 0.88 0.01 0.03], 'callback', @show_core2, 'fontsize', size_font, 'value', 0);
core_check(3)               = uicontrol(fgui(2), 'style', 'checkbox', 'units', 'normalized', 'position', [0.66 0.88 0.01 0.03], 'callback', @show_core3, 'fontsize', size_font, 'value', 0);
nearest_check               = uicontrol(fgui(2), 'style', 'checkbox', 'units', 'normalized', 'position', [0.825 0.925 0.01 0.03], 'fontsize', size_font, 'value', 1);

figure(fgui(1))

linkprop(layer_list(:, 1), {'value' 'string'});
linkprop(layer_list(:, 2), {'value' 'string'});

%% Clear plots

    function clear_plots(source, eventdata)
        if (any(p_bed(:, curr_rad)) && any(ishandle(p_bed(:, curr_rad))))
            delete(p_bed((logical(p_bed(:, curr_rad)) & ishandle(p_bed(:, curr_rad))), curr_rad))
        end
        if (logical(p_bedflat(curr_rad)) && ishandle(p_bedflat(curr_rad)))
            delete(p_bedflat(curr_rad))
        end
        if (any(p_coreflat{curr_rad}) && any(ishandle(p_coreflat{curr_rad})))
            delete(p_coreflat{curr_rad}(logical(p_coreflat{curr_rad}) & ishandle(p_coreflat{curr_rad})))
        end
        if (any(p_corenameflat{curr_rad}) && any(ishandle(p_corenameflat{curr_rad})))
            delete(p_corenameflat{curr_rad}(logical(p_corenameflat{curr_rad}) & ishandle(p_corenameflat{curr_rad})))
        end
        if (logical(p_data(curr_rad)) && ishandle(p_data(curr_rad)))
            delete(p_data(curr_rad))
        end
        for ii = 1:2
            for jj = 1:3
                if (any(p_int1{ii, jj}) && any(ishandle(p_int1{ii, jj})))
                    delete(p_int1{ii, jj}(logical(p_int1{ii, jj}) & ishandle(p_int1{ii, jj})))
                end
            end
            if (any(p_core{ii, curr_rad}) && any(ishandle(p_core{ii, curr_rad})))
                delete(p_core{ii, curr_rad}(logical(p_core{ii, curr_rad}) & ishandle(p_core{ii, curr_rad})))
            end
            if (any(p_corename{ii, curr_rad}) && any(ishandle(p_corename{ii, curr_rad})))
                delete(p_corename{ii, curr_rad}(logical(p_corename{ii, curr_rad}) & ishandle(p_corename{ii, curr_rad})))
            end
            if (any(p_int2{ii, curr_rad}) && any(ishandle(p_int2{ii, curr_rad})))
                delete(p_int2{ii, curr_rad}(logical(p_int2{ii, curr_rad}) & ishandle(p_int2{ii, curr_rad})))
            end
            if (any(p_pk{ii, curr_rad}) && any(ishandle(p_pk{ii, curr_rad})))
                delete(p_pk{ii, curr_rad}(logical(p_pk{ii, curr_rad}) & ishandle(p_pk{ii, curr_rad})))
            end
        end
        if (any(p_pkflat{curr_rad}) && any(ishandle(p_pkflat{curr_rad})))
            delete(p_pkflat{curr_rad}(logical(p_pkflat{curr_rad}) & ishandle(p_pkflat{curr_rad})))
        end
        if (any(p_surf(:, curr_rad)) && any(ishandle(p_surf(:, curr_rad))))
            delete(p_surf((logical(p_surf(:, curr_rad)) & ishandle(p_surf(:, curr_rad))), curr_rad))
        end
        set(file_box(1 + curr_rad), 'string', '')
        if (curr_rad == 1)
            set(file_box(1), 'string', '')
        end
        set([pk_check(:, curr_rad); data_check(curr_rad); int_check'], 'value', 0)
        set(disp_group(curr_rad), 'selectedobject', disp_check(curr_rad, 1))
        set([layer_list(:, curr_rad); intnum_list; data_list(curr_rad)], 'string', 'N/A', 'value', 1)
        axes(ax(curr_ax))
    end

%% Clear data and picks

    function clear_data(source, eventdata)
        [bed_avail(curr_rad), data_done(curr_rad), flat_done(curr_rad), pk_done(curr_rad), surf_avail(curr_rad)] ...
                            = deal(false);
        [amp_mean{curr_rad}, colors{curr_rad}, depth{curr_rad}, depth_bed{curr_rad}, depth_bed_flat{curr_rad}, depth_flat{curr_rad}, depth_layer_flat{curr_rad}, depth_mat{curr_rad}, elev{curr_rad}, ...
            file_pk_short{curr_rad}, ind_decim{1, curr_rad}, ind_decim{2, curr_rad}, ind_corr{curr_rad}, ind_int_core{curr_rad}, layer_str{curr_rad}, pk{curr_rad}, p_core{1, curr_rad}, p_core{2, curr_rad}, ...
            p_corename{1, curr_rad}, p_corename{2, curr_rad}, p_coreflat{curr_rad}, p_corenameflat{curr_rad}, p_int1{1, 1}, p_int1{1, 2}, p_int1{1, 3}, p_int1{2, 1}, p_int1{2, 2}, p_int1{2, 3}, p_int2{1, curr_rad}, ...
            p_int2{2, curr_rad}, p_pk{1, curr_rad}, p_pk{2, curr_rad}, p_pkflat{curr_rad}, twtt{curr_rad}] ...
                            = deal([]);
        [curr_layer(curr_rad), curr_trans(curr_rad), curr_year(curr_rad), dt(curr_rad), num_data(curr_rad), num_decim(:, curr_rad), num_sample(curr_rad), p_bed(curr_rad), p_bedflat(curr_rad), p_data(curr_rad), ...
            p_surf(curr_rad)] ...
                            = deal(0);
        [curr_ind_int, ii, ind_x_pk, ind_y_pk, jj, num_int, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8] ...
                            = deal(0);
    end

%% Load intersection data

    function load_int(source, eventdata)
        [curr_gui, curr_ax] = deal(1);
        if int_done
            set(status_box(1), 'string', 'Transect intersections already loaded.')
            return
        end
        if ispc
            if exist('\\melt\icebridge\data\mat\', 'dir')
                path_ref    = '\\melt\icebridge\data\mat\';
            end
        else
            if exist('/Volumes/melt/icebridge/data/mat/', 'dir')
                path_ref    = '/Volumes/melt/icebridge/data/mat/';
            end
        end
        if (~isempty(path_ref) && exist([path_ref 'xy_int.mat'], 'file'))
            file_ref       = 'xy_int.mat';
        else
            % Dialog box to choose picks file to load
            if ~isempty(path_ref)
                [file_ref, path_ref] = uigetfile('*.mat', 'Load intersection data (xy_int.mat):', path_ref);
            elseif ~isempty(path_pk{curr_rad})
                [file_ref, path_ref] = uigetfile('*.mat', 'Load intersection data (xy_int.mat):', path_pk{curr_rad});
            elseif ~isempty(path_data{curr_rad})
                [file_ref, path_ref] = uigetfile('*.mat', 'Load intersection data (xy_int.mat):', path_data{curr_rad});
            else
                [file_ref, path_ref] = uigetfile('*.mat', 'Load intersection data (xy_int.mat):');
            end
            if isnumeric(file_ref)
                [file_ref, path_ref] = deal('');
            end
        end
        if ~isempty(file_ref)
            tmp1            = load([path_ref file_ref]);
            try
                [name_trans, int_all, num_year, num_trans] ...
                            = deal(tmp1.name_trans, tmp1.int_all, tmp1.num_year, tmp1.num_trans);
            catch
               set(status_box(1), 'string', 'Chosen file does not contain expected intersection data.')
               return
            end
            int_done        = true;
            set(status_box(1), 'string', 'Intersection data loaded.')
        else
            set(status_box(1), 'string', 'No intersection data loaded.')
        end
    end

%% Load core intersection data

    function load_core(source, eventdata)
        
        [curr_gui, curr_ax] = deal(1);
        if core_done
            set(status_box(1), 'string', 'Core intersections already loaded.')
            return
        end
        if ~int_done
            set(status_box(1), 'string', 'Load transect intersections first.')
            return
        end
        
        if (~isempty(path_ref) && exist([path_ref 'core_int.mat'], 'file'))
            [file_core, path_core] ...
                            = deal('core_int.mat', path_ref);
        elseif ~isempty(path_core)
            [file_core, path_core] = uigetfile('*.mat', 'Load core intersections (core_int.mat):', path_core);
        elseif ~isempty(path_ref)
            [file_core, path_core] = uigetfile('*.mat', 'Load core intersections (core_int.mat):', path_ref);
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
        
        if ~isempty(file_core)
            
            set(status_box(1), 'string', 'Loading core intersections...')
            pause(0.1)
            
            % load core intersection file
            tmp1        = load([path_core file_core]);
            try
                [int_core, name_core, rad_threshold, x_core, y_core] ...
                        = deal(tmp1.int_core, tmp1.name_core, tmp1.rad_threshold, tmp1.x_core, tmp1.y_core);
            catch % give up, force restart
                set(status_box(1), 'string', [file_core ' does not contain the expected variables. Try again.'])
                return
            end
            
            core_done   = true;
            set(status_box(1), 'string', 'Core intersection data loaded.')
            return
            
        else
            set(status_box(1), 'string', 'No core intersections loaded.')
        end
    end

%% Load core breakout

    function load_core_breakout(source, eventdata)
        
        for ii = 1:2
            for jj = 1:2
                if (any(p_core{jj, ii}) && any(ishandle(p_core{jj, ii})))
                    delete(p_core{jj, ii}(logical(p_core{jj, ii}) & ishandle(p_core{jj, ii})))
                end
                if (any(p_corename{jj, ii}) && any(ishandle(p_corename{jj, ii})))
                    delete(p_corename{jj, ii}(logical(p_corename{jj, ii}) & ishandle(p_corename{jj, ii})))
                end
            end
            if (any(p_coreflat{ii}) && any(ishandle(p_coreflat{ii})))
                delete(p_coreflat{ii}(logical(p_coreflat{ii}) & ishandle(p_coreflat{ii})))
            end
            if (any(p_corenameflat{ii}) && any(ishandle(p_corenameflat{ii})))
                delete(p_corenameflat{ii}(logical(p_corenameflat{ii}) & ishandle(p_corenameflat{ii})))
            end
        end
        
        for ii = 1:2
            if ~isempty(int_core{curr_year(ii)}{curr_trans(ii)})
                
                ind_int_core{ii}    = [];
                for jj = 1:size(int_core{curr_year(ii)}{curr_trans(ii)}, 1)
                    try %#ok<TRYNC>
                        [tmp1, tmp2]= min(sqrt(((pk{ii}.x - int_core{curr_year(ii)}{curr_trans(ii)}(jj, 4)) .^ 2) + ((pk{ii}.y - int_core{curr_year(ii)}{curr_trans(ii)}(jj, 5)) .^ 2)));
                        if (tmp1 < rad_threshold)
                            ind_int_core{ii} ...
                                    = [ind_int_core{ii} tmp2];
                        end
                    end
                end
                
                if isempty(ind_int_core{ii})
                    continue
                else
                    num_int_core(ii)= length(ind_int_core{ii});
                end
                
                for jj = 1:2
                    [p_core{jj, ii}, p_corename{jj, ii}] = deal(zeros(1, num_int_core(ii)));
                    for kk = 1:num_int_core(ii)
                        if (jj == 1)
                            axes(ax(1)) %#ok<LAXES>
                            p_core{jj, ii}(kk) = plot3(repmat(x_core(int_core{curr_year(ii)}{curr_trans(ii)}(kk, 3)), 1, 2), repmat(y_core(int_core{curr_year(ii)}{curr_trans(ii)}(kk, 3)), 1, 2), ...
                                                       [elev_min_ref elev_max_ref], 'color', 'k', 'linewidth', 2, 'visible', 'off');
                            p_corename{jj, ii}(kk) = text((x_core(int_core{curr_year(ii)}{curr_trans(ii)}(kk, 3)) + 1), (y_core(int_core{curr_year(ii)}{curr_trans(ii)}(kk, 3)) + 1), (elev_max_ref - 50), ...
                                                          name_core{int_core{curr_year(ii)}{curr_trans(ii)}(kk, 3)}, 'color', 'k', 'fontsize', size_font, 'visible', 'off');
                        else
                            axes(ax(ii + 1)) %#ok<LAXES>
                            p_core{jj, ii}(kk) = plot(repmat(pk{ii}.dist_lin(ind_int_core{ii}(kk)), 1, 2), [elev_min_ref elev_max_ref], 'color', [0.5 0.5 0.5], 'linewidth', 2, 'visible', 'off');
                            p_corename{jj, ii}(kk) = text((pk{ii}.dist_lin(ind_int_core{ii}(kk)) + 1), (elev_max_ref - 50), name_core{int_core{curr_year(ii)}{curr_trans(ii)}(kk, 3)}, ...
                                                          'color', [0.5 0.5 0.5], 'fontsize', size_font, 'visible', 'off');
                        end
                    end
                end
                
                if flat_done(ii)
                    [p_coreflat{ii}, p_corenameflat{ii}] = deal(zeros(1, num_int_core(ii)));
                    axes(ax(ii + 1)) %#ok<LAXES>
                    for jj = 1:num_int_core(ii)
                        p_coreflat{ii}(jj) = plot(repmat(pk{ii}.dist_lin(ind_int_core{ii}(jj)), 1, 2), [depth_min_ref depth_max_ref], 'color', [0.5 0.5 0.5], 'linewidth', 2, 'visible', 'off');
                        p_corenameflat{ii}(jj) = text((pk{ii}.dist_lin(ind_int_core{ii}(jj)) + 1), (depth_min_ref + 50), name_core{int_core{curr_year(ii)}{curr_trans(ii)}(jj, 3)}, ...
                                                       'color', [0.5 0.5 0.5], 'fontsize', size_font, 'visible', 'off');
                    end
                end
            end
        end
        
        set(status_box, 'string', ['Core intersections loaded. ' num2str(num_int_core(1)) '/' num2str(num_int_core(2)) ' for these transects within ' num2str(rad_threshold) ' km.'])        
        core_done       = true;
        set(core_check, 'value', 1)
        show_core3
        show_core2
        show_core1
    end

%% Load picks for this transect

    function load_pk1(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                        = deal(1, 1, 1, 2);
        load_pk
    end

    function load_pk2(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                        = deal(1, 1, 2, 1);
        load_pk
    end

    function load_pk(source, eventdata)
        
        if ~int_done
            set(status_box(1), 'string', 'Load transect intersection data before loading master picks.')
            return
        end
        if ~core_done
            set(status_box(1), 'string', 'Load core intersection data before loading picks.')
            return
        end
        if ((curr_rad == 2) && ~pk_done(1))
            set(status_box(1), 'string', 'Load master transect before intersecting transect.')
            return
        end
        
        % Dialog box to choose picks file to load
        if ~isempty(path_pk{curr_rad})
            [file_pk{curr_rad}, path_pk{curr_rad}] ...
                            = uigetfile('*.mat', 'Load merged picks:', path_pk{curr_rad});
        elseif ~isempty(path_data{curr_rad})
            [file_pk{curr_rad}, path_pk{curr_rad}] ...
                            = uigetfile('*.mat', 'Load merged picks:', path_data{curr_rad});
        elseif ~isempty(path_pk{curr_rad_alt})
            [file_pk{curr_rad}, path_pk{curr_rad}] ...
                            = uigetfile('*.mat', 'Load merged picks:', path_pk{curr_rad_alt});
        elseif ~isempty(path_data{curr_rad_alt})
            [file_pk{curr_rad}, path_pk{curr_rad}] ...
                            = uigetfile('*.mat', 'Load merged picks:', path_data{curr_rad_alt});
        elseif ~isempty(path_ref)
            [file_pk{curr_rad}, path_pk{curr_rad}] ...
                            = uigetfile('*.mat', 'Load merged picks:', path_ref);
        else
            [file_pk{curr_rad}, path_pk{curr_rad}] ...
                            = uigetfile('*.mat', 'Load merged picks:');
        end
        if isnumeric(file_pk{curr_rad})
            [file_pk{curr_rad}, path_pk{curr_rad}] ...
                            = deal('');
        end
        
        if ~isempty(file_pk{curr_rad})
            
            % get rid of radargrams if present
            clear_plots
            clear_data
            if ((curr_rad == 1) && pk_done(2))
                curr_rad    = 2;
                clear_plots
                clear_data
                set(int_list, 'string', 'N/A', 'value', 1)
                curr_rad    = 1;
            end
            
            axes(ax(curr_ax))
            
            pause(0.1)
            
            if (curr_rad == 2)
                tmp1        = get(int_list, 'string');
                if ~strcmp(tmp1{get(int_list, 'value')}(1:11), file_pk{curr_rad}(1:11))
                    set(status_box(1), 'string', 'Chosen picks filename does not match selected transect. Continue loading? Y: yes; otherwise: no...')
                    waitforbuttonpress
                    if ~strcmpi(get(fgui(1), 'currentcharacter'), 'Y')
                        set(status_box(1), 'string', 'Loading of picks file cancelled by user.')
                        return
                    end
                end
            end
            
            % load picks files
            set(status_box(1), 'string', ['Loading ' file_pk{curr_rad}(1:(end - 4)) '...'])
            pause(0.1)
            tmp1            = load([path_pk{curr_rad} file_pk{curr_rad}]);
            try
                pk{curr_rad}= tmp1.pk;
                tmp1        = 0;
                if ~isfield(pk{curr_rad}, 'merge_flag')
                    set(status_box(1), 'string', 'Load merged picks files only.')
                    return
                end
                if ((curr_rad == 1) && ~isfield(pk{curr_rad}, 'ind_layer'))
                    set(status_box(1), 'string', 'Primary picks must already be tied to a core-intersecting transect.')
                    return
                end
                % check to see if surface and bed picks are available
                if isfield(pk{curr_rad}, 'elev_surf')
                    surf_avail(curr_rad) ...
                            = true;
                end
                if isfield(pk{curr_rad}, 'elev_bed')
                    bed_avail(curr_rad) ...
                            = true;
                end
                if isfield(pk{curr_rad}, 'poly_flat_merge')
                    flat_done(curr_rad) ...
                            = true;
                    set(disp_check(curr_rad, 2), 'visible', 'on')
                else
                    set(disp_check(curr_rad, 2), 'visible', 'off')
                end
            catch % give up, force restart
                set(status_box(1), 'string', [file_pk{curr_rad} ' does not contain a pk structure. Try again.'])
                return
            end
            
            if ~any(surf_avail(curr_rad))
                set(status_box(1), 'string', 'Surface pick not included in pick file, which will create problems.')
            end
            
            % add master picks matrix if not already present
            if ((curr_rad == 2) && ~isfield(pk{curr_rad}, 'ind_layer'))
                pk{curr_rad}.ind_layer ...
                            = [];
            end
            
            % extract date and best name from pk files
            [tmp1, tmp2]    = strtok(file_pk{curr_rad}, '_');
            file_pk_short{curr_rad} ...
                            = [tmp1 tmp2(1:3)];
            if ~strcmp(tmp2(4), '_')
                file_pk_short{curr_rad} ...
                            = [file_pk_short{curr_rad} tmp2(4)];
            end
            if (curr_rad == 1)
                set(file_box([1 2]), 'string', file_pk_short{curr_rad})
            else
                set(file_box(3), 'string', file_pk_short{curr_rad})
            end
            
            % decimated vectors for display
            if (decim(1, curr_rad) > 1)
                ind_decim{1, curr_rad} ...
                            = (1 + ceil(decim(1, curr_rad) / 2)):decim(1, curr_rad):(pk{curr_rad}.num_trace_tot - ceil(decim(1, curr_rad) / 2));
            else
                ind_decim{1, curr_rad} ...
                            = 1:pk{curr_rad}.num_trace_tot;
            end
            num_decim(1, curr_rad) ...
                            = length(ind_decim{1, curr_rad});
            if (decim(2, curr_rad) > 1)
                ind_decim{2, curr_rad} ...
                            = (1 + ceil(decim(2, curr_rad) / 2)):decim(2, curr_rad):(pk{curr_rad}.num_trace_tot - ceil(decim(2, curr_rad) / 2));
            else
                ind_decim{2, curr_rad} ...
                            = 1:pk{curr_rad}.num_trace_tot;
            end
            num_decim(2, curr_rad) ...
                            = length(ind_decim{2, curr_rad});
            
            % determine current year/transect
            tmp1            = file_pk_short{curr_rad};
            if isnan(str2double(tmp1(end)))
                tmp1        = tmp1(1:(end - 1));
            end
            for ii = 1:num_year
                for jj = 1:num_trans(ii)
                    if strcmp(tmp1, name_trans{ii}{jj})
                        break
                    end
                end
                if strcmp(tmp1, name_trans{ii}{jj})
                    break
                end
            end
            [curr_year(curr_rad), curr_trans(curr_rad)] ...
                            = deal(ii, jj);
            
            % find intersections for primary transect
            if (curr_rad == 1)
                tmp1        = find((int_all(:, 1) == curr_year(curr_rad)) & (int_all(:, 2) == curr_trans(curr_rad)));
                tmp2        = find((int_all(:, 5) == curr_year(curr_rad)) & (int_all(:, 6) == curr_trans(curr_rad)));
                tmp3        = [];
                for ii = 1:length(tmp1)
                    tmp3    = [tmp3; name_trans{int_all(tmp1(ii), 5)}(int_all(tmp1(ii), 6))]; %#ok<AGROW>
                end
                for ii = 1:length(tmp2)
                    tmp3    = [tmp3; name_trans{int_all(tmp2(ii), 1)}(int_all(tmp2(ii), 2))]; %#ok<AGROW>
                end
                tmp4        = unique(tmp3);
                for ii = 1:length(tmp4)
                    tmp4{ii}= [tmp4{ii} ' (' num2str(length(find(strcmp(tmp4{ii}, tmp3)))) ')'];
                end
                set(int_list, 'string', tmp4, 'value', 1)
            end
            
            % figure out intersections
            if (curr_rad == 2)
                
                curr_ind_int= [];
                
                tmp1        = 1:decim(1, 1):pk{1}.num_trace_tot;
                tmp1(end)   = pk{1}.num_trace_tot;
                tmp2        = 1:decim(1, 2):pk{2}.num_trace_tot;
                tmp2(end)   = pk{2}.num_trace_tot;
                
                [~, ~, tmp3, tmp4] ...
                            = intersecti(pk{1}.x(tmp1), pk{1}.y(tmp1), pk{2}.x(tmp2), pk{2}.y(tmp2)); % intersection of two vectors using linear interpolation
                for ii = 1:length(tmp3)
                    if (tmp3(ii) == 1)
                        tmp5= 1:tmp1(2);
                    elseif (tmp3(ii) == length(tmp1))
                        tmp5= tmp1(end - 1):tmp1(end);
                    else
                        tmp5= tmp1(tmp3(ii) - 1):tmp1(tmp3(ii) + 1);
                    end
                    if (tmp4(ii) == 1)
                        tmp6= 1:tmp2(2);
                    elseif (tmp4(ii) == length(tmp2))
                        tmp6= tmp2(end - 1):tmp2(end);
                    else
                        tmp6= tmp2(tmp4(ii) - 1):tmp2(tmp4(ii) + 1);
                    end
                    [~, ~, tmp7, tmp8] ...
                            = intersecti(pk{1}.x(tmp5), pk{1}.y(tmp5), pk{2}.x(tmp6), pk{2}.y(tmp6));
                    if ~isempty(tmp3)
                        [tmp7, tmp8] ...
                            = deal((tmp7 + tmp5(1) - 1), (tmp8 + tmp6(1) - 1));
                        if ~isempty(find((diff(pk{1}.dist(tmp7)) < 1), 1))
                            tmp7 = tmp7(diff(pk{1}.dist(tmp7)) >= 1); % ensure density of consecutive intersections is not too great
                        end
                        if ~isempty(find((diff(pk{2}.dist(tmp8)) < 1), 1))
                            tmp8 = tmp8(diff(pk{2}.dist(tmp8)) >= 1);
                        end
                        curr_ind_int ...
                            = [curr_ind_int; tmp7 tmp8]; %#ok<AGROW>
                    end
                end
                
                tmp1        = [];
                
                % ensure density of consecutive intersections is not too great (1 per 2.5 km)
                if ~isempty(find(((abs(diff(pk{1}.dist(curr_ind_int(:, 1)))) < 2.5) & (abs(diff(pk{2}.dist(curr_ind_int(:, 2)))) < 2.5)), 1))
                    tmp1    = [tmp1 find(((abs(diff(pk{1}.dist(curr_ind_int(:, 1)))) < 2.5) & (abs(diff(pk{2}.dist(curr_ind_int(:, 2)))) < 2.5)))];
                end
                
                for ii = 1:size(curr_ind_int, 1)
                    % check for large jumps in transects (>5 times median distance change per trace
                    tmp2    = (curr_ind_int(ii, 1) - 1):(curr_ind_int(ii, 1) + 1);
                    tmp2    = tmp2((tmp2 > 0) & (tmp2 <= (pk{1}.num_trace_tot - 1)));
                    if any(diff(pk{1}.dist(tmp2)) > (5 * median(diff(pk{1}.dist))))
                        tmp1= [tmp1 ii]; %#ok<AGROW>
                    end
                    tmp3    = (curr_ind_int(ii, 2) - 1):(curr_ind_int(ii, 2) + 1);
                    tmp3    = tmp3((tmp3 > 0) & (tmp3 <= (pk{2}.num_trace_tot - 1)));
                    if any(diff(pk{2}.dist(tmp3)) > (5 * median(diff(pk{2}.dist))))
                        tmp1= [tmp1 ii]; %#ok<AGROW>
                    end
                    % transect intersection angle must be greater than 2 degrees
                    if (abs(atand(diff(pk{1}.y(tmp2([1 end]))) ./ diff(pk{1}.x(tmp2([1 end])))) - atand(diff(pk{2}.y(tmp3([1 end]))) ./ diff(pk{2}.x(tmp3([1 end]))))) < 2.5)
                        tmp1= [tmp1 ii]; %#ok<AGROW>
                    end
                end
                
                % trim intersections
                if ~isempty(unique(tmp1))
                    curr_ind_int ...
                            = curr_ind_int(setdiff(1:size(curr_ind_int, 1), unique(tmp1)), :); % the transect intersection where data actually exist
                end
                num_int     = size(curr_ind_int, 1);
                
                if isempty(curr_ind_int)
                    set(status_box(1), 'string', 'Selected merged picks file does not intersect primary transect.')
                    return
                end
                
                set(intnum_list, 'string', num2cell(1:num_int), 'value', 1)
            end
            
            set(data_list(curr_rad), 'string', pk{curr_rad}.file_block, 'value', 1)
            axes(ax(curr_ax))
            
            % display merged picks
            colors{curr_rad}= repmat(colors_def, ceil(pk{curr_rad}.num_layer / size(colors_def, 1)), 1); % extend predefined color pattern
            colors{curr_rad}= colors{curr_rad}(1:pk{curr_rad}.num_layer, :);
            [p_pk{1, curr_rad}, p_pk{2, curr_rad}] ...
                            = deal(zeros(1, pk{curr_rad}.num_layer));
            layer_str{curr_rad} ...
                            = num2cell(1:pk{curr_rad}.num_layer);
            for ii = 1:pk{curr_rad}.num_layer %#ok<*FXUP>
                if ~any(~isnan(pk{curr_rad}.elev_smooth(ii, ind_decim{1, curr_rad})))
                    p_pk{1, curr_rad}(ii) ...
                            = plot3(0, 0, 0, 'w.', 'markersize', 1, 'visible', 'off');
                    layer_str{curr_rad}{ii} ...
                            = [num2str(ii) ' H'];
                else
                    p_pk{1, curr_rad}(ii) ...
                            = plot3(pk{curr_rad}.x(ind_decim{1, curr_rad}(~isnan(pk{curr_rad}.elev_smooth(ii, ind_decim{1, curr_rad})))), ...
                                    pk{curr_rad}.y(ind_decim{1, curr_rad}(~isnan(pk{curr_rad}.elev_smooth(ii, ind_decim{1, curr_rad})))), ...
                                    pk{curr_rad}.elev_smooth(ii, ind_decim{1, curr_rad}(~isnan(pk{curr_rad}.elev_smooth(ii, ind_decim{1, curr_rad})))), ...
                                    '.', 'color', colors{curr_rad}(ii, :), 'markersize', 12, 'visible', 'off');
                end
            end
            
            set(layer_list(:, curr_rad), 'string', layer_str{curr_rad}, 'value', 1)
            axes(ax(curr_ax))
            
            % calculate bed depth (for flattened projection)
            if (surf_avail(curr_rad) && bed_avail(curr_rad))
                depth_bed{curr_rad} ...
                            = pk{curr_rad}.elev_surf(ind_decim{2, curr_rad}) - pk{curr_rad}.elev_bed(ind_decim{2, curr_rad});
            end
            
            % display surface and bed
            if surf_avail(curr_rad)
                p_surf(curr_gui, curr_rad) ...
                            = plot3(pk{curr_rad}.x(ind_decim{1, curr_rad}(~isnan(pk{curr_rad}.elev_surf(ind_decim{1, curr_rad})))), ...
                                    pk{curr_rad}.y(ind_decim{1, curr_rad}(~isnan(pk{curr_rad}.elev_surf(ind_decim{1, curr_rad})))), ...
                                    pk{curr_rad}.elev_surf(ind_decim{1, curr_rad}(~isnan(pk{curr_rad}.elev_surf(ind_decim{1, curr_rad})))), 'k.', 'markersize', 12, 'visible', 'off');
            end
            if bed_avail(curr_rad)
                p_bed(curr_gui, curr_rad) ...
                            = plot3(pk{curr_rad}.x(ind_decim{1, curr_rad}(~isnan(pk{curr_rad}.elev_bed(ind_decim{1, curr_rad})))), ...
                                    pk{curr_rad}.y(ind_decim{1, curr_rad}(~isnan(pk{curr_rad}.elev_bed(ind_decim{1, curr_rad})))), ...
                                    pk{curr_rad}.elev_bed(ind_decim{1, curr_rad}(~isnan(pk{curr_rad}.elev_bed(ind_decim{1, curr_rad})))), 'k.', 'markersize', 12, 'visible', 'off');
            end
            
            % display picks, surface and bed in 2D GUI
            axes(ax(curr_ax + curr_rad))
            for ii = 1:pk{curr_rad}.num_layer
                if ~any(~isnan(pk{curr_rad}.elev_smooth(ii, ind_decim{2, curr_rad})))
                    p_pk{2, curr_rad}(ii) ...
                            = plot(0, 0, 'w.', 'markersize', 1, 'visible', 'off');
                else
                    p_pk{2, curr_rad}(ii) ...
                            = plot(pk{curr_rad}.dist_lin(ind_decim{2, curr_rad}(~isnan(pk{curr_rad}.elev_smooth(ii, ind_decim{2, curr_rad})))), ...
                                   pk{curr_rad}.elev_smooth(ii, ind_decim{2, curr_rad}(~isnan(pk{curr_rad}.elev_smooth(ii, ind_decim{2, curr_rad})))), ...
                                   '.', 'color', colors{curr_rad}(ii, :), 'markersize', 12, 'visible', 'off');
                end
            end
            if surf_avail(curr_rad)
                p_surf(2, curr_rad) ...
                            = plot(pk{curr_rad}.dist_lin(ind_decim{2, curr_rad}(~isnan(pk{curr_rad}.elev_surf(ind_decim{2, curr_rad})))), ...
                                   pk{curr_rad}.elev_surf(ind_decim{2, curr_rad}(~isnan(pk{curr_rad}.elev_surf(ind_decim{2, curr_rad})))), 'g--', 'linewidth', 2, 'visible', 'off');
                if any(isnan(pk{curr_rad}.elev_surf(ind_decim{2, curr_rad}(~isnan(pk{curr_rad}.elev_surf(ind_decim{2, curr_rad}))))))
                    set(p_surf(2, curr_rad), 'marker', '.', 'linestyle', 'none', 'markersize', 12)
                end
            end
            if bed_avail(curr_rad)
                p_bed(2, curr_rad) ...
                            = plot(pk{curr_rad}.dist_lin(ind_decim{2, curr_rad}(~isnan(pk{curr_rad}.elev_bed(ind_decim{2, curr_rad})))), ...
                                   pk{curr_rad}.elev_bed(ind_decim{2, curr_rad}(~isnan(pk{curr_rad}.elev_bed(ind_decim{2, curr_rad})))), 'g--', 'linewidth', 2, 'visible', 'off');
                if any(isnan(pk{curr_rad}.elev_bed(ind_decim{2, curr_rad}(~isnan(pk{curr_rad}.elev_bed(ind_decim{2, curr_rad}))))))
                    set(p_bed(2, curr_rad), 'marker', '.', 'linestyle', 'none', 'markersize', 12)
                end
            end
            axes(ax(curr_ax))
            
            % get new min/max limits for 3D display
            if (curr_rad == 2)
                [tmp1, tmp2, tmp3, tmp4, tmp5] ...
                            = deal([pk{1}.x pk{2}.x], [pk{1}.y pk{2}.y], [pk{1}.elev_surf pk{2}.elev_surf], [pk{1}.elev_bed pk{2}.elev_bed], [pk{1}.elev_smooth(:); pk{2}.elev_smooth(:)]);
                [x_min_ref, x_max_ref, x_min, x_max] ...
                            = deal(min(tmp1), max(tmp1), min(tmp1), max(tmp1));
                [y_min_ref, y_max_ref, y_min, y_max] ...
                            = deal(min(tmp2), max(tmp2), min(tmp2), max(tmp2));
                [dist_min_ref(curr_rad), dist_max_ref(curr_rad), dist_min(curr_rad), dist_max(curr_rad)] ...
                            = deal(min(pk{curr_rad}.dist_lin), max(pk{curr_rad}.dist_lin), min(pk{curr_rad}.dist_lin), max(pk{curr_rad}.dist_lin));
                if surf_avail(curr_rad)
                    [elev_max_ref, elev_max(1:2)] ...
                            = deal(max(tmp3) + (0.1 * (max(tmp3) - min(tmp3))));
                else
                    [elev_max_ref, elev_max(1:2)] ...
                            = deal(max(tmp5) + (0.1 * (max(tmp5) - min(tmp5))));
                end
                if bed_avail(curr_rad)
                    [elev_min_ref, elev_min(1:2)] ...
                            = deal(min(tmp4) - (0.1 * (max(tmp4) - min(tmp4))));
                else
                    [elev_min_ref, elev_min(1:2)] ...
                            = deal(min(tmp5) - (0.1 * (max(tmp5) - min(tmp5))));
                end
            else
                [x_min_ref, x_max_ref, x_min, x_max] ...
                            = deal(min(pk{curr_rad}.x), max(pk{curr_rad}.x), min(pk{curr_rad}.x), max(pk{curr_rad}.x));
                [y_min_ref, y_max_ref, y_min, y_max] ...
                            = deal(min(pk{curr_rad}.y), max(pk{curr_rad}.y), min(pk{curr_rad}.y), max(pk{curr_rad}.y));
                [dist_min_ref(curr_rad), dist_max_ref(curr_rad), dist_min(curr_rad), dist_max(curr_rad)] ...
                            = deal(min(pk{curr_rad}.dist_lin), max(pk{curr_rad}.dist_lin), min(pk{curr_rad}.dist_lin), max(pk{curr_rad}.dist_lin));
                if surf_avail(curr_rad)
                    [elev_max_ref, elev_max(1:2)] ...
                            = deal(max(pk{curr_rad}.elev_surf) + (0.1 * (max(pk{curr_rad}.elev_surf) - min(pk{curr_rad}.elev_surf))));
                else
                    [elev_max_ref, elev_max(1:2)] ...
                            = deal(max(pk{curr_rad}.elev_smooth(:)) + (0.1 * (max(pk{curr_rad}.elev_smooth(:)) - min(pk{curr_rad}.elev_smooth(:)))));
                end
                if bed_avail(curr_rad)
                    [elev_min_ref, elev_min(1:2)] ...
                            = deal(min(pk{curr_rad}.elev_bed) - (0.1 * (max(pk{curr_rad}.elev_bed) - min(pk{curr_rad}.elev_bed))));
                else
                    [elev_min_ref, elev_min(1:2)] ...
                            = deal(min(pk{curr_rad}.elev_smooth(:)) - (0.1 * (max(pk{curr_rad}.elev_smooth(:)) - min(pk{curr_rad}.elev_smooth(:)))));
                end
            end
            set(x_min_slide, 'min', x_min_ref, 'max', x_max_ref, 'value', x_min_ref)
            set(x_max_slide, 'min', x_min_ref, 'max', x_max_ref, 'value', x_max_ref)
            set(y_min_slide, 'min', y_min_ref, 'max', y_max_ref, 'value', y_min_ref)
            set(y_max_slide, 'min', y_min_ref, 'max', y_max_ref, 'value', y_max_ref)
            set(dist_min_slide(curr_rad), 'min', dist_min_ref(curr_rad), 'max', dist_max_ref(curr_rad), 'value', dist_min_ref(curr_rad))
            set(dist_max_slide(curr_rad), 'min', dist_min_ref(curr_rad), 'max', dist_max_ref(curr_rad), 'value', dist_max_ref(curr_rad))
            set(z_min_slide, 'min', elev_min_ref, 'max', elev_max_ref, 'value', elev_min_ref)
            set(z_max_slide, 'min', elev_min_ref, 'max', elev_max_ref, 'value', elev_max_ref)
            set(x_min_edit, 'string', sprintf('%4.1f', x_min_ref))
            set(x_max_edit, 'string', sprintf('%4.1f', x_max_ref))
            set(y_min_edit, 'string', sprintf('%4.1f', y_min_ref))
            set(y_max_edit, 'string', sprintf('%4.1f', y_max_ref))
            set(dist_min_edit(curr_rad), 'string', sprintf('%4.1f', dist_min_ref(curr_rad)))
            set(dist_max_edit(curr_rad), 'string', sprintf('%4.1f', dist_max_ref(curr_rad)))
            set(z_min_edit, 'string', sprintf('%4.0f', elev_min_ref))
            set(z_max_edit, 'string', sprintf('%4.0f', elev_max_ref))
            update_x_range
            update_y_range
            update_z_range
            
            pk_done(curr_rad) ...
                            = true;
            set(pk_check(:, curr_rad), 'value', 1)
            show_pk
            
            % intersection plots
            if (curr_rad == 2)
                for ii = 1:2
                    for jj = 1:3
                        p_int1{ii, jj} = zeros(1, num_int);
                    end
                end
                axes(ax(1))
                for ii = 1:num_int
                    p_int1{1, 1}(ii) = plot3(repmat(pk{1}.x(curr_ind_int(ii, 1)), 1, 2), repmat(pk{1}.y(curr_ind_int(ii, 1)), 1, 2), [elev_min_ref elev_max_ref], 'm--', 'linewidth', 2, 'visible', 'off');
                end
                axes(ax(2))
                for ii = 1:num_int
                    p_int1{1, 2}(ii) = plot(repmat(pk{1}.dist_lin(curr_ind_int(ii, 1)), 1, 2), [elev_min_ref elev_max_ref], 'm--', 'linewidth', 2, 'visible', 'off');
                end
                if (any(p_int2{1, 1}) && any(ishandle(p_int2{1, 1})))
                    delete(p_int2{1, 1}(logical(p_int2{1, 1}) & ishandle(p_int2{1, 1})))
                end
                if (any(p_int2{2, 1}) && any(ishandle(p_int2{2, 1})))
                    delete(p_int2{2, 1}(logical(p_int2{2, 1}) & ishandle(p_int2{2, 1})))
                end
                p_int2{1, 1} = zeros(1, pk{2}.num_layer);
                for ii = 1:pk{2}.num_layer
                    p_int2{1, 1}(ii) = plot(pk{1}.dist_lin(curr_ind_int(:, 1)), pk{2}.elev_smooth(ii, curr_ind_int(:, 2)), 'ko', 'markersize', 8, 'markerfacecolor', colors{2}(ii, :), 'visible', 'off');
                end
                axes(ax(3))
                for ii = 1:num_int
                    p_int1{1, 3}(ii) = plot(repmat(pk{2}.dist_lin(curr_ind_int(ii, 2)), 1, 2), [elev_min_ref elev_max_ref], 'm--', 'linewidth', 2, 'visible', 'off');
                end
                p_int2{1, 2} = zeros(1, pk{1}.num_layer);
                for ii = 1:pk{1}.num_layer
                    p_int2{1, 2}(ii) = plot(pk{2}.dist_lin(curr_ind_int(:, 2)), pk{1}.elev_smooth(ii, curr_ind_int(:, 1)), 'ko', 'markersize', 8, 'markerfacecolor', colors{1}(ii, :), 'visible', 'off');
                end
                if ~isempty(pk{2}.ind_layer)
                    for ii = 1:size(pk{2}.ind_layer, 1)
                        if ((pk{2}.ind_layer(ii, 2) == curr_year(1)) && (pk{2}.ind_layer(ii, 3) == curr_trans(1))) % match associated with current transect
                            set(p_int2{1, 1}(pk{2}.ind_layer(ii, 1)), 'marker', '^', 'markerfacecolor', colors{1}(pk{2}.ind_layer(ii, 4), :))
                            set(p_int2{1, 2}(pk{2}.ind_layer(ii, 4)), 'marker', '^', 'markerfacecolor', colors{1}(pk{2}.ind_layer(ii, 4), :))
                            set(p_pk{1, 2}(pk{2}.ind_layer(ii, 1)), 'color', colors{1}(pk{2}.ind_layer(ii, 4), :))
                            set(p_pk{2, 2}(pk{2}.ind_layer(ii, 1)), 'color', colors{1}(pk{2}.ind_layer(ii, 4), :))
                        else
                            set(p_int2{1, 1}(pk{2}.ind_layer(ii, 1)), 'marker', 's')
                            set(p_int2{1, 2}(pk{2}.ind_layer(ii, 4)), 'marker', 's')
                        end
                    end
                end
                set(int_check, 'value', 1)
                show_int1
                show_int2
                show_int3
            end
            
            % show gui2 display items
            switch curr_rad
                case 1
                    [curr_gui, curr_ax] ...
                            = deal(2);
                case 2
                    [curr_gui, curr_ax] ...
                            = deal(2, 3);
            end
            axes(ax(curr_ax))
            update_z_range
            update_dist_range
            show_pk
            [curr_gui, curr_ax] ...
                            = deal(1);
            axes(ax(curr_ax))
            
            if (core_done && all(pk_done))
                load_core_breakout
            end
            
            if ((curr_rad == 1) && strcmp(disp_type{2}, disp_type{3}))
                linkaxes(ax(2:3), 'y')
                h_link1   = linkprop(z_min_slide(2:3), {'value' 'min' 'max'});
                h_link2   = linkprop(z_max_slide(2:3), {'value' 'min' 'max'});
                h_link3   = linkprop(z_min_edit(2:3), 'string');
                h_link4   = linkprop(z_max_edit(2:3), 'string');
            end
            
            if (curr_rad == 2)
                reset_xz1
            end
            
            set(status_box(1), 'string', ['Loaded ' file_pk_short{curr_rad} '.'])
            
        else
            set(status_box(1), 'string', 'No picks loaded.')
        end
    end

%% Load radar data

    function load_data1(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 2, 1, 2);
        load_data
    end

    function load_data2(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 3, 2, 1);
        load_data
    end

    function load_data(source, eventdata) %#ok<*INUSD>
        
        if ~pk_done(curr_rad)
            set(status_box(2), 'string', 'Load picks before data.')
            return
        end
        
        if ispc
            if (~isempty(path_pk{curr_rad}) && exist([path_pk{curr_rad} '..\block\' file_pk_short{curr_rad}], 'dir'))
                path_data{curr_rad} = [path_pk{curr_rad}(1:strfind(path_pk{curr_rad}, '\merge')) 'block\' file_pk_short{curr_rad} '\'];
            end
        else
            if (~isempty(path_pk{curr_rad}) && exist([path_pk{curr_rad} '../block/' file_pk_short{curr_rad}], 'dir'))
                path_data{curr_rad} = [path_pk{curr_rad}(1:strfind(path_pk{curr_rad}, '/merge')) 'block/' file_pk_short{curr_rad} '/'];
            end
        end
        
        if (~isempty(path_data{curr_rad}) && exist([path_data{curr_rad} pk{curr_rad}.file_block{1} '.mat'], 'file'))
            file_data{curr_rad} = pk{curr_rad}.file_block;
            for ii = 1:length(file_data{curr_rad})
                file_data{curr_rad}{ii} ...
                                = [file_data{curr_rad}{ii} '.mat'];
            end
        elseif ~isempty(path_data{curr_rad}) % Dialog box to choose radar data file to load
            [file_data{curr_rad}, path_data{curr_rad}] ...
                                = uigetfile('*.mat', 'Load radar data:', path_data{curr_rad}, 'multiselect', 'on');
        elseif ~isempty(path_data{curr_rad_alt})
            [file_data{curr_rad}, path_data{curr_rad}] ...
                                = uigetfile('*.mat', 'Load radar data:', path_data{curr_rad_alt}, 'multiselect', 'on');
        elseif ~isempty(path_pk{curr_rad})
            [file_data{curr_rad}, path_data{curr_rad}] ...
                                = uigetfile('*.mat', 'Load radar data:', path_pk{curr_rad}, 'multiselect', 'on');
        elseif ~isempty(path_pk{curr_rad_alt})
            [file_data{curr_rad}, path_data{curr_rad}] ...
                                = uigetfile('*.mat', 'Load radar data:', path_pk{curr_rad_alt}, 'multiselect', 'on');
        else
            [file_data{curr_rad}, path_data{curr_rad}] ...
                                = uigetfile('*.mat', 'Load radar data:', 'multiselect', 'on');
        end
        
        if isnumeric(file_data{curr_rad})
            [file_data{curr_rad}, path_data{curr_rad}] ...
                                = deal('');
        end
        
        if ~isempty(file_data{curr_rad})
            
            pause(0.1)
            
            num_data(curr_rad) ...
                            = length(file_data{curr_rad});
            if (num_data(curr_rad) ~= length(pk{curr_rad}.file_block))
                set(status_box(2), 'string', ['Number of data blocks selected (' num2str(num_data(curr_rad)) ') does not match number of blocks in intersection list (' num2str(length(pk{curr_rad}.file_block)) ').'])
                return
            end
            if ~strcmp(file_data{curr_rad}, [pk{curr_rad}.file_block{1} '.mat'])
                set(status_box(2), 'string', 'Correct intersecting blocks may not have been selected. Try again.')
                return
            end
            
            % get rid of all data handle if somehow present
            if data_done(curr_rad)
                data_done(curr_rad) ...
                                = false;
            end
            
            load_data_breakout
            
        else
            set(status_box(2), 'string', 'No radar data loaded.')
        end
    end

    function load_data_breakout(source, eventdata)
        
        axes(ax(curr_ax))
        
        % attempt to load data
        set(status_box(2), 'string', 'Loading radar data...')
        pause(0.1)
        
        for ii = 1:num_data(curr_rad)
            
            set(status_box(2), 'string', ['Loading ' file_data{curr_rad}{ii}(1:(end - 4)) ' (' num2str(ii) ' / ' num2str(num_data(curr_rad)) ')...'])
            pause(0.1)
            tmp1            = load([path_data{curr_rad} file_data{curr_rad}{ii}]);
            try
                tmp1        = tmp1.block;
            catch %#ok<*CTCH>
                set(status_box(2), 'string', [file_data{curr_rad}{ii} ' does not contain a block structure. Try again.'])
                return
            end
            
            if (ii == 1) % start decimation
                
                [dt(curr_rad), twtt{curr_rad}, num_sample(curr_rad)] ...
                            = deal(tmp1.dt, tmp1.twtt, tmp1.num_sample);
                amp_mean{curr_rad} ...
                            = NaN(num_sample(curr_rad), num_decim(2, curr_rad));
                tmp2        = (1 + ceil(decim(2, curr_rad) / 2)):decim(2, curr_rad):(pk{curr_rad}.num_trace(ii) - ceil(decim(2, curr_rad) / 2));
                tmp3        = floor(decim(2, curr_rad) / 2);
                for jj = 1:length(tmp2)
                    amp_mean{curr_rad}(:, jj) ...
                            = nanmean(tmp1.amp(:, (tmp2(jj) - tmp3):(tmp2(jj) + tmp3)), 2);
                end
                
            else % middle/end
                
                tmp2        = repmat((pk{curr_rad}.ind_trace_start(ii) + pk{curr_rad}.ind_overlap(ii, 1)), 1, 2);
                tmp2(2)     = tmp2(2) + pk{curr_rad}.num_trace(ii) - pk{curr_rad}.ind_overlap(ii, 1) - 1;
                if (ii < num_data(curr_rad))
                    tmp3    = (find((tmp2(1) > ind_decim{2, curr_rad}), 1, 'last') + 1):find((tmp2(2) <= ind_decim{2, curr_rad}), 1);
                else
                    tmp3    = (find((tmp2(1) > ind_decim{2, curr_rad}), 1, 'last') + 1):num_decim(2, curr_rad);
                end
                if (ii == 2)
                    tmp3    = [(tmp3(1) - 2) (tmp3(1) - 1) tmp3]; %#ok<AGROW> % correct weirdness for first block
                end
                tmp2        = ind_decim{2, curr_rad}(tmp3) - pk{curr_rad}.ind_trace_start(ii) + 1;
                tmp4        = tmp2 - floor(decim(2, curr_rad) / 2);
                tmp5        = tmp2 + floor(decim(2, curr_rad) / 2);
                tmp4(tmp4 < 1) ...
                            = 1;
                tmp5(tmp5 > pk{curr_rad}.num_trace(ii)) ...
                            = pk{curr_rad}.num_trace(ii);
                for jj = 1:length(tmp3)
                    amp_mean{curr_rad}(:, tmp3(jj)) ...
                            = nanmean(tmp1.amp(:, tmp4(jj):tmp5(jj)), 2);
                end
                
            end
            
            tmp1            = 0;
            
        end
        
        % convert to dB
        amp_mean{curr_rad}  = 10 .* log10(abs(amp_mean{curr_rad}));
        num_sample(curr_rad)= size(amp_mean{curr_rad}, 1);
        depth{curr_rad}     = ((speed_ice / 2) .* twtt{curr_rad}); % simple monotonically increasing depth vector
        tmp1                = NaN(size(amp_mean{curr_rad}), 'single');
        if surf_avail(curr_rad)
            tmp2            = interp1(twtt{curr_rad}, 1:num_sample(curr_rad), pk{curr_rad}.twtt_surf(ind_decim{2, curr_rad}), 'nearest', 'extrap'); % surface traveltime indices
            if any(isnan(tmp2))
                tmp2(isnan(tmp2)) ...
                            = round(interp1(find(~isnan(tmp2)), tmp2(~isnan(tmp2)), find(isnan(tmp2)), 'linear', 'extrap'));
            end
        else
            tmp2            = ones(1, num_decim(2, curr_rad));
        end
        tmp3                = pk{curr_rad}.elev_surf(ind_decim{2, curr_rad});
        if any(isnan(tmp3))
            tmp3(isnan(tmp3)) ...
                            = interp1(find(~isnan(tmp3)), tmp3(~isnan(tmp3)), find(isnan(tmp3)), 'linear', 'extrap');
        end
        for ii = 1:num_decim(2, curr_rad)
            tmp1(1:(num_sample(curr_rad) - tmp2(ii) + 1), ii) ...
                            = amp_mean{curr_rad}(tmp2(ii):num_sample(curr_rad), ii); % shift data up to surface
        end
        [amp_mean{curr_rad}, ind_corr{curr_rad}] ...
                            = topocorr(tmp1, depth{curr_rad}, tmp3); % topographically correct data
        tmp1                = 0;
        amp_mean{curr_rad}  = flipud(amp_mean{curr_rad}); % flip for axes
        ind_corr{curr_rad}  = max(ind_corr{curr_rad}) - ind_corr{curr_rad} + 1;
        depth{curr_rad}     = (speed_ice / 2) .* (0:dt(curr_rad):((num_sample(curr_rad) - 1) * dt(curr_rad)))'; % simple monotonically increasing depth vector
        elev{curr_rad}      = flipud(max(pk{curr_rad}.elev_surf(ind_decim{2, curr_rad})) - depth{curr_rad}); % elevation vector
        depth{curr_rad}     = depth{curr_rad}(1:(num_sample(curr_rad) - max(ind_corr{curr_rad})));
        
        if flat_done(curr_rad)
            
            set(status_box(2), 'string', 'Polynomials available so now flattening...')
            pause(0.1)
            
            tic
            
            % fix polynomials to current decimation vector
            if (num_decim(2, curr_rad) ~= size(pk{curr_rad}.poly_flat_merge, 2))
                pk{curr_rad}.poly_flat_merge ...
                            = interp2(pk{curr_rad}.poly_flat_merge, linspace(1, num_decim(2, curr_rad), num_decim(2, curr_rad)), (1:size(pk{curr_rad}.poly_flat_merge, 1))');
            end
            
            depth_mat{curr_rad} ...
                            = single(depth{curr_rad}(:, ones(1, num_decim(2, curr_rad)))); % depth matrix
            depth_flat{curr_rad} ...
                            = ((depth_mat{curr_rad} .^ 3) .* pk{curr_rad}.poly_flat_merge(ones((num_sample(curr_rad) - max(ind_corr{curr_rad})), 1), :)) + ...
                              ((depth_mat{curr_rad} .^ 2) .* pk{curr_rad}.poly_flat_merge((2 .* ones((num_sample(curr_rad) - max(ind_corr{curr_rad})), 1)), :)) + ...
                              (depth_mat{curr_rad} .* (pk{curr_rad}.poly_flat_merge((3 .* ones((num_sample(curr_rad) - max(ind_corr{curr_rad})), 1)), :))) + ...
                              pk{curr_rad}.poly_flat_merge((4 .* ones((num_sample(curr_rad) - max(ind_corr{curr_rad})), 1)), :);
            depth_mat{curr_rad} ...
                            = 0;
            depth_flat{curr_rad}(depth_flat{curr_rad} < 0) ...
                            = 0;
            depth_flat{curr_rad}(depth_flat{curr_rad} > depth{curr_rad}(end)) ...
                            = depth{curr_rad}(end);
            set(status_box(2), 'string', ['Calculated remapping in ' num2str(toc, '%.0f') ' s...'])
            pause(0.1)
            
            tic
            
            % flattened radargram based on the polyfits
            tmp1            = flipud(amp_mean{curr_rad});
            [amp_flat{curr_rad}, tmp2] ...
                            = deal(NaN((num_sample(curr_rad) - max(ind_corr{curr_rad})), num_decim(2, curr_rad), 'single'));
            for ii = 1:num_decim(2, curr_rad)
                tmp2(:, ii) = tmp1(ind_corr{curr_rad}(ii):(end - (max(ind_corr{curr_rad}) - ind_corr{curr_rad}(ii)) - 1), ii); % grab starting at surface
            end
            tmp1            = tmp2;
            tmp2            = find(sum(~isnan(depth_flat{curr_rad})));
            if parallel_check
                tmp1        = tmp1(:, tmp2);
                tmp3        = depth_flat{curr_rad}(:, tmp2);
                tmp4        = amp_flat{curr_rad}(:, tmp2);
                pctRunOnAll warning('off', 'MATLAB:interp1:NaNinY')
                parfor ii = 1:length(tmp2)
                    tmp4(:, ii) = interp1(depth{curr_rad}, tmp1(:, ii), tmp3(:, ii)); %#ok<PFBNS>
                end
                pctRunOnAll warning('on', 'MATLAB:interp1:NaNinY')
                amp_flat{curr_rad}(:, tmp2) ...
                            = tmp4;
            else
                warning('off', 'MATLAB:interp1:NaNinY')
                for ii = 1:length(tmp2)
                    amp_flat{curr_rad}(:, tmp2(ii)) ...
                            = interp1(depth{curr_rad}, tmp1(:, tmp2(ii)), depth_flat{curr_rad}(:, tmp2(ii)), 'linear');
                end
                warning('on', 'MATLAB:interp1:NaNinY')
            end
            tmp1            = 0;
            
            set(status_box(2), 'string', ['Flattened amplitude in ' num2str(toc, '%.0f') ' s...'])
            pause(0.1)
            
            tic
            
            % flatten bed pick
            if bed_avail(curr_rad)
                depth_bed_flat{curr_rad} ...
                            = NaN(1, num_decim(2, curr_rad));
                if parallel_check
                    tmp3    = depth_flat{curr_rad}(:, tmp2);
                    tmp4    = depth_bed_flat{curr_rad}(tmp2);
                    tmp5    = depth_bed{curr_rad}(tmp2);
                    parfor ii = 1:length(tmp2)
                        [~,  tmp1] ...
                            = unique(tmp3(:, ii), 'last');
                        tmp1= intersect((1 + find(diff(tmp3(:, ii)) > 0)), tmp1);
                        if (length(tmp1) > 1)
                            tmp4(ii) ...
                            = interp1(tmp3(tmp1, ii), depth{curr_rad}(tmp1), tmp5(ii), 'nearest', 'extrap'); %#ok<PFBNS>
                        end
                    end
                    depth_bed_flat{curr_rad}(tmp2) ...
                            = tmp4;
                else
                    for ii = 1:length(tmp2)
                        [~, tmp1] ...
                            = unique(depth_flat{curr_rad}(:, tmp2(ii)), 'last');
                        tmp1= intersect((1 + find(diff(depth_flat{curr_rad}(:, ii)) > 0)), tmp1);
                        if (length(tmp1) > 1)
                            depth_bed_flat{curr_rad}(tmp2(ii)) ...
                            = interp1(depth_flat{curr_rad}(tmp1, tmp2(ii)), depth{curr_rad}(tmp1), depth_bed{curr_rad}(tmp2(ii)), 'nearest', 'extrap');
                        end
                    end
                end
                depth_bed_flat{curr_rad}((depth_bed_flat{curr_rad} < 0) | (depth_bed_flat{curr_rad} > depth{curr_rad}(end))) ...
                            = NaN;
            end
            
            % flatten merged layer depths
            warning('off', 'MATLAB:interp1:NaNinY')
            depth_layer_flat{curr_rad} ...
                            = NaN(pk{curr_rad}.num_layer, num_decim(2, curr_rad));
            for ii = 1:length(tmp2)
                [~, tmp1]   = unique(depth_flat{curr_rad}(:, tmp2(ii)), 'last');
                tmp1        = intersect((1 + find(diff(depth_flat{curr_rad}(:, tmp2(ii))) > 0)), tmp1);
                depth_layer_flat{curr_rad}(~isnan(pk{curr_rad}.depth(:, ind_decim{2, curr_rad}(tmp2(ii)))), tmp2(ii)) ...
                            = interp1(depth_flat{curr_rad}(tmp1, tmp2(ii)), depth{curr_rad}(tmp1), pk{curr_rad}.depth(~isnan(pk{curr_rad}.depth(:, ind_decim{2, curr_rad}(tmp2(ii)))), ...
                                      ind_decim{2, curr_rad}(tmp2(ii))), 'nearest', 'extrap');
            end
            warning('on', 'MATLAB:interp1:NaNinY')
            
            % plot flat layers
            axes(ax(curr_ax))
            p_pkflat{curr_rad} ...
                            = zeros(1, pk{curr_rad}.num_layer);
            for ii = 1:pk{curr_rad}.num_layer
                if ~isempty(find(~isnan(depth_layer_flat{curr_rad}(ii, :)), 1))
                    p_pkflat{curr_rad}(ii) ...
                            = plot(pk{curr_rad}.dist_lin(ind_decim{2, curr_rad}(~isnan(depth_layer_flat{curr_rad}(ii, :)))), depth_layer_flat{curr_rad}(ii, ~isnan(depth_layer_flat{curr_rad}(ii, :))), '.', ...
                                   'markersize', 12, 'color', colors{curr_rad}(ii, :), 'visible', 'off');
                else
                    p_pkflat{curr_rad}(ii) ...
                            = plot(0, 0, '.', 'markersize', 1, 'color', colors{curr_rad}(ii, :), 'visible', 'off');
                end
            end

            % assign traveltime and distance reference values/sliders based on data
            [elev_min_ref, db_min_ref(curr_rad), elev_max_ref, db_max_ref(curr_rad), elev_min(2), db_min(curr_rad), elev_max(2), db_max(curr_rad), depth_min_ref, depth_max_ref, depth_min, depth_max] ...
                            = deal(min(elev{curr_rad}), min(amp_mean{curr_rad}(~isinf(amp_mean{curr_rad}(:)))), max(elev{curr_rad}), max(amp_mean{curr_rad}(~isinf(amp_mean{curr_rad}(:)))), min(elev{curr_rad}), ...
                                   min(amp_mean{curr_rad}(~isinf(amp_mean{curr_rad}(:)))), max(elev{curr_rad}), max(amp_mean{curr_rad}(~isinf(amp_mean{curr_rad}(:)))), min(depth{curr_rad}), max(depth{curr_rad}), ...
                                   min(depth{curr_rad}), max(depth{curr_rad}));
            set(cb_min_slide(curr_rad), 'min', db_min_ref(curr_rad), 'max', db_max_ref(curr_rad), 'value', db_min_ref(curr_rad))
            set(cb_max_slide(curr_rad), 'min', db_min_ref(curr_rad), 'max', db_max_ref(curr_rad), 'value', db_max_ref(curr_rad))
            set(cb_min_edit(curr_rad), 'string', sprintf('%3.0f', db_min_ref(curr_rad)))
            set(cb_max_edit(curr_rad), 'string', sprintf('%3.0f', db_max_ref(curr_rad)))
            update_z_range
            
            if (curr_rad == 2)
                axes(ax(1))
                axes(ax(2))
                for ii = 1:num_int
                    p_int1{2, 2}(ii) = plot(repmat(pk{1}.dist_lin(curr_ind_int(ii, 1)), 1, 2), [depth_min_ref depth_max_ref], 'm--', 'linewidth', 2, 'visible', 'off');
                end
                axes(ax(3))
                for ii = 1:num_int
                    p_int1{2, 3}(ii) = plot(repmat(pk{2}.dist_lin(curr_ind_int(ii, 2)), 1, 2), [depth_min_ref depth_max_ref], 'm--', 'linewidth', 2, 'visible', 'off');
                end
                if (any(p_int2{2, 1}) && any(ishandle(p_int2{2, 1})))
                    delete(p_int2{2, 1}(logical(p_int2{2, 1}) & ishandle(p_int2{2, 1})))
                end
                p_int2{2, 1} = zeros(1, pk{2}.num_layer);
                for ii = 1:pk{2}.num_layer
                    p_int2{2, 1}(ii) = plot(pk{1}.dist_lin(curr_ind_int(:, 1)), pk{2}.depth_smooth(ii, curr_ind_int(:, 2)), 'ko', 'markersize', 8, 'markerfacecolor', colors{2}(ii, :), 'visible', 'off');
                end
                p_int2{2, 2} = zeros(1, pk{1}.num_layer);
                for ii = 1:pk{1}.num_layer
                    p_int2{2, 2}(ii) = plot(pk{2}.dist_lin(curr_ind_int(:, 2)), pk{1}.depth_smooth(ii, curr_ind_int(:, 1)), 'ko', 'markersize', 8, 'markerfacecolor', colors{1}(ii, :), 'visible', 'off');
                end
                if ~isempty(pk{2}.ind_layer)
                    for ii = 1:size(pk{2}.ind_layer, 1)
                        if ((pk{2}.ind_layer(ii, 2) == curr_year(1)) && (pk{2}.ind_layer(ii, 3) == curr_trans(1))) % match associated with current transect
                            set(p_int2{2, 1}(pk{2}.ind_layer(ii, 1)), 'marker', '^', 'markerfacecolor', colors{1}(pk{2}.ind_layer(ii, 4), :))
                            set(p_int2{2, 2}(pk{2}.ind_layer(ii, 4)), 'marker', '^', 'markerfacecolor', colors{1}(pk{2}.ind_layer(ii, 4), :))
                            set(p_pkflat{2}(pk{2}.ind_layer(ii, 1)), 'color', colors{1}(pk{2}.ind_layer(ii, 4), :))
                        else
                            set(p_int2{2, 1}(pk{2}.ind_layer(ii, 1)), 'marker', 's')
                            set(p_int2{2, 2}(pk{2}.ind_layer(ii, 4)), 'marker', 's')
                        end
                    end
                end
                set(int_check, 'value', 1)
                show_int1
                show_int2
                show_int3
            end            
            
        end
                
        if all(pk_done)
            for ii = 1:3
                if (any(p_int1{1, ii}) && any(ishandle(p_int1{1, ii})))
                    if (ii == 1)
                        set(p_int1{1, ii}(logical(p_int1{1, ii}) & ishandle(p_int1{1, ii})), 'zdata', [elev_min_ref elev_max_ref])
                    else
                        set(p_int1{1, ii}(logical(p_int1{1, ii}) & ishandle(p_int1{1, ii})), 'ydata', [elev_min_ref elev_max_ref])
                    end
                end
            end
            for ii = 2:3
                if (any(p_int1{2, ii}) && any(ishandle(p_int1{2, ii})))
                    set(p_int1{2, ii}(logical(p_int1{2, ii}) & ishandle(p_int1{2, ii})), 'ydata', [depth_min_ref depth_max_ref])
                end
            end
        end
        
        if core_done
            for ii = 1:2
                if (any(p_core{ii, curr_rad}) && any(ishandle(p_core{ii, curr_rad})))
                    if (ii == 1)
                        set(p_core{ii, curr_rad}(logical(p_core{ii, curr_rad}) & ishandle(p_core{ii, curr_rad})), 'zdata', [elev_min_ref elev_max_ref])
                    else
                        set(p_core{ii, curr_rad}(logical(p_core{ii, curr_rad}) & ishandle(p_core{ii, curr_rad})), 'ydata', [elev_min_ref elev_max_ref])                        
                    end
                end
                for jj = 1:length(p_corename{ii, curr_rad})
                    if (logical(p_corename{ii, curr_rad}(jj)) && ishandle(p_corename{ii, curr_rad}(jj)))
                        tmp1 = get(p_corename{ii, curr_rad}(jj), 'position');
                        set(p_corename{ii, curr_rad}(jj), 'position', [tmp1(1:2) (elev_max_ref - 50)])
                    end
                end
                if (any(p_coreflat{ii}) && any(ishandle(p_coreflat{ii})))
                    set(p_coreflat{ii}(logical(p_coreflat{ii}) & ishandle(p_coreflat{ii})), 'ydata', [depth_min_ref depth_max_ref])
                end
                for jj = 1:length(p_corenameflat{ii})
                    if (logical(p_corenameflat{ii}(jj)) && ishandle(p_corenameflat{ii}(jj)))
                        tmp1 = get(p_corenameflat{ii}(jj), 'position');
                        set(p_corenameflat{ii}(jj), 'position', [tmp1(1) (depth_min_ref + 50) 0])
                    end
                end
            end
        end
        
        % plot data
        data_done(curr_rad) = true;
        set(data_check(curr_rad), 'value', 1)
        plot_db
        set(status_box(2), 'string', 'Transect radar data loaded.')
    end

%% Choose the current layer

    function choose_layer1(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                        = deal(1, 1, 1, 2);
        choose_layer
    end

    function choose_layer2(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                        = deal(1, 1, 2, 1);
        choose_layer
    end

    function choose_layer3(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                        = deal(2, 2, 1, 2);
        choose_layer
    end

    function choose_layer4(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                        = deal(2, 3, 2, 1);
        choose_layer
    end

    function choose_layer(source, eventdata)
        curr_layer(curr_rad)= get(layer_list(curr_gui, curr_rad), 'value');
        if pk_done(curr_rad)
            set(layer_list(:, curr_rad), 'value', curr_layer(curr_rad))
            if (any(p_pk{1, curr_rad}) && any(ishandle(p_pk{1, curr_rad})))
                set(p_pk{1, curr_rad}(logical(p_pk{1, curr_rad}) & ishandle(p_pk{1, curr_rad})), 'markersize', 12)
            end
            if (any(p_pk{2, curr_rad}) && any(ishandle(p_pk{2, curr_rad})))
                set(p_pk{2, curr_rad}(logical(p_pk{2, curr_rad}) & ishandle(p_pk{2, curr_rad})), 'markersize', 12)
            end
            if (any(p_int2{1, curr_rad_alt}) && any(ishandle(p_int2{1, curr_rad_alt})))
                set(p_int2{1, curr_rad_alt}(logical(p_int2{1, curr_rad_alt}) & ishandle(p_int2{1, curr_rad_alt})), 'markersize', 8)
            end
            if (any(p_int2{2, curr_rad_alt}) && any(ishandle(p_int2{2, curr_rad_alt})))
                set(p_int2{2, curr_rad_alt}(logical(p_int2{2, curr_rad_alt}) & ishandle(p_int2{2, curr_rad_alt})), 'markersize', 8)
            end
            if (logical(p_pk{1, curr_rad}(curr_layer(curr_rad))) && ishandle(p_pk{1, curr_rad}(curr_layer(curr_rad))))
                set(p_pk{1, curr_rad}(curr_layer(curr_rad)), 'markersize', 24)
            end
            if (logical(p_pk{2, curr_rad}(curr_layer(curr_rad))) && ishandle(p_pk{2, curr_rad}(curr_layer(curr_rad))))
                set(p_pk{2, curr_rad}(curr_layer(curr_rad)), 'markersize', 24)
            end
            if all(pk_done)
                if (logical(p_int2{1, curr_rad_alt}(curr_layer(curr_rad))) && ishandle(p_int2{1, curr_rad_alt}(curr_layer(curr_rad))))
                    set(p_int2{1, curr_rad_alt}(curr_layer(curr_rad)), 'markersize', 16)
                end
                if (logical(p_int2{2, curr_rad_alt}(curr_layer(curr_rad))) && ishandle(p_int2{2, curr_rad_alt}(curr_layer(curr_rad))))
                    set(p_int2{2, curr_rad_alt}(curr_layer(curr_rad)), 'markersize', 16)
                end
            end
            if (flat_done(curr_rad) && data_done(curr_rad))
                if (any(p_pkflat{curr_rad}) && any(ishandle(p_pkflat{curr_rad})))
                    set(p_pkflat{curr_rad}(logical(p_pkflat{curr_rad}) & ishandle(p_pkflat{curr_rad})), 'markersize', 12)
                end
                if (logical(p_pkflat{curr_rad}(curr_layer(curr_rad))) && ishandle(p_pkflat{curr_rad}(curr_layer(curr_rad))))
                    set(p_pkflat{curr_rad}(curr_layer(curr_rad)), 'markersize', 24)
                end
            end
        end
        if ischar(tmp1)
            if strcmp(tmp1, 'rechoose')
                return
            end
        end
        if (all(pk_done) && ~isempty(pk{2}.ind_layer))
            switch curr_rad
                case 1
                    if ~isempty(find(pk{2}.ind_layer(:, 4) == curr_layer(1), 1))
                        curr_layer(2) = pk{2}.ind_layer(find(pk{2}.ind_layer(:, 4) == curr_layer(1), 1), 1);
                        if (any(p_pk{1, 2}) && any(ishandle(p_pk{1, 2})))
                            set(p_pk{1, 2}(logical(p_pk{1, 2}) & ishandle(p_pk{1, 2})), 'markersize', 12)
                        end
                        if (any(p_pk{2, 2}) && any(ishandle(p_pk{2, 2})))
                            set(p_pk{2, 2}(logical(p_pk{2, 2}) & ishandle(p_pk{1, 2})), 'markersize', 12)
                        end
                        if (logical(p_pk{1, 2}(curr_layer(2))) && ishandle(p_pk{1, 2}(curr_layer(2))))
                            set(p_pk{1, 2}(curr_layer(2)), 'markersize', 24)
                        end
                        if (logical(p_pk{2, 2}(curr_layer(2))) && ishandle(p_pk{2, 2}(curr_layer(2))))
                            set(p_pk{2, 2}(curr_layer(2)), 'markersize', 24)
                        end
                        if (any(p_int2{1, 2}) && any(ishandle(p_int2{1, 2})))
                            set(p_int2{1, 2}(logical(p_int2{1, 2}) & ishandle(p_int2{1, 2})), 'markersize', 8)
                        end
                        if (any(p_int2{2, 2}) && any(ishandle(p_int2{2, 2})))
                            set(p_int2{2, 2}(logical(p_int2{2, 2}) & ishandle(p_int2{2, 2})), 'markersize', 8)
                        end
                        if (logical(p_int2{1, 2}(curr_layer(1))) && ishandle(p_int2{1, 2}(curr_layer(1))))
                            set(p_int2{1, 2}(curr_layer(1)), 'markersize', 16)
                        end
                        if (logical(p_int2{2, 2}(curr_layer(1))) && ishandle(p_int2{2, 2}(curr_layer(1))))
                            set(p_int2{2, 2}(curr_layer(1)), 'markersize', 16)
                        end
                        if (flat_done(2) && data_done(2))
                            if (any(p_pkflat{2}) && any(ishandle(p_pkflat{2})))
                                set(p_pkflat{2}(logical(p_pkflat{2}) & ishandle(p_pkflat{2})), 'markersize', 12)
                            end
                            if (logical(p_pkflat{2}(curr_layer(2))) && ishandle(p_pkflat{2}(curr_layer(2))))
                                set(p_pkflat{2}(curr_layer(2)), 'markersize', 24)
                            end
                        end
                        tmp1 = 'matched';
                    end
                    set(layer_list(:, 2), 'value', curr_layer(2))
                case 2
                    if ~isempty(find(pk{2}.ind_layer(:, 1) == curr_layer(2), 1))
                        curr_layer(1) = pk{2}.ind_layer(find(pk{2}.ind_layer(:, 1) == curr_layer(2), 1), 4);
                        if (any(p_pk{1, 1}) && any(ishandle(p_pk{1, 2})))
                            set(p_pk{1, 1}(logical(p_pk{1, 1}) & ishandle(p_pk{1, 1})), 'markersize', 12)
                        end
                        if (any(p_pk{2, 1}) && any(ishandle(p_pk{2, 1})))
                            set(p_pk{2, 1}(logical(p_pk{2, 1}) & ishandle(p_pk{1, 1})), 'markersize', 12)
                        end
                        if (logical(p_pk{1, 1}(curr_layer(1))) && ishandle(p_pk{1, 1}(curr_layer(1))))
                            set(p_pk{1, 1}(curr_layer(1)), 'markersize', 24)
                        end
                        if (logical(p_pk{2, 1}(curr_layer(1))) && ishandle(p_pk{2, 1}(curr_layer(1))))
                            set(p_pk{2, 1}(curr_layer(1)), 'markersize', 24)
                        end
                        if (any(p_int2{1, 1}) && any(ishandle(p_int2{1, 1})))
                            set(p_int2{1, 1}(logical(p_int2{1, 1}) & ishandle(p_int2{1, 1})), 'markersize', 8)
                        end
                        if (any(p_int2{2, 1}) && any(ishandle(p_int2{2, 1})))
                            set(p_int2{2, 1}(logical(p_int2{2, 1}) & ishandle(p_int2{2, 1})), 'markersize', 8)
                        end
                        if (logical(p_int2{1, 1}(curr_layer(2))) && ishandle(p_int2{1, 1}(curr_layer(2))))
                            set(p_int2{1, 1}(curr_layer(2)), 'markersize', 16)
                        end
                        if (logical(p_int2{2, 1}(curr_layer(2))) && ishandle(p_int2{2, 1}(curr_layer(2))))
                            set(p_int2{2, 1}(curr_layer(2)), 'markersize', 16)
                        end
                        if (flat_done(1) && data_done(1))
                            if (any(p_pkflat{1}) && any(ishandle(p_pkflat{1})))
                                set(p_pkflat{1}(logical(p_pkflat{1}) & ishandle(p_pkflat{1})), 'markersize', 12)
                            end
                            if (logical(p_pkflat{1}(curr_layer(1))) && ishandle(p_pkflat{1}(curr_layer(1))))
                                set(p_pkflat{1}(curr_layer(1)), 'markersize', 24)
                            end
                        end
                        tmp1 = 'matched';
                    end
                    set(layer_list(:, 1), 'value', curr_layer(1))
            end
        end
        if (get(nearest_check, 'value') && all(pk_done) && ~strcmp(tmp1, 'matched'))
            tmp2            = 'nearest';
            if any(~isnan(pk{curr_rad_alt}.elev_smooth(:, curr_ind_int(curr_int, curr_rad_alt))))
                [~, curr_layer(curr_rad_alt)] ...
                            = min(abs(pk{curr_rad}.elev_smooth(curr_layer(curr_rad), curr_ind_int(curr_int, curr_rad)) - pk{curr_rad_alt}.elev_smooth(:, curr_ind_int(curr_int, curr_rad_alt))));
                [curr_rad, curr_rad_alt] ...
                            = deal(curr_rad_alt, curr_rad);
                set(layer_list(curr_gui, curr_rad), 'value', curr_layer(curr_rad))
                tmp1        = 'rechoose';
                choose_layer
            end
        end
    end

%% Choose a layer manually

    function choose_pk1(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                        = deal(2, 2, 1, 2);
        choose_pk
    end

    function choose_pk2(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                        = deal(2, 3, 2, 1);
        choose_pk
    end

    function choose_pk(source, eventdata)
        if ~pk_done(curr_rad)
            set(status_box(curr_gui), 'string', 'No layers to focus on.')
            return
        end
        if ~get(pk_check(2, curr_rad), 'value')
            set(pk_check(2, curr_rad), 'value', 1)
            show_pk
        end
        set(status_box(2), 'string', 'Choose a layer to highlight...')
        [ind_x_pk, ind_y_pk]= ginput(1);
        switch disp_type{curr_ax}
            case 'amp.'
                [tmp1, tmp2]= unique(pk{curr_rad}.elev_smooth(:, interp1(pk{curr_rad}.dist_lin(ind_decim{curr_gui, curr_rad}), ind_decim{curr_gui, curr_rad}, ind_x_pk, 'nearest', 'extrap')));
            case 'flat'
                [tmp1, tmp2]= unique(pk{curr_rad}.depth_smooth(:, interp1(pk{curr_rad}.dist_lin(ind_decim{curr_gui, curr_rad}), ind_decim{curr_gui, curr_rad}, ind_x_pk, 'nearest', 'extrap')));
        end
        tmp1                = tmp1(~isnan(tmp1));
        if (length(tmp1) > 1)
            curr_layer(curr_rad) ...
                            = interp1(tmp1, tmp2(~isnan(tmp1)), ind_y_pk, 'nearest', 'extrap');
        elseif (length(tmp2) > 1)
            set(status_box(2), 'string', 'Could not determine which layer to select. Try again.')
            return
        else
            curr_layer(curr_rad) ...
                            = tmp2;
        end
        set(layer_list(:, curr_rad), 'value', curr_layer(curr_rad))
        choose_layer
        set(status_box(curr_gui), 'string', ['Layer #' num2str(curr_layer(curr_rad)) ' chosen.'])
    end

%% Switch to previous layer in list

    function pk_last1(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                        = deal(1, 1, 1, 2);
        pk_last
    end

    function pk_last2(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                        = deal(1, 1, 2, 1);
        pk_last
    end

    function pk_last3(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                        = deal(2, 2, 1, 2);
        pk_last
    end

    function pk_last4(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                        = deal(2, 3, 2, 1);
        pk_last
    end

    function pk_last(source, eventdata)
        if (pk_done(curr_rad) && (curr_layer(curr_rad) > 1))
            curr_layer(curr_rad) ...
                            = curr_layer(curr_rad) - 1;
            set(layer_list(:, curr_rad), 'value', curr_layer(curr_rad))
            if (any(p_pk{1, curr_rad}) && any(ishandle(p_pk{1, curr_rad})))
                set(p_pk{1, curr_rad}(logical(p_pk{1, curr_rad}) & ishandle(p_pk{1, curr_rad})), 'markersize', 12)
            end
            if (logical(p_pk{1, curr_rad}(curr_layer(curr_rad))) && ishandle(p_pk{1, curr_rad}(curr_layer(curr_rad))))
                set(p_pk{1, curr_rad}(curr_layer(curr_rad)), 'markersize', 24)
            end
            if (any(p_pk{2, curr_rad}) && any(ishandle(p_pk{2, curr_rad})))
                set(p_pk{2, curr_rad}(logical(p_pk{1, curr_rad}) & ishandle(p_pk{2, curr_rad})), 'markersize', 12)
            end
            if (logical(p_pk{2, curr_rad}(curr_layer(curr_rad))) && ishandle(p_pk{2, curr_rad}(curr_layer(curr_rad))))
                set(p_pk{2, curr_rad}(curr_layer(curr_rad)), 'markersize', 24)
            end
            if (flat_done(curr_rad) && data_done(curr_rad))
                if (any(p_pkflat{curr_rad}) && any(ishandle(p_pkflat{curr_rad})))
                    set(p_pkflat{curr_rad}(logical(p_pkflat{curr_rad}) & ishandle(p_pkflat{curr_rad})), 'markersize', 12)
                end
                if (logical(p_pkflat{curr_rad}(curr_layer(curr_rad))) && ishandle(p_pkflat{curr_rad}(curr_layer(curr_rad))))
                    set(p_pkflat{curr_rad}(curr_layer(curr_rad)), 'markersize', 24)
                end
            end
        end
    end

%% Switch to next layer in the list

    function pk_next1(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                        = deal(1, 1, 1, 2);
        pk_next
    end

    function pk_next2(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                        = deal(1, 1, 2, 1);
        pk_next
    end

    function pk_next3(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                        = deal(2, 2, 1, 2);
        pk_next
    end

    function pk_next4(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                        = deal(2, 3, 2, 1);
        pk_next
    end

    function pk_next(source, eventdata)
        if (pk_done(curr_rad) && (curr_layer(curr_rad) < pk{curr_rad}.num_layer))
            curr_layer(curr_rad) ...
                            = curr_layer(curr_rad) + 1;
            set(layer_list(:, curr_rad), 'value', curr_layer(curr_rad))
            if (any(p_pk{1, curr_rad}) && any(ishandle(p_pk{1, curr_rad})))
                set(p_pk{1, curr_rad}(logical(p_pk{1, curr_rad}) & ishandle(p_pk{1, curr_rad})), 'markersize', 12)
            end
            if (logical(p_pk{1, curr_rad}(curr_layer(curr_rad))) && ishandle(p_pk{1, curr_rad}(curr_layer(curr_rad))))
                set(p_pk{1, curr_rad}(curr_layer(curr_rad)), 'markersize', 24)
            end
            if (any(p_pk{2, curr_rad}) && any(ishandle(p_pk{2, curr_rad})))
                set(p_pk{2, curr_rad}(logical(p_pk{1, curr_rad}) & ishandle(p_pk{2, curr_rad})), 'markersize', 12)
            end
            if (logical(p_pk{2, curr_rad}(curr_layer(curr_rad))) && ishandle(p_pk{2, curr_rad}(curr_layer(curr_rad))))
                set(p_pk{2, curr_rad}(curr_layer(curr_rad)), 'markersize', 24)
            end
            if (flat_done(curr_rad) && data_done(curr_rad))
                if (any(p_pkflat{curr_rad}) && any(ishandle(p_pkflat{curr_rad})))
                    set(p_pkflat{curr_rad}(logical(p_pkflat{curr_rad}) & ishandle(p_pkflat{curr_rad})), 'markersize', 12)
                end
                if (logical(p_pkflat{curr_rad}(curr_layer(curr_rad))) && ishandle(p_pkflat{curr_rad}(curr_layer(curr_rad))))
                    set(p_pkflat{curr_rad}(curr_layer(curr_rad)), 'markersize', 24)
                end
            end
        end
    end

%% Match two intersecting layers

    function pk_match(source, eventdata)
        
        curr_gui            = 2;
        
        if ~all(pk_done)
            set(status_box(2), 'string', 'Picks must be loaded for both master and intersecting transects.')
            return
        end
        
        if (logical(p_pk{1, 2}(curr_layer(2))) && ishandle(p_pk{1, 2}(curr_layer(2))))
           set(p_pk{1, 2}(curr_layer(2)), 'color', colors{1}(curr_layer(1), :))
        end
        if (logical(p_pk{2, 2}(curr_layer(2))) && ishandle(p_pk{2, 2}(curr_layer(2))))
           set(p_pk{2, 2}(curr_layer(2)), 'color', colors{1}(curr_layer(1), :))
        end
        if (logical(p_int2{1, 1}(curr_layer(2))) && ishandle(p_int2{1, 1}(curr_layer(2))))
           set(p_int2{1, 1}(curr_layer(2)), 'markerfacecolor', colors{1}(curr_layer(1), :))
        end
        if (logical(p_int2{2, 1}(curr_layer(2))) && ishandle(p_int2{2, 1}(curr_layer(2))))
           set(p_int2{2, 1}(curr_layer(2)), 'markerfacecolor', colors{1}(curr_layer(1), :))
        end
        if (logical(p_int2{1, 2}(curr_layer(1))) && ishandle(p_int2{1, 2}(curr_layer(1))))
           set(p_int2{1, 2}(curr_layer(1)), 'markerfacecolor', colors{1}(curr_layer(1), :))
        end
        if (logical(p_int2{2, 2}(curr_layer(1))) && ishandle(p_int2{2, 2}(curr_layer(1))))
           set(p_int2{2, 2}(curr_layer(1)), 'markerfacecolor', colors{1}(curr_layer(1), :))
        end
        
        pause(0.1)
        
        set(status_box(2), 'string', 'Matching correct? (Y: yes; otherwise: cancel)...')
        
        waitforbuttonpress
        
        if ~strcmpi(get(fgui(2), 'currentcharacter'), 'Y')
            set(p_pk{1, 2}(curr_layer(2)), 'color', colors{2}(curr_layer(2), :))
            set(p_pk{2, 2}(curr_layer(2)), 'color', colors{2}(curr_layer(2), :))
            set(p_int2{1, 1}(curr_layer(2)), 'markerfacecolor', colors{2}(curr_layer(2), :))
            set(p_int2{2, 1}(curr_layer(2)), 'markerfacecolor', colors{2}(curr_layer(2), :))
            set(p_int2{1, 2}(curr_layer(1)), 'markerfacecolor', colors{1}(curr_layer(1), :))
            set(p_int2{2, 2}(curr_layer(1)), 'markerfacecolor', colors{1}(curr_layer(1), :))
            set(status_box(2), 'string', 'Layer matching cancelled by user.')
            return
        end
        
        set([p_int2{1, 1}(curr_layer(2)) p_int2{2, 1}(curr_layer(2)) p_int2{1, 2}(curr_layer(1)) p_int2{2, 2}(curr_layer(1))], 'marker', '^')
        
        set(status_box(2), 'string', ['Matching master transect layer #' num2str(curr_layer(1)) ' with intersecting transect layer # ' num2str(curr_layer(2)) '...'])
        pause(0.1)
        
        pk{2}.ind_layer     = [pk{2}.ind_layer; [curr_layer(2) curr_year(1) curr_trans(1) curr_layer(1) (pk{1}.ind_layer(curr_layer(1), end) + 1)]];
        
        set(status_box(2), 'string', ['Intersecting layer # ' num2str(curr_layer(2)) ' matched to master layer #' num2str(curr_layer(1)) '.'])
    end

%% Match two intersecting layers

    function pk_unmatch(source, eventdata)
        
        curr_gui            = 2;
        
        if ~all(pk_done)
            set(status_box(2), 'string', 'Picks must be loaded for both master and intersecting transects.')
            return
        end
        
        set(status_box(2), 'string', 'Unmatch current pair? (Y: yes; otherwise: cancel)...')
        
        waitforbuttonpress
        
        if strcmpi(get(fgui(2), 'currentcharacter'), 'Y')
            
            set(status_box(2), 'string', ['Unmatching master transect layer #' num2str(curr_layer(1)) ' from intersecting transect layer # ' num2str(curr_layer(2)) '...'])
            pause(0.1)
            
            if ~isempty(find(((pk{2}.ind_layer(:, 1) == curr_layer(2)) & (pk{2}.ind_layer(:, 4) == curr_layer(1))), 1))
                pk{2}.ind_layer = pk{2}.ind_layer(setdiff(1:size(pk{2}.ind_layer, 1), find(((pk{2}.ind_layer(:, 1) == curr_layer(2)) & (pk{2}.ind_layer(:, 4) == curr_layer(1))))), :);
            else
                set(status_box(2), 'string', 'Current layer pair not matched. Unmatching cancelled.')
                return
            end
            
            set(p_pk{1, 2}(curr_layer(2)), 'color', colors{2}(curr_layer(2), :))
            set(p_pk{2, 2}(curr_layer(2)), 'color', colors{2}(curr_layer(2), :))
            set(p_int2{1, 1}(curr_layer(2)), 'markerfacecolor', colors{2}(curr_layer(2), :))
            set(p_int2{2, 1}(curr_layer(2)), 'markerfacecolor', colors{2}(curr_layer(2), :))
            set(p_int2{1, 2}(curr_layer(1)), 'markerfacecolor', colors{1}(curr_layer(1), :))
            set(p_int2{2, 2}(curr_layer(1)), 'markerfacecolor', colors{1}(curr_layer(1), :))            
            if ~isempty(find((pk{2}.ind_layer(:, 1) == curr_layer(2)), 1))
                set(p_int2{1, 1}(curr_layer(2)), 'marker', 's')
                set(p_int2{2, 1}(curr_layer(2)), 'marker', 's')
            else
                set(p_int2{1, 1}(curr_layer(2)), 'marker', 'o')
                set(p_int2{2, 1}(curr_layer(2)), 'marker', 'o')
            end
            if isempty(find((pk{2}.ind_layer(:, 4) == curr_layer(1)), 1))
                set(p_int2{1, 2}(curr_layer(1)), 'marker', 'o')
                set(p_int2{2, 2}(curr_layer(1)), 'marker', 'o')
            end
            
            set(status_box(2), 'string', ['Intersecting layer # ' num2str(curr_layer(2)) ' unmatched from master layer #' num2str(curr_layer(1)) '.'])
            
        else
            set(statux_box(2), 'string', 'Unmatching cancelled by user.')
            return
        end
    end

%% Save merged picks

    function pk_save(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(1, 1, 2, 1);
        % want everything done before saving
        if ~pk_done(2)
            set(status_box(1), 'string', 'Intersecting picks not loaded yet.')
            return
        end
        if isempty(pk{2}.ind_layer)
            set(status_box(1), 'string', 'No layers matched. No need to save intersecting picks.')
            return
        end
        
        pk{2}.ind_layer     = sortrows(pk{2}.ind_layer);
        tmp1                = pk{2}.ind_layer(:, [1 4]);
        [~, tmp1]           = unique(tmp1, 'rows', 'first');
        pk{2}.ind_layer     = pk{2}.ind_layer(tmp1, :);
        pk{2}               = orderfields(pk{2});
        
        if ~isempty(path_pk{2})
            [file_pk{2}, path_pk{2}] = uiputfile('*.mat', 'Save picks:', [path_pk{2} file_pk{2}]);
        elseif ~isempty(path_data{2})
            [file_pk{2}, path_pk{2}] = uiputfile('*.mat', 'Save picks:', [path_data{2} file_pk{2}]);
        elseif ~isempty(path_pk{1})
            [file_pk{2}, path_pk{2}] = uiputfile('*.mat', 'Save picks:', [path_pk{1} file_pk{2}]);
        elseif ~isempty(path_data{1})
            [file_pk{2}, path_pk{2}] = uiputfile('*.mat', 'Save picks:', [path_data{1} file_pk{2}]);
        else
            [file_pk{2}, path_pk{2}] = uiputfile('*.mat', 'Save picks:', file_pk{2});
        end
        
        if ~file_pk{2}
            [file_pk{2}, path_pk{2}] ...
                            = deal('');
        else
            set(status_box(1), 'string', 'Saving intersecting picks...')
            pause(0.1)
            tmp1            = pk;
            pk              = pk{2};
            save([path_pk{2} file_pk{2}], '-v7.3', 'pk')
            pk              = tmp1;
            set(status_box(1), 'string', ['Intersecting picks saved as ' file_pk{2}(1:(end - 4)) ' in ' path_pk{2} '.'])
        end
    end

%% Update minimum dB/x/y/z/dist

    function slide_db_min1(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 2, 1, 2);
        slide_db_min
    end

    function slide_db_min2(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 3, 2, 1);
        slide_db_min
    end

    function slide_db_min(source, eventdata)
        if (get(cb_min_slide(curr_rad), 'value') < db_max(curr_rad))
            if get(cbfix_check1(curr_rad), 'value')
                tmp1        = db_max(curr_rad) - db_min(curr_rad);
            end
            db_min(curr_rad) = get(cb_min_slide(curr_rad), 'value');
            if get(cbfix_check1(curr_rad), 'value')
                db_max(curr_rad) = db_min(curr_rad) + tmp1;
                if (db_max(curr_rad) > db_max_ref(curr_rad))
                    db_max(curr_rad) = db_max_ref(curr_rad);
                    db_min(curr_rad) = db_max(curr_rad) - tmp1;
                    if (db_min(curr_rad) < db_min_ref(curr_rad))
                        db_min(curr_rad) = db_min_ref(curr_rad);
                    end
                    if (db_min(curr_rad) < get(cb_min_slide(curr_rad), 'min'))
                        set(cb_min_slide(curr_rad), 'value', get(cb_min_slide(curr_rad), 'min'))
                    else
                        set(cb_min_slide(curr_rad), 'value', db_min(curr_rad))
                    end
                end
                set(cb_max_edit(curr_rad), 'string', sprintf('%3.0f', db_max(curr_rad)))
                if (db_max(curr_rad) > get(cb_max_slide(curr_rad), 'max'))
                    set(cb_max_slide(curr_rad), 'value', get(cb_max_slide(curr_rad), 'max'))
                else
                    set(cb_max_slide(curr_rad), 'value', db_max(curr_rad))
                end
            end
            set(cb_min_edit(curr_rad), 'string', sprintf('%3.0f', db_min(curr_rad)))
            update_db_range
        else
            if (db_min(curr_rad) < get(cb_min_slide(curr_rad), 'min'))
                set(cb_min_slide(curr_rad), 'value', get(cb_min_slide(curr_rad), 'min'))
            else
                set(cb_min_slide(curr_rad), 'value', db_min(curr_rad))
            end
        end
        set(cb_min_slide(curr_rad), 'enable', 'off')
        drawnow
        set(cb_min_slide(curr_rad), 'enable', 'on')
    end

    function slide_dist_min1(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 2, 1, 2);
        slide_dist_min
    end

    function slide_dist_min2(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 3, 2, 1);
        slide_dist_min
    end

    function slide_dist_min(source, eventdata)
        if (get(dist_min_slide(curr_rad), 'value') < dist_max(curr_rad))
            if get(distfix_check(curr_rad), 'value')
                tmp1        = dist_max(curr_rad) - dist_min(curr_rad);
            end
            dist_min(curr_rad) = get(dist_min_slide(curr_rad), 'value');
            if get(distfix_check(curr_rad), 'value')
                dist_max(curr_rad) = dist_min(curr_rad) + tmp1;
                if (dist_max(curr_rad) > dist_max_ref(curr_rad))
                    dist_max(curr_rad) = dist_max_ref(curr_rad);
                    dist_min(curr_rad) = dist_max(curr_rad) - tmp1;
                    if (dist_min(curr_rad) < dist_min_ref(curr_rad))
                        dist_min(curr_rad) = dist_min_ref(curr_rad);
                    end
                    if (dist_min(curr_rad) < get(dist_min_slide(curr_rad), 'min'))
                        set(dist_min_slide(curr_rad), 'value', get(dist_min_slide(curr_rad), 'min'))
                    else
                        set(dist_min_slide(curr_rad), 'value', dist_min(curr_rad))
                    end
                end
                set(dist_max_edit(curr_rad), 'string', sprintf('%3.0f', dist_max(curr_rad)))
                if (dist_max(curr_rad) > get(dist_max_slide(curr_rad), 'max'))
                    set(dist_max_slide(curr_rad), 'value', get(dist_max_slide(curr_rad), 'max'))
                else
                    set(dist_max_slide(curr_rad), 'value', dist_max(curr_rad))
                end
            end
            set(dist_min_edit(curr_rad), 'string', sprintf('%3.0f', dist_min(curr_rad)))
            update_dist_range
        else
            if (dist_min(curr_rad) < get(dist_min_slide(curr_rad), 'min'))
                set(dist_min_slide(curr_rad), 'value', get(dist_min_slide(curr_rad), 'min'))
            else
                set(dist_min_slide(curr_rad), 'value', dist_min(curr_rad))
            end
        end
        set(dist_min_slide(curr_rad), 'enable', 'off')
        drawnow
        set(dist_min_slide(curr_rad), 'enable', 'on')
    end

    function slide_x_min(source, eventdata)
        [curr_gui, curr_ax] = deal(1);
        if (get(x_min_slide, 'value') < x_max)
            if get(xfix_check, 'value')
                tmp1        = x_max - x_min;
            end
            x_min           = get(x_min_slide, 'value');
            if get(xfix_check, 'value')
                x_max       = x_min + tmp1;
                if (x_max > x_max_ref)
                    x_max   = x_max_ref;
                    x_min   = x_max - tmp1;
                    if (x_min < x_min_ref)
                        x_min = x_min_ref;
                    end
                    if (x_min < get(x_min_slide, 'min'))
                        set(x_min_slide, 'value', get(x_min_slide, 'min'))
                    else
                        set(x_min_slide, 'value', x_min)
                    end
                end
                if (x_max > get(x_max_slide, 'max'))
                    set(x_max_slide, 'value', get(x_max_slide, 'max'))
                else
                    set(x_max_slide, 'value', x_max)
                end
                set(x_max_edit, 'string', sprintf('%4.1f', x_max))
            end
            set(x_min_edit, 'string', sprintf('%4.1f', x_min))
            update_x_range
        else
            if (x_min < get(x_min_slide, 'min'))
                set(x_min_slide, 'value', get(x_min_slide, 'min'))
            else
                set(x_min_slide, 'value', x_min)
            end
        end
        set(x_min_slide, 'enable', 'off')
        drawnow
        set(x_min_slide, 'enable', 'on')
    end

    function slide_y_min(source, eventdata)
        [curr_gui, curr_ax] = deal(1);
        if (get(y_min_slide, 'value') < y_max)
            if get(yfix_check, 'value')
                tmp1        = y_max - y_min;
            end
            y_min           = get(y_min_slide, 'value');
            if get(yfix_check, 'value')
                y_max       = y_min + tmp1;
                if (y_max > y_max_ref)
                    y_max   = y_max_ref;
                    y_min   = y_max - tmp1;
                    if (y_min < y_min_ref)
                        y_min = y_min_ref;
                    end
                    if (y_min < get(y_min_slide, 'min'))
                        set(y_min_slide, 'value', get(y_min_slide, 'min'))
                    else
                        set(y_min_slide, 'value', y_min)
                    end
                end
                if (y_max > get(y_max_slide, 'max'))
                    set(y_max_slide, 'value', get(y_max_slide, 'max'))
                else
                    set(y_max_slide, 'value', y_max)
                end
                set(y_max_edit, 'string', sprintf('%4.1f', y_max))
            end
            set(y_min_edit, 'string', sprintf('%4.1f', y_min))
            update_y_range
        else
            if (y_min < get(y_min_slide, 'min'))
                set(y_min_slide, 'value', get(y_min_slide, 'min'))
            else
                set(y_min_slide, 'value', y_min)
            end
        end
        set(y_min_slide, 'enable', 'off')
        drawnow
        set(y_min_slide, 'enable', 'on')
    end

    function slide_z_min1(source, eventdata)
        [curr_gui, curr_ax] = deal(1);
        slide_z_min
    end

    function slide_z_min2(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 2, 1, 2);
        slide_z_min
    end

    function slide_z_min3(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 3, 2, 1);
        slide_z_min
    end

    function slide_z_min(source, eventdata)
        switch disp_type{curr_ax}
            case 'amp.'
                if (get(z_min_slide(curr_ax), 'value') < elev_max(curr_gui))
                    if (flat_done(curr_rad) && data_done(curr_rad))
                        tmp2        = [elev_min(curr_gui) elev_max(curr_gui)];
                    end
                    if get(zfix_check(curr_ax), 'value')
                        tmp1        = elev_max(curr_gui) - elev_min(curr_gui);
                    end
                    elev_min(curr_gui) = get(z_min_slide(curr_ax), 'value');
                    if get(zfix_check(curr_ax), 'value')
                        elev_max(curr_gui) = elev_min(curr_gui) + tmp1;
                        if (elev_max(curr_gui) > elev_max_ref)
                            elev_max(curr_gui) = elev_max_ref;
                            elev_min(curr_gui) = elev_max(curr_gui) - tmp1;
                            if (elev_min(curr_gui) < elev_min_ref)
                                elev_min(curr_gui) = elev_min_ref;
                            end
                            if (elev_min(curr_gui) < get(z_min_slide(curr_ax), 'min'))
                                set(z_min_slide(curr_ax), 'value', get(z_min_slide(curr_ax), 'min'))
                            else
                                set(z_min_slide(curr_ax), 'value', elev_min(curr_gui))
                            end
                        end
                        set(z_max_edit(curr_ax), 'string', sprintf('%4.0f', elev_max(curr_gui)))
                        if (elev_max(curr_gui) > get(z_max_slide(curr_ax), 'max'))
                            set(z_max_slide(curr_ax), 'value', get(z_max_slide(curr_ax), 'max'))
                        else
                            set(z_max_slide(curr_ax), 'value', elev_max(curr_gui))
                        end
                    end
                    if (flat_done(curr_rad) && data_done(curr_rad))
                        if get(zfix_check(curr_ax), 'value')
                            depth_min = depth_min - (elev_max(curr_gui) - tmp2(2));
                            depth_max = depth_max - (elev_min(curr_gui) - tmp2(1));
                            if (depth_min < depth_min_ref)
                                depth_min = depth_min_ref;
                            end
                            if (depth_max > depth_max_ref)
                                depth_max = depth_max_ref;
                            end
                        elseif ((depth_max - (elev_min(curr_gui) - tmp2(1))) > depth_min)
                            depth_max   = (depth_max - (elev_min(curr_gui) - tmp2(1)));
                        end
                    end
                    set(z_min_edit(curr_ax), 'string', sprintf('%4.0f', elev_min(curr_gui)))
                    update_z_range
                else
                    if (elev_min(curr_gui) < get(z_min_slide(curr_ax), 'min'))
                        set(z_min_slide(curr_ax), 'value', get(z_min_slide(curr_ax), 'min'))
                    else
                        set(z_min_slide(curr_ax), 'value', elev_min(curr_gui))
                    end
                end
            case 'flat'
                if ((depth_max_ref - (get(z_min_slide(curr_ax), 'value') - depth_min_ref)) > depth_min)
                    tmp2            = [depth_min depth_max];
                    if get(zfix_check(curr_ax), 'value')
                        tmp1        = depth_max - depth_min;
                    end
                    depth_max = depth_max_ref - (get(z_min_slide(curr_ax), 'value') - depth_min_ref);
                    if get(zfix_check(curr_ax), 'value')
                        depth_min    = depth_max - tmp1;
                        if (depth_min < depth_min_ref)
                            depth_min = depth_min_ref;
                            depth_max = depth_min + tmp1;
                            if (depth_max > depth_max_ref)
                                depth_max = depth_max_ref;
                            end
                            if ((depth_max_ref - (depth_max - depth_min_ref)) < get(z_min_slide(curr_ax), 'min'))
                                set(z_min_slide(curr_ax), 'value', get(z_min_slide(curr_ax), 'min'))
                            else
                                set(z_min_slide(curr_ax), 'value', (depth_max_ref - (depth_max - depth_min_ref)))
                            end
                            set(z_min_slide(curr_ax), 'value', (depth_max_ref - (depth_max - depth_min_ref)))
                        end
                        set(z_min_edit(curr_ax), 'string', sprintf('%4.0f', depth_max))
                        if ((depth_max_ref - (depth_min - depth_min_ref)) > get(z_max_slide(curr_ax), 'max'))
                            set(z_max_slide(curr_ax), 'value', get(z_max_slide(curr_ax), 'max'))
                        else
                            set(z_max_slide(curr_ax), 'value', (depth_max_ref - (depth_min - depth_min_ref)))
                        end
                    end
                    if get(zfix_check(curr_ax), 'value')
                        elev_min(curr_gui) = elev_min(curr_gui) - (depth_max - tmp2(2));
                        elev_max(curr_gui) = elev_max(curr_gui) - (depth_min - tmp2(1));
                        if (elev_min(curr_gui) < elev_min_ref)
                            elev_min(curr_gui) = elev_min_ref;
                        end
                        if (elev_max(curr_gui) > elev_max_ref)
                            elev_max(curr_gui) = elev_max_ref;
                        end
                    elseif (elev_min(curr_gui) - (depth_max - tmp2(2)) < elev_max)
                        elev_min        = elev_min(curr_gui) - (depth_max - tmp2(2));
                    end
                    set(z_min_edit(curr_ax), 'string', sprintf('%4.0f', depth_max))
                    update_z_range
                else
                    if ((depth_max_ref - (depth_max - depth_min_ref)) < get(z_min_slide(curr_ax), 'min'))
                        set(z_min_slide(curr_ax), 'value', get(z_min_slide(curr_ax), 'min'))
                    else
                        set(z_min_slide(curr_ax), 'value', (depth_max_ref - (depth_max - depth_min_ref)))
                    end
                end
        end
        set(z_min_slide(curr_ax), 'enable', 'off')
        drawnow
        set(z_min_slide(curr_ax), 'enable', 'on')
    end

%% Update maximum dB/x/y/z/dist

    function slide_db_max1(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 2, 1, 2);
        slide_db_max
    end

    function slide_db_max2(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 3, 2, 1);
        slide_db_max
    end

    function slide_db_max(source, eventdata)
        if (get(cb_max_slide(curr_rad), 'value') > db_min(curr_rad))
            if get(cbfix_check1(curr_rad), 'value')
                tmp1        = db_max(curr_rad) - db_min(curr_rad);
            end
            db_max(curr_rad) = get(cb_max_slide(curr_rad), 'value');
            if get(cbfix_check1(curr_rad), 'value')
                db_min(curr_rad) = db_max(curr_rad) - tmp1;
                if (db_min(curr_rad) < db_min_ref(curr_rad))
                    db_min(curr_rad) = db_min_ref(curr_rad);
                    db_max(curr_rad) = db_min(curr_rad) + tmp1;
                    if (db_max(curr_rad) > db_max_ref(curr_rad))
                        db_max(curr_rad) = db_max_ref(curr_rad);
                    end
                    if (db_max(curr_rad) > get(cb_max_slide(curr_rad), 'max'))
                        set(cb_max_slide(curr_rad), 'value', get(cb_max_slide(curr_rad), 'max'))
                    else
                        set(cb_max_slide(curr_rad), 'value', db_max(curr_rad))
                    end
                end
                set(cb_min_edit(curr_rad), 'string', sprintf('%3.0f', db_min(curr_rad)))
                if (db_min(curr_rad) < get(cb_min_slide(curr_rad), 'min'))
                    set(cb_min_slide(curr_rad), 'value', get(cb_min_slide(curr_rad), 'min'))
                else
                    set(cb_min_slide(curr_rad), 'value', db_min(curr_rad))
                end
            end
            set(cb_max_edit(curr_rad), 'string', sprintf('%3.0f', db_max(curr_rad)))
            update_db_range
        else
            if (db_max(curr_rad) > get(cb_max_slide(curr_rad), 'max'))
                set(cb_max_slide(curr_rad), 'value', get(cb_max_slide(curr_rad), 'max'))
            else
                set(cb_max_slide(curr_rad), 'value', db_max(curr_rad))
            end
        end
        set(cb_max_slide(curr_rad), 'enable', 'off')
        drawnow
        set(cb_max_slide(curr_rad), 'enable', 'on')
    end

    function slide_dist_max1(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 2, 1, 2);
        slide_dist_max
    end

    function slide_dist_max2(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 3, 2, 1);
        slide_dist_max
    end

    function slide_dist_max(source, eventdata)
        if (get(dist_max_slide(curr_rad), 'value') > dist_min(curr_rad))
            if get(distfix_check(curr_rad), 'value')
                tmp1        = dist_max(curr_rad) - dist_min(curr_rad);
            end
            dist_max(curr_rad) = get(dist_max_slide(curr_rad), 'value');
            if get(distfix_check(curr_rad), 'value')
                dist_min(curr_rad) = dist_max(curr_rad) - tmp1;
                if (dist_min(curr_rad) < dist_min_ref(curr_rad))
                    dist_min(curr_rad) = dist_min_ref(curr_rad);
                    dist_max(curr_rad) = dist_min(curr_rad) + tmp1;
                    if (dist_max(curr_rad) > dist_max_ref(curr_rad))
                        dist_max(curr_rad) = dist_max_ref(curr_rad);
                    end
                    if (dist_max(curr_rad) > get(dist_max_slide(curr_rad), 'max'))
                        set(dist_max_slide(curr_rad), 'value', get(dist_max_slide(curr_rad), 'max'))
                    else
                        set(dist_max_slide(curr_rad), 'value', dist_max(curr_rad))
                    end
                end
                set(dist_min_edit(curr_rad), 'string', sprintf('%3.0f', dist_min(curr_rad)))
                if (dist_min(curr_rad) < get(dist_min_slide(curr_rad), 'min'))
                    set(dist_min_slide(curr_rad), 'value', get(dist_min_slide(curr_rad), 'min'))
                else
                    set(dist_min_slide(curr_rad), 'value', dist_min(curr_rad))
                end
            end
            set(dist_max_edit(curr_rad), 'string', sprintf('%3.0f', dist_max(curr_rad)))
            update_dist_range
        else
            if (dist_max(curr_rad) > get(dist_max_slide(curr_rad), 'max'))
                set(dist_max_slide(curr_rad), 'value', get(dist_max_slide(curr_rad), 'max'))
            else
                set(dist_max_slide(curr_rad), 'value', dist_max(curr_rad))
            end
        end
        set(dist_max_slide(curr_rad), 'enable', 'off')
        drawnow
        set(dist_max_slide(curr_rad), 'enable', 'on')
    end

    function slide_x_max(source, eventdata)
        [curr_gui, curr_ax] = deal(1);
        if (get(x_max_slide, 'value') > x_min)
            if get(xfix_check, 'value')
                tmp1        = x_max - x_min;
            end
            x_max           = get(x_max_slide, 'value');
            if get(xfix_check, 'value')
                x_min       = x_max - tmp1;
                if (x_min < x_min_ref)
                    x_min   = x_min_ref;
                    x_max   = x_min + tmp1;
                    if (x_max > x_max_ref)
                        x_max = x_max_ref;
                    end
                    if (x_max > get(x_max_slide, 'max'))
                        set(x_max_slide, 'value', get(x_max_slide, 'max'))
                    else
                        set(x_max_slide, 'value', x_max)
                    end
                end
                if (x_min < get(x_min_slide, 'min'))
                    set(x_min_slide, 'value', get(x_min_slide, 'min'))
                else
                    set(x_min_slide, 'value', x_min)
                end
                set(x_min_edit, 'string', sprintf('%4.1f', x_min))
            end
            set(x_max_edit, 'string', sprintf('%4.1f', x_max))
            update_x_range
        else
            if (x_max > get(x_max_slide, 'max'))
                set(x_max_slide, 'value', get(x_max_slide, 'max'))
            else
                set(x_max_slide, 'value', x_max)
            end
        end
        set(x_max_slide, 'enable', 'off')
        drawnow
        set(x_max_slide, 'enable', 'on')
    end

    function slide_y_max(source, eventdata)
        [curr_gui, curr_ax] = deal(1);
        if (get(y_max_slide, 'value') > y_min)
            if get(yfix_check, 'value')
                tmp1        = y_max - y_min;
            end
            y_max           = get(y_max_slide, 'value');
            if get(yfix_check, 'value')
                y_min       = y_max - tmp1;
                if (y_min < y_min_ref)
                    y_min   = y_min_ref;
                    y_max   = y_min + tmp1;
                    if (y_max > y_max_ref)
                        y_max = y_max_ref;
                    end
                    if (y_max > get(y_max_slide, 'max'))
                        set(y_max_slide, 'value', get(y_max_slide, 'max'))
                    else
                        set(y_max_slide, 'value', y_max)
                    end
                end
                if (y_min < get(y_min_slide, 'min'))
                    set(y_min_slide, 'value', get(y_min_slide, 'min'))
                else
                    set(y_min_slide, 'value', y_min)
                end
                set(y_min_edit, 'string', sprintf('%4.1f', y_min))
            end
            set(y_max_edit, 'string', sprintf('%4.1f', y_max))
            update_y_range
        else
            if (y_max > get(y_max_slide, 'max'))
                set(y_max_slide, 'value', get(y_max_slide, 'max'))
            else
                set(y_max_slide, 'value', y_max)
            end
        end
        set(y_max_slide, 'enable', 'off')
        drawnow
        set(y_max_slide, 'enable', 'on')
    end

    function slide_z_max1(source, eventdata)
        [curr_gui, curr_ax] = deal(1);
        slide_z_max
    end

    function slide_z_max2(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 2, 1, 2);
        slide_z_max
    end

    function slide_z_max3(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 3, 2, 1);
        slide_z_max
    end

    function slide_z_max(source, eventdata)
        switch disp_type{curr_ax}
            case 'amp.'
                if (get(z_max_slide(curr_ax), 'value') > elev_min(curr_gui))
                    if (flat_done(curr_rad) && data_done(curr_rad))
                        tmp2        = [elev_min(curr_gui) elev_max(curr_gui)];
                    end
                    if get(zfix_check(curr_ax), 'value')
                        tmp1        = elev_max(curr_gui) - elev_min(curr_gui);
                    end
                    elev_max(curr_gui) = get(z_max_slide(curr_ax), 'value');
                    if get(zfix_check(curr_ax), 'value')
                        elev_min(curr_gui) = elev_max(curr_gui) - tmp1;
                        if (elev_min(curr_gui) < elev_min_ref)
                            elev_min(curr_gui) = elev_min_ref;
                            elev_max(curr_gui) = elev_min(curr_gui) + tmp1;
                            if (elev_max(curr_gui) > get(z_max_slide(curr_ax), 'max'))
                                set(z_max_slide(curr_ax), 'value', get(z_max_slide(curr_ax), 'max'))
                            else
                                set(z_max_slide(curr_ax), 'value', elev_max(curr_gui))
                            end
                        end
                        set(z_min_edit(curr_ax), 'string', sprintf('%4.0f', elev_min(curr_gui)))
                        if (elev_min(curr_gui) < get(z_min_slide(curr_ax), 'min'))
                            set(z_min_slide(curr_ax), 'value', get(z_min_slide(curr_ax), 'min'))
                        else
                            set(z_min_slide(curr_ax), 'value', elev_min(curr_gui))
                        end
                    end
                    if (flat_done(curr_rad) && data_done(curr_rad))
                        if get(zfix_check(curr_gui), 'value')
                            depth_min = depth_min - (elev_max(curr_gui) - tmp2(2));
                            depth_max = depth_max - (elev_min(curr_gui) - tmp2(1));
                            if (depth_min < depth_min_ref)
                                depth_min = depth_min_ref;
                            end
                            if (depth_max > depth_max_ref)
                                depth_max = depth_max_ref;
                            end
                        elseif (depth_min - (elev_max(curr_gui) - tmp2(2)) < depth_max)
                            depth_min   = depth_min - (elev_max(curr_gui) - tmp2(2));
                        end
                    end
                    set(z_max_edit(curr_ax), 'string', sprintf('%4.0f', elev_max(curr_gui)))
                    update_z_range
                else
                    if (elev_max(curr_gui) > get(z_max_slide(curr_ax), 'max'))
                        set(z_max_slide(curr_ax), 'value', get(z_max_slide(curr_ax), 'max'))
                    else
                        set(z_max_slide(curr_ax), 'value', elev_max(curr_gui))
                    end
                end
            case 'flat'
                if ((depth_max_ref - (get(z_max_slide(curr_ax), 'value') - depth_min_ref)) < depth_max)
                    tmp2            = [depth_min depth_max];
                    if get(zfix_check(curr_ax), 'value')
                        tmp1        = depth_max - depth_min;
                    end
                    depth_min = depth_max_ref - (get(z_max_slide(curr_ax), 'value') - depth_min_ref);
                    if get(zfix_check(curr_ax), 'value')
                        depth_max = depth_min + tmp1;
                        if (depth_max > depth_max_ref)
                            depth_max = depth_max_ref;
                            depth_min = depth_max - tmp1;
                            if ((depth_max_ref - (depth_min - depth_min_ref)) > get(z_max_slide(curr_ax), 'max'))
                                set(z_max_slide(curr_ax), 'min', get(z_max_slide(curr_ax), 'max'))
                            else
                                set(z_max_slide(curr_ax), 'value', (depth_max_ref - (depth_min - depth_min_ref)))
                            end
                        end
                        set(z_min_edit(curr_ax), 'string', sprintf('%4.0f', depth_max))
                        if ((depth_max_ref - (depth_max - depth_min_ref)) < get(z_min_slide(curr_ax), 'min'))
                            set(z_min_slide(curr_ax), 'value', get(z_min_slide(curr_ax), 'min'))
                        else
                            set(z_min_slide(curr_ax), 'value', (depth_max_ref - (depth_max - depth_min_ref)))
                        end
                    end
                    if get(zfix_check(curr_ax), 'value')
                        elev_min(curr_gui) = elev_min(curr_gui) - (depth_max - tmp2(2));
                        elev_max(curr_gui) = elev_max(curr_gui) - (depth_min - tmp2(1));
                        if (elev_min(curr_gui) < elev_min_ref)
                            elev_min(curr_gui) = elev_min_ref;
                        end
                        if (elev_max(curr_gui) > elev_max_ref)
                            elev_max(curr_gui) = elev_max_ref;
                        end
                    elseif (elev_min(curr_gui) - (depth_max - tmp2(2)) < elev_max)
                        elev_min = elev_min(curr_gui) - (depth_max - tmp2(2));
                    end
                    set(z_max_edit(curr_ax), 'string', sprintf('%4.0f', depth_min))
                    update_z_range
                else
                    if ((depth_max_ref - (depth_min - depth_min_ref)) > get(z_max_slide(curr_ax), 'max'))
                        set(z_max_slide(curr_ax), 'value', get(z_max_slide(curr_ax), 'max'))
                    else
                        set(z_max_slide(curr_ax), 'value', (depth_max_ref - (depth_min - depth_min_ref)))
                    end
                end
        end
        set(z_max_slide(curr_ax), 'enable', 'off')
        drawnow
        set(z_max_slide(curr_ax), 'enable', 'on')
    end

%% Reset minimum dB/x/y/z/dist

    function reset_db_min1(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 2, 1, 2);
        reset_db_min
    end

    function reset_db_min2(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 3, 2, 1);
        reset_db_min
    end

    function reset_db_min(source, eventdata)
        if (db_min_ref(curr_rad) < get(cb_min_slide(curr_rad), 'min'))
            set(cb_min_slide(curr_rad), 'value', get(cb_min_slide(curr_rad), 'min'))
        else
            set(cb_min_slide(curr_rad), 'value', db_min_ref(curr_rad))
        end
        set(cb_min_edit(curr_rad), 'string', num2str(db_min_ref(curr_rad)))
        db_min(curr_rad)    = db_min_ref(curr_rad);
        update_db_range
    end

    function reset_dist_min1(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 2, 1, 2);
        reset_dist_min
    end

    function reset_dist_min2(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 3, 2, 1);
        reset_dist_min
    end

    function reset_dist_min(source, eventdata)
        if (dist_min_ref(curr_rad) < get(dist_min_slide(curr_rad), 'min'))
            set(dist_min_slide(curr_rad), 'value', get(dist_min_slide(curr_rad), 'min'))
        else
            set(dist_min_slide(curr_rad), 'value', dist_min_ref(curr_rad))
        end
        set(dist_min_edit(curr_rad), 'string', sprintf('%3.1f', dist_min_ref(curr_rad)))
        dist_min(curr_rad)  = dist_min_ref(curr_rad);
        update_dist_range
    end

    function reset_x_min(source, eventdata)
        [curr_gui, curr_ax] = deal(1);
        if (x_min_ref < get(x_min_slide, 'min'))
            set(x_min_slide, 'value', get(x_min_slide, 'min'))
        else
            set(x_min_slide, 'value', x_min_ref)
        end
        set(x_min_edit, 'string', sprintf('%4.1f', x_min_ref))
        x_min              = x_min_ref;
        update_x_range
    end

    function reset_y_min(source, eventdata)
        [curr_gui, curr_ax] = deal(1);
        if (y_min_ref < get(y_min_slide, 'min'))
            set(y_min_slide, 'value', get(y_min_slide, 'min'))
        else
            set(y_min_slide, 'value', y_min_ref)
        end
        set(y_min_edit, 'string', sprintf('%4.1f', y_min_ref))
        y_min              = y_min_ref;
        update_y_range
    end

    function reset_z_min1(source, eventdata)
        [curr_gui, curr_ax] = deal(1);
        reset_z_min
    end

    function reset_z_min2(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 2, 1, 2);
        reset_z_min
    end

    function reset_z_min3(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 3, 2, 1);
        reset_z_min
    end

    function reset_z_min(source, eventdata)
        elev_min(curr_gui)   = elev_min_ref;
        depth_max = depth_max_ref;
        switch disp_type{curr_ax}
            case 'amp.'
                if (elev_min_ref < get(z_min_slide(curr_ax), 'min'))
                    set(z_min_slide(curr_ax), 'value', get(z_min_slide(curr_ax), 'min'))
                else
                    set(z_min_slide(curr_ax), 'value', elev_min_ref)
                end
                set(z_min_edit(curr_ax), 'string', sprintf('%4.0f', elev_min_ref))
            case 'flat'
                if (depth_min_ref < get(z_min_slide(curr_ax), 'min'))
                    set(z_min_slide(curr_ax), 'value', get(z_min_slide(curr_ax), 'min'))
                else
                    set(z_min_slide(curr_ax), 'value', depth_min_ref)
                end
                set(z_min_edit(curr_ax), 'string', sprintf('%4.0f', depth_max_ref))
        end
        update_z_range
    end

%% Reset maximum dB/x/y/z

    function reset_db_max1(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 2, 1, 2);
        reset_db_max
    end

    function reset_db_max2(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 3, 2, 1);
        reset_db_max
    end

    function reset_db_max(source, eventdata)
        if (db_max_ref(curr_rad) > get(cb_max_slide(curr_rad), 'max'))
            set(cb_max_slide(curr_rad), 'value', get(cb_max_slide(curr_rad), 'max'))
        else
            set(cb_max_slide(curr_rad), 'value', db_max_ref(curr_rad))
        end
        set(cb_max_edit(curr_rad), 'string', num2str(db_max_ref(curr_rad)))
        db_max(curr_rad)    = db_max_ref(curr_rad);
        update_db_range
    end

    function reset_dist_max1(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 2, 1, 2);
        reset_dist_max
    end

    function reset_dist_max2(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 3, 2, 1);
        reset_dist_max
    end

    function reset_dist_max(source, eventdata)
        if (dist_max_ref(curr_rad) > get(dist_max_slide(curr_rad), 'max'))
            set(dist_max_slide(curr_rad), 'value', get(dist_max_slide(curr_rad), 'max'))
        else
            set(dist_max_slide(curr_rad), 'value', dist_max_ref(curr_rad))
        end
        set(dist_max_edit(curr_rad), 'string', sprintf('%3.1f', dist_max_ref(curr_rad)))
        dist_max(curr_rad)  = dist_max_ref(curr_rad);
        update_dist_range
    end

    function reset_x_max(source, eventdata)
        [curr_gui, curr_ax] = deal(1);
        if (x_max_ref > get(x_max_slide, 'max'))
            set(x_max_slide, 'value', get(x_max_slide, 'max'))
        else
            set(x_max_slide, 'value', x_max_ref)
        end
        set(x_max_edit, 'string', sprintf('%4.1f', x_max_ref))
        x_max              = x_max_ref;
        update_x_range
    end

    function reset_y_max(source, eventdata)
        [curr_gui, curr_ax] = deal(1);
        if (y_max_ref > get(y_max_slide, 'max'))
            set(y_max_slide, 'value', get(y_max_slide, 'max'))
        else
            set(y_max_slide, 'value', y_max_ref)
        end
        set(y_max_edit, 'string', sprintf('%4.1f', y_max_ref))
        y_max              = y_max_ref;
        update_y_range
    end

    function reset_z_max1(source, eventdata)
        [curr_gui, curr_ax] = deal(1);
        reset_z_max
    end

    function reset_z_max2(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 2, 1, 2);
        reset_z_max
    end

    function reset_z_max3(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 3, 2, 1);
        reset_z_max
    end

    function reset_z_max(source, eventdata)
        elev_max(curr_gui)  = elev_max_ref;
        depth_min = depth_min_ref;
        switch disp_type{curr_ax}
            case 'amp.'
                if (elev_max_ref > get(z_max_slide(curr_ax), 'max'))
                    set(z_max_slide(curr_ax), 'value', get(z_max_slide(curr_ax), 'max'))
                else
                    set(z_max_slide(curr_ax), 'value', elev_max_ref)
                end
                set(z_max_edit(curr_ax), 'string', sprintf('%4.0f', elev_max_ref))
            case 'flat'
                if (depth_max_ref > get(z_max_slide(curr_ax), 'max'))
                    set(z_max_slide(curr_ax), 'value', get(z_max_slide(curr_ax), 'max'))
                else
                set(z_max_slide(curr_ax), 'value', depth_max_ref)
                end
                set(z_max_edit(curr_ax), 'string', sprintf('%4.0f', depth_min_ref))
        end
        update_z_range
    end

%% Reset all x/y/z

    function reset_xyz(source, eventdata)
        [curr_gui, curr_ax] = deal(1);
        reset_x_min
        reset_x_max
        reset_y_min
        reset_y_max
        reset_z_min
        reset_z_max
    end

    function reset_xz1(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 2, 1, 2);
        reset_xz
    end

    function reset_xz2(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 3, 2, 1);
        reset_xz
    end

    function reset_xz(source, eventdata)
        reset_dist_min
        reset_dist_max
        reset_z_min
        reset_z_max
    end

%% Update dB/x/y/z/elevation/depth range

    function update_db_range(source, eventdata)
        axes(ax(curr_ax))
        caxis([db_min(curr_rad) db_max(curr_rad)])
    end

    function update_dist_range(source, eventdata)
        axes(ax(curr_ax))
        xlim([dist_min(curr_rad) dist_max(curr_rad)])
    end

    function update_x_range(source, eventdata)
        axes(ax(curr_ax))
        xlim([x_min x_max])
    end

    function update_y_range(source, eventdata)
        axes(ax(curr_ax))
        ylim([y_min y_max])
    end

    function update_z_range(source, eventdata)
        axes(ax(curr_ax))
        switch disp_type{curr_ax}
            case 'amp.'
                if (curr_gui == 1)
                    zlim([elev_min(curr_gui) elev_max(curr_gui)])
                else
                    ylim([elev_min(curr_gui) elev_max(curr_gui)])
                end
            case 'flat'
                ylim([depth_min depth_max])
        end
        narrow_cb
    end

%% Adjust slider limits after panning or zooming

    function panzoom1
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 2, 1, 2);
        panzoom
    end

    function panzoom2
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 3, 2, 1);
        panzoom
    end

    function panzoom(source, eventdata)
        tmp1                = get(ax(curr_ax), 'xlim');
        if (tmp1(1) < dist_min_ref(curr_rad))
            reset_dist_min
        else
            if (tmp1(1) < get(dist_min_slide(curr_rad), 'min'))
                set(dist_min_slide(curr_rad), 'value', get(dist_min_slide(curr_rad), 'min'))
            else
                set(dist_min_slide(curr_rad), 'value', tmp1(1))
            end
            set(dist_min_edit(curr_rad), 'string', sprintf('%3.1f', tmp1(1)))
            dist_min(curr_rad) = tmp1(1);
        end
        if (tmp1(2) > dist_max_ref(curr_rad))
            reset_dist_max
        else
            if (tmp1(2) > get(dist_max_slide(curr_rad), 'max'))
                set(dist_max_slide(curr_rad), 'value', get(dist_max_slide(curr_rad), 'max'))
            else
                set(dist_max_slide(curr_rad), 'value', tmp1(2))
            end
            set(dist_max_edit(curr_rad), 'string', sprintf('%3.1f', tmp1(2)))
            dist_max(curr_rad) = tmp1(2);
        end
        tmp1                = get(ax(curr_ax), 'ylim');
        switch disp_type(curr_ax)
            case 'amp.'
                if (flat_done(curr_rad) && data_done(curr_rad))
                    tmp2    = [elev_min(curr_gui) elev_max(curr_gui)];
                end
                if (tmp1(1) < elev_min_ref)
                    reset_y_min
                else
                    if (tmp1(1) < get(z_min_slide(curr_ax), 'min'))
                        set(z_min_slide(curr_ax), 'value', get(z_min_slide(curr_ax), 'min'))
                    else
                        set(z_min_slide(curr_ax), 'value', tmp1(1))
                    end
                    set(z_min_edit(curr_ax), 'string', sprintf('%4.0f', tmp1(1)))
                    elev_min(curr_gui) = tmp1(1);
                end
                if (tmp1(2) > elev_max_ref)
                    reset_y_max
                else
                    if (tmp1(2) > get(z_max_slide(curr_ax), 'max'))
                        set(z_max_slide(curr_ax), 'value', get(z_max_slide(curr_ax), 'max'))
                    else
                        set(z_max_slide(curr_ax), 'value', tmp1(2))
                    end
                    set(z_max_edit(curr_ax), 'string', sprintf('%4.0f', tmp1(2)))
                    elev_max(curr_gui) = tmp1(2);
                end
                if (flat_done(curr_rad) && data_done(curr_rad))
                    depth_min = depth_min - (elev_max(curr_gui) - tmp2(2));
                    depth_max = depth_max - (elev_min(curr_gui) - tmp2(1));
                    if (depth_min < depth_min_ref)
                        depth_min = depth_min_ref;
                    end
                    if (depth_max > depth_max_ref)
                        depth_max = depth_max_ref;
                    end
                end
            case 'flat'
                tmp2        = [depth_min depth_max];
                if (tmp1(1) < depth_min_ref)
                    reset_y_min
                else
                    if (tmp1(1) < get(z_min_slide(curr_ax), 'min'))
                        set(z_min_slide(curr_ax), 'value', get(z_min_slide(curr_ax), 'min'))
                    elseif (tmp1(1) > get(z_min_slide(curr_ax), 'max'))
                        set(z_min_slide(curr_ax), 'value', get(z_min_slide(curr_ax), 'max'))
                    else
                        set(z_min_slide(curr_ax), 'value', tmp1(1))
                    end
                    set(z_min_edit(curr_ax), 'string', sprintf('%4.0f', tmp1(1)))
                    depth_min ...
                            = tmp1(1);
                end
                if (tmp1(2) > depth_max_ref)
                    reset_y_max
                else
                    if (tmp1(2) < get(z_max_slide(curr_ax), 'min'))
                        set(z_max_slide(curr_ax), 'value', get(z_max_slide(curr_ax), 'min'))
                    elseif (tmp1(2) > get(y_max_slide, 'max'))
                        set(z_max_slide(curr_ax), 'value', get(z_max_slide(curr_ax), 'max'))
                    else
                        set(z_max_slide(curr_ax), 'value', tmp1(2))
                    end
                    set(z_max_edit(curr_ax), 'string', sprintf('%4.0f', tmp1(2)))
                    depth_max ...
                            = tmp1(2);
                end
                elev_min(curr_gui) = elev_min(curr_gui) - (depth_max - tmp2(2));
                elev_max(curr_gui) = elev_max(curr_gui) - (depth_min - tmp2(1));
                if (elev_min(curr_gui) < elev_min_ref)
                    elev_min(curr_gui) = elev_min_ref;
                end
                if (elev_max(curr_gui) > elev_max_ref)
                    elev_max(curr_gui) = elev_max_ref;
                end
        end
        narrow_cb
    end

%% Plot data (amplitude)

    function plot_db(source, eventdata)
        if ~data_done(curr_rad)
            set(status_box(2), 'string', 'Data not loaded yet.')
            return
        end
        if (logical(p_data(curr_rad)) && ishandle(p_data(curr_rad)))
            delete(p_data(curr_rad))
        end
        if (logical(p_bedflat(curr_rad)) && ishandle(p_bedflat(curr_rad)))
            set(p_bedflat(curr_rad), 'visible', 'off')
        end
        axes(ax(curr_ax))
        if strcmp(disp_type{2}, disp_type{3})
            linkaxes(ax(2:3), 'y')
            h_link1         = linkprop(z_min_slide(2:3), {'value' 'min' 'max'});
            h_link2         = linkprop(z_max_slide(2:3), {'value' 'min' 'max'});
            h_link3         = linkprop(z_min_edit(2:3), 'string');
            h_link4         = linkprop(z_max_edit(2:3), 'string');
        else
            linkaxes(ax(2:3), 'off')
            removeprop(h_link1, 'value');
            removeprop(h_link1, 'min');
            removeprop(h_link1, 'max');
            removeprop(h_link2, 'value');
            removeprop(h_link2, 'min');
            removeprop(h_link2, 'max');
            removeprop(h_link3, 'string');
            removeprop(h_link4, 'string');
        end
        axis xy
        set(z_min_slide(curr_ax), 'min', elev_min_ref, 'max', elev_max_ref, 'value', (elev_max_ref - (elev_max(2) - elev_min_ref)))
        set(z_max_slide(curr_ax), 'min', elev_min_ref, 'max', elev_max_ref, 'value', (elev_max_ref - (elev_min(2) - elev_min_ref)))
        set(z_min_edit(curr_ax), 'string', sprintf('%4.0f', elev_max(2)))
        set(z_max_edit(curr_ax), 'string', sprintf('%4.0f', elev_min(2)))
        ylim([elev_min(2) elev_max(2)])
        if data_done(curr_rad)
            p_data(curr_rad)= imagesc(pk{curr_rad}.dist_lin(ind_decim{2, curr_rad}), elev{curr_rad}, amp_mean{curr_rad}, [db_min(curr_rad) db_max(curr_rad)]);
        end
        disp_type{curr_ax}  = 'amp.';
        narrow_cb
        show_data
        show_core
        show_pk
        show_int
    end

%% Plot layer-flattened radargram

    function plot_flat(source, eventdata)
        if (~flat_done(curr_rad) || ~data_done(curr_rad))
            set(disp_group(curr_rad), 'selectedobject', disp_check(curr_rad, 1))
            disp_type{curr_ax} = 'amp.';
            plot_db
            return
        end
        if (logical(p_data(curr_rad)) && ishandle(p_data(curr_rad))) % get rid of old plotted data
            delete(p_data(curr_rad))
        end
        if (logical(p_bedflat(curr_rad)) && ishandle(p_bedflat(curr_rad)))
            delete(p_bedflat(curr_rad))
        end
        if (logical(p_bed(2, curr_rad)) && ishandle(p_bed(2, curr_rad)))
            set(p_bed(curr_rad), 'visible', 'off')
        end
        if (logical(p_surf(2, curr_rad)) && ishandle(p_surf(2, curr_rad)))
            set(p_surf(2, curr_rad), 'visible', 'off')
        end
        axes(ax(curr_ax))
        if (get(cmap_list, 'value') ~= 1)
            set(cmap_list, 'value', 1)
            change_cmap
        end
        if strcmp(disp_type{2}, disp_type{3})
            linkaxes(ax(2:3), 'y')
            h_link1         = linkprop(z_min_slide(2:3), {'value' 'min' 'max'});
            h_link2         = linkprop(z_max_slide(2:3), {'value' 'min' 'max'});
            h_link3         = linkprop(z_min_edit(2:3), 'string');
            h_link4         = linkprop(z_max_edit(2:3), 'string');
        else
            linkaxes(ax(2:3), 'off')
            removeprop(h_link1, 'value');
            removeprop(h_link1, 'min');
            removeprop(h_link1, 'max');
            removeprop(h_link2, 'value');
            removeprop(h_link2, 'min');
            removeprop(h_link2, 'max');
            removeprop(h_link3, 'string');
            removeprop(h_link4, 'string');
        end
        axis ij
        set(z_min_slide(curr_ax), 'min', depth_min_ref, 'max', depth_max_ref, 'value', (depth_max_ref - (depth_max - depth_min_ref)))
        set(z_max_slide(curr_ax), 'min', depth_min_ref, 'max', depth_max_ref, 'value', (depth_max_ref - (depth_min - depth_min_ref)))
        set(z_min_edit(curr_ax), 'string', sprintf('%4.0f', depth_max))
        set(z_max_edit(curr_ax), 'string', sprintf('%4.0f', depth_min))
        ylim([depth_min depth_max])
        p_data(curr_rad)    = imagesc(pk{curr_rad}.dist_lin(ind_decim{2, curr_rad}), depth{curr_rad}, amp_flat{curr_rad}, [db_min(curr_rad) db_max(curr_rad)]);
        if (bed_avail(curr_rad) && any(~isnan(depth_bed_flat{curr_rad})))
            if any(isnan(depth_bed_flat{curr_rad}))
                p_bedflat(curr_rad) ...
                            = plot(pk{curr_rad}.dist_lin(ind_decim{2, curr_rad}(~isnan(depth_bed_flat{curr_rad}))), depth_bed_flat{curr_rad}(~isnan(depth_bed_flat{curr_rad})), 'g.', 'markersize', 12);
            else
                p_bedflat(curr_rad) ...
                            = plot(pk.dist_lin(ind_decim{2, curr_rad}), depth_bed_flat{curr_rad}, 'g--', 'linewidth', 2);
            end
        end
        disp_type{curr_ax}  = 'flat';
        narrow_cb
        show_data
        show_core
        show_pk
        show_int
    end

%% Show radar data

    function show_data1(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 2, 1, 2);
       show_data
    end

    function show_data2(source, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 3, 2, 1);
       show_data
    end

    function show_data(source, eventdata)
        if data_done(curr_rad)
            if (get(data_check(curr_rad), 'value') && logical(p_data(curr_rad)) && ishandle(p_data(curr_rad)))
                set(p_data(curr_rad), 'visible', 'on')
                set(p_surf(curr_gui, curr_rad), 'visible', 'on')
                set(p_bed(curr_gui, curr_rad), 'visible', 'on')
            elseif (logical(p_data(curr_rad)) && ishandle(p_data(curr_rad)))
                set(p_data(curr_rad), 'visible', 'off')
                set(p_surf(curr_gui, curr_rad), 'visible', 'off')
                set(p_bed(curr_gui, curr_rad), 'visible', 'off')
            end
        elseif get(data_check(curr_rad), 'value')
            set(data_check(curr_rad), 'value', 0)
        end
    end

%% Show picked layers

    function show_pk1(source, eventdata)
       [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(1, 1, 1, 2);
       show_pk
    end

    function show_pk2(source, eventdata)
       [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(1, 1, 2, 1);
       show_pk
    end

    function show_pk3(source, eventdata)
       [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 2, 1, 2);
       show_pk
    end

    function show_pk4(source, eventdata)
       [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 3, 2, 1);
       show_pk
    end

    function show_pk(source, eventdata)
        if pk_done(curr_rad)
            if get(pk_check(curr_gui, curr_rad), 'value')
                switch disp_type{curr_ax}
                    case 'amp.'
                        if (any(p_pk{curr_gui, curr_rad}) && any(ishandle(p_pk{curr_gui, curr_rad})))
                            set(p_pk{curr_gui, curr_rad}(logical(p_pk{curr_gui, curr_rad}) & ishandle(p_pk{curr_gui, curr_rad})), 'visible', 'on')
                            uistack(p_pk{curr_gui, curr_rad}(logical(p_pk{curr_gui, curr_rad}) & ishandle(p_pk{curr_gui, curr_rad})), 'top')
                            if (any(p_pkflat{curr_rad}) && any(ishandle(p_pkflat{curr_rad})))
                                set(p_pkflat{curr_rad}(logical(p_pkflat{curr_rad}) & ishandle(p_pkflat{curr_rad})), 'visible', 'off')
                            end
                        end
                        if (logical(p_surf(curr_gui, curr_rad)) && ishandle(p_surf(curr_gui, curr_rad)))
                            set(p_surf(curr_gui, curr_rad), 'visible', 'on')
                        end
                        if (logical(p_bed(curr_gui, curr_rad)) && ishandle(p_bed(curr_gui, curr_rad)))
                            set(p_bed(curr_gui, curr_rad), 'visible', 'on')
                        end
                    case 'flat'
                        if (any(p_pkflat{curr_rad}) && any(ishandle(p_pkflat{curr_rad})))
                            set(p_pkflat{curr_rad}(logical(p_pkflat{curr_rad}) & ishandle(p_pkflat{curr_rad})), 'visible', 'on')
                            uistack(p_pkflat{curr_rad}(logical(p_pkflat{curr_rad}) & ishandle(p_pkflat{curr_rad})), 'top')
                        end
                        if (any(p_pk{curr_gui, curr_rad}) && any(ishandle(p_pk{curr_gui, curr_rad})))
                            set(p_pk{curr_gui, curr_rad}(logical(p_pk{curr_gui, curr_rad}) & ishandle(p_pk{curr_gui, curr_rad})), 'visible', 'off')
                        end
                end
            else
                if (any(p_pk{curr_gui, curr_rad}) && any(ishandle(p_pk{curr_gui, curr_rad})))
                    set(p_pk{curr_gui, curr_rad}(logical(p_pk{curr_gui, curr_rad}) & ishandle(p_pk{curr_gui, curr_rad})), 'visible', 'off')
                end
                if (flat_done(curr_rad) && data_done(curr_rad))
                    set(p_pkflat{curr_rad}(logical(p_pkflat{curr_rad}) & ishandle(p_pkflat{curr_rad})), 'visible', 'off')
                end
                if (logical(p_surf(curr_gui, curr_rad)) && ishandle(p_surf(curr_gui, curr_rad)))
                    set(p_surf(curr_gui, curr_rad), 'visible', 'off')
                end
                if (logical(p_bed(curr_gui, curr_rad)) && ishandle(p_bed(curr_gui, curr_rad)))
                    set(p_bed(curr_gui, curr_rad), 'visible', 'off')
                end
            end
        elseif get(pk_check(curr_gui, curr_rad), 'value')
            set(pk_check(curr_gui, curr_rad), 'value', 0)
        end
    end

%% Show intersections

    function show_int1(source, eventdata)
       [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(1, 1, 1, 2);
       show_int
    end

    function show_int2(source, eventdata)
       [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 2, 1, 2);
       show_int
    end

    function show_int3(source, eventdata)
       [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 3, 2, 1);
       show_int
    end

    function show_int(source, eventdata)
        if pk_done(2)
            if get(int_check(curr_ax), 'value')
                switch disp_type{curr_ax}
                    case 'amp.'
                        if (any(p_int1{1, curr_ax}) && any(ishandle(p_int1{1, curr_ax})))
                            set(p_int1{1, curr_ax}(logical(p_int1{1, curr_ax}) & ishandle(p_int1{1, curr_ax})), 'visible', 'on')
                            uistack(p_int1{1, curr_ax}(logical(p_int1{1, curr_ax}) & ishandle(p_int1{1, curr_ax})), 'top')
                        end
                        if (any(p_int2{1, curr_rad}) && any(ishandle(p_int2{1, curr_rad})))
                            set(p_int2{1, curr_rad}(logical(p_int2{1, curr_rad}) & ishandle(p_int2{1, curr_rad})), 'visible', 'on')
                            uistack(p_int2{1, curr_rad}(logical(p_int2{1, curr_rad}) & ishandle(p_int2{1, curr_rad})), 'top')
                        end
                        if (any(p_int1{2, curr_ax}) && any(ishandle(p_int1{2, curr_ax})))
                            set(p_int1{2, curr_ax}(logical(p_int1{2, curr_ax}) & ishandle(p_int1{2, curr_ax})), 'visible', 'off')
                        end
                        if (any(p_int2{2, curr_rad}) && any(ishandle(p_int2{1, curr_rad})))
                            set(p_int2{2, curr_rad}(logical(p_int2{2, jj}) & ishandle(p_int2{2, curr_rad})), 'visible', 'off')
                        end
                case 'flat'
                        if (any(p_int1{2, curr_ax}) && any(ishandle(p_int1{2, curr_ax})))
                            set(p_int1{2, curr_ax}(logical(p_int1{2, curr_ax}) & ishandle(p_int1{2, curr_ax})), 'visible', 'on')
                            uistack(p_int1{2, curr_ax}(logical(p_int1{2, curr_ax}) & ishandle(p_int1{2, curr_ax})), 'top')
                        end
                        if (any(p_int2{2, curr_rad}) && any(ishandle(p_int2{1, curr_rad})))
                            set(p_int2{2, curr_rad}(logical(p_int2{2, jj}) & ishandle(p_int2{2, curr_rad})), 'visible', 'on')
                            uistack(p_int2{2, curr_rad}(logical(p_int2{2, curr_rad}) & ishandle(p_int2{2, curr_rad})), 'top')
                        end
                        if (any(p_int1{1, curr_ax}) && any(ishandle(p_int1{1, curr_ax})))
                            set(p_int1{1, curr_ax}(logical(p_int1{1, curr_ax}) & ishandle(p_int1{1, curr_ax})), 'visible', 'off')
                        end
                        if (any(p_int2{1, curr_rad}) && any(ishandle(p_int2{1, curr_rad})))
                            set(p_int2{1, curr_rad}(logical(p_int2{1, curr_rad}) & ishandle(p_int2{1, curr_rad})), 'visible', 'off')
                        end
                end
            else
                for ii = 1:2
                    if (any(p_int1{ii, curr_ax}) && any(ishandle(p_int1{ii, curr_ax})))
                        set(p_int1{ii, curr_ax}(logical(p_int1{ii, curr_ax}) & ishandle(p_int1{ii, curr_ax})), 'visible', 'off')
                    end
                    for jj = 1:2
                        if (any(p_int2{ii, jj}) && any(ishandle(p_int2{ii, jj})))
                            set(p_int2{ii, jj}(logical(p_int2{ii, jj}) & ishandle(p_int2{ii, jj})), 'visible', 'off')
                        end
                    end
                end
            end
        elseif get(int_check(curr_ax), 'value')
            set(int_check(curr_ax), 'value', 0)
        end
    end

%% Show core intersections

    function show_core1(source, eventdata)
       [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(1, 1, 1, 2);
       show_core
    end

    function show_core2(source, eventdata)
       [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 2, 1, 2);
       show_core
    end

    function show_core3(source, eventdata)
       [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 3, 2, 1);
       show_core
    end

    function show_core(source, eventdata)
        if (core_done && ~isempty(ind_int_core{curr_rad}))
            if get(core_check(curr_ax), 'value')
                switch disp_type{curr_ax}
                    case 'amp.'
                        if (any(p_core{curr_gui, curr_rad}) && any(ishandle(p_core{curr_gui, curr_rad})))
                            set(p_core{curr_gui, curr_rad}(logical(p_core{curr_gui, curr_rad}) & ishandle(p_core{curr_gui, curr_rad})), 'visible', 'on')
                            uistack(p_core{curr_gui, curr_rad}(logical(p_core{curr_gui, curr_rad}) & ishandle(p_core{curr_gui, curr_rad})), 'top')
                        end
                        if (any(p_corename{curr_gui, curr_rad}) && any(ishandle(p_corename{curr_gui, curr_rad})))
                            set(p_corename{curr_gui, curr_rad}(logical(p_corename{curr_gui, curr_rad}) & ishandle(p_corename{curr_gui, curr_rad})), 'visible', 'on')
                            uistack(p_corename{curr_gui, curr_rad}(logical(p_corename{curr_gui, curr_rad}) & ishandle(p_corename{curr_gui, curr_rad})), 'top')
                        end
                        if (any(p_coreflat{curr_rad}) && any(ishandle(p_coreflat{curr_rad})))
                            set(p_coreflat{curr_rad}(logical(p_coreflat{curr_rad}) & ishandle(p_coreflat{curr_rad})), 'visible', 'off')
                        end
                        if (any(p_corenameflat{curr_rad}) && any(ishandle(p_corenameflat{curr_rad})))
                            set(p_corenameflat{curr_rad}(logical(p_corenameflat{curr_rad}) & ishandle(p_corenameflat{curr_rad})), 'visible', 'off')
                        end
                    case 'flat'
                        if (any(p_coreflat{curr_rad}) && any(ishandle(p_coreflat{curr_rad})))
                            set(p_coreflat{curr_rad}(logical(p_coreflat{curr_rad}) & ishandle(p_coreflat{curr_rad})), 'visible', 'on')
                        end
                        if (any(p_corenameflat{curr_rad}) && any(ishandle(p_corenameflat{curr_rad})))
                            set(p_corenameflat{curr_rad}(logical(p_corenameflat{curr_rad}) & ishandle(p_corenameflat{curr_rad})), 'visible', 'on')
                        end
                        if (any(p_core{curr_gui, curr_rad}) && any(ishandle(p_core{curr_gui, curr_rad})))
                            set(p_core{curr_gui, curr_rad}(logical(p_core{curr_gui, curr_rad}) & ishandle(p_core{curr_gui, curr_rad})), 'visible', 'off')
                        end
                        if (any(p_corename{curr_gui, curr_rad}) && any(ishandle(p_corename{curr_gui, curr_rad})))
                            set(p_corename{curr_gui, curr_rad}(logical(p_corename{curr_gui, curr_rad}) & ishandle(p_corename{curr_gui, curr_rad})), 'visible', 'off')
                        end
                end
            else
                if (any(p_core{curr_gui, curr_rad}) && any(ishandle(p_core{curr_gui, curr_rad})))
                    set(p_core{curr_gui, curr_rad}(logical(p_core{curr_gui, curr_rad}) & ishandle(p_core{curr_gui, curr_rad})), 'visible', 'off')
                end
                if (any(p_corename{curr_gui, curr_rad}) && any(ishandle(p_corename{curr_gui, curr_rad})))
                    set(p_corename{curr_gui, curr_rad}(logical(p_corename{curr_gui, curr_rad}) & ishandle(p_corename{curr_gui, curr_rad})), 'visible', 'off')
                end
                if (any(p_coreflat{curr_rad}) && any(ishandle(p_coreflat{curr_rad})))
                    set(p_coreflat{curr_rad}(logical(p_coreflat{curr_rad}) & ishandle(p_coreflat{curr_rad})), 'visible', 'off')
                end
                if (any(p_corenameflat{curr_rad}) && any(ishandle(p_corenameflat{curr_rad})))
                    set(p_corenameflat{curr_rad}(logical(p_corenameflat{curr_rad}) & ishandle(p_corenameflat{curr_rad})), 'visible', 'off')
                end
            end
        elseif get(core_check(curr_ax), 'value')
            set(core_check(curr_ax), 'value', 0)
        end
    end

%% Adjust number of indices to display

    function adj_decim1(source, eventdata)
       [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(1, 1, 1, 2);
        adj_decim
    end

    function adj_decim2(source, eventdata)
       [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(1, 1, 2, 1);
        adj_decim
    end

    function adj_decim3(source, eventdata)
       [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 2, 1, 2);
        adj_decim
    end

    function adj_decim4(source, eventdata)
       [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 3, 2, 1);
        adj_decim
    end

    function adj_decim(source, eventdata)
        decim(curr_gui, curr_rad) ...
                            = abs(round(str2double(get(decim_edit(curr_gui, curr_rad), 'string'))));
        if pk_done(curr_rad)
            set(status_box(curr_gui), 'string', 'Updating picks to new decimation...')
            pause(0.1)
            if (decim(curr_gui, curr_rad) > 1)
                ind_decim{curr_gui, curr_rad} ...
                            = (1 + ceil(decim(curr_gui, curr_rad) / 2)):decim(curr_gui, curr_rad):(pk{curr_rad}.num_trace_tot - ceil(decim(curr_gui, curr_rad) / 2));
            else
                ind_decim{curr_gui, curr_rad} ...
                            = 1:pk{curr_rad}.num_trace_tot;
            end
            num_decim(curr_gui, curr_rad) ...
                            = length(ind_decim{curr_gui, curr_rad});
            if (any(p_pk{curr_gui, curr_rad}) && any(ishandle(p_pk{curr_gui, curr_rad})))
                delete(p_pk{curr_gui, curr_rad}(logical(p_pk{curr_gui, curr_rad}) & ishandle(p_pk{curr_gui, curr_rad})))
            end
            layer_str{curr_rad} ...
                            = num2cell(1:pk{curr_rad}.num_layer);
            switch curr_gui
                case 1
                    for ii = 1:pk{curr_rad}.num_layer
                        if ~any(~isnan(pk{curr_rad}.elev_smooth(ii, ind_decim{curr_gui, curr_rad})))
                            p_pk{curr_gui, curr_rad}(ii) = plot(0, 0, 'w.', 'markersize', 1);
                            layer_str{curr_rad}{ii} = [num2str(ii) ' H'];
                        else
                            p_pk{curr_gui, curr_rad}(ii) = plot3(pk{curr_rad}.x(ind_decim{curr_gui, curr_rad}(~isnan(pk{curr_rad}.elev_smooth(ii, ind_decim{curr_gui, curr_rad})))), ...
                                                                 pk{curr_rad}.y(ind_decim{curr_gui, curr_rad}(~isnan(pk{curr_rad}.elev_smooth(ii, ind_decim{curr_gui, curr_rad})))), ...
                                                                 pk{curr_rad}.elev_smooth(ii, ind_decim{curr_gui, curr_rad}(~isnan(pk{curr_rad}.elev_smooth(ii, ind_decim{curr_gui, curr_rad})))), ...
                                                                 '.', 'color', colors{curr_rad}(ii, :), 'markersize', 12, 'visible', 'off');
                        end
                    end
                case 2
                    for ii = 1:pk{curr_rad}.num_layer
                        if ~any(~isnan(pk{curr_rad}.elev_smooth(ii, ind_decim{curr_gui, curr_rad})))
                            p_pk{curr_gui, curr_rad}(ii) = plot(0, 0, 'w.', 'markersize', 1);
                            layer_str{curr_rad}{ii} = [num2str(ii) ' H'];
                        else
                            p_pk{curr_gui, curr_rad}(ii) = plot(pk{curr_rad}.dist_lin(ind_decim{curr_gui, curr_rad}(~isnan(pk{curr_rad}.elev_smooth(ii, ind_decim{curr_gui, curr_rad})))), ...
                                                                pk{curr_rad}.elev_smooth(ii, ind_decim{curr_gui, curr_rad}(~isnan(pk{curr_rad}.elev_smooth(ii, ind_decim{curr_gui, curr_rad})))), ...
                                                                '.', 'color', colors{curr_rad}(ii, :), 'markersize', 12, 'visible', 'off');
                        end
                    end
            end
            % calculate bed depth (for flattened projection)
            if (surf_avail(curr_rad) && bed_avail(curr_rad))
                depth_bed{curr_rad} ...
                            = pk{curr_rad}.elev_surf(ind_decim{2, curr_rad}) - pk{curr_rad}.elev_bed(ind_decim{2, curr_rad});
            end
            if curr_layer(curr_rad)
                set(layer_list(:, curr_rad), 'string', layer_str{curr_rad}, 'value', curr_layer(curr_rad))
            else
                set(layer_list(:, curr_rad), 'string', layer_str{curr_rad}, 'value', 1)
            end
            show_pk
        end
        if ((curr_gui == 2) && data_done(curr_rad))
            set(status_box(curr_gui), 'string', 'Updating data to new decimation...')
            pause(0.1)
            load_data_breakout
        end
        set(decim_edit(curr_gui, curr_rad), 'string', num2str(decim(curr_gui, curr_rad)))
        set(status_box(curr_gui), 'string', ['Decimation number set to 1/' num2str(decim(curr_gui, curr_rad)) ' indice(s).'])
    end

%% Change colormap

    function change_cmap(source, eventdata)
        axes(ax(curr_ax))
        colormap(cmaps{get(cmap_list, 'value')})
    end

%% Change intersection

    function change_int(source, eventdata)
        curr_gui            = 2;
        curr_int            = get(intnum_list, 'value');
        if (num_int && (curr_ind_int(1) ~= 0))
            for ii = 1:2
                for jj = 1:3
                    if (any(p_int1{ii, jj}) && any(ishandle(p_int1{ii, jj})))
                        set(p_int1{ii, jj}(logical(p_int1{ii, jj}) & ishandle(p_int1{ii, jj})), 'linewidth', 2)
                    end
                    if (logical(p_int1{ii, jj}(curr_int)) && ishandle(p_int1{ii, jj}(curr_int)))
                        set(p_int1{ii, jj}(curr_int), 'linewidth', 4)
                    end
                end
            end
            [dist_min(1), dist_max(1), dist_min(2), dist_max(2)] ...
                            = deal((pk{1}.dist_lin(curr_ind_int(curr_int, 1)) - 10), (pk{1}.dist_lin(curr_ind_int(curr_int, 1)) + 10), ...
                                   (pk{2}.dist_lin(curr_ind_int(curr_int, 2)) - 10), (pk{2}.dist_lin(curr_ind_int(curr_int, 2)) + 10));
            if (dist_min(1) < dist_min_ref(1))
                dist_min(1) = dist_min_ref(1);
            end
            if (dist_min(2) < dist_min_ref(2))
                dist_min(2) = dist_min_ref(2);
            end
            if (dist_max(1) > dist_max_ref(1))
                dist_max(1) = dist_max_ref(1);
            end
            if (dist_max(2) > dist_max_ref(2))
                dist_max(2) = dist_max_ref(2);
            end
            axes(ax(2))
            xlim([dist_min(1) dist_max(1)])
            narrow_cb
            axes(ax(3))
            xlim([dist_min(2) dist_max(2)])
            narrow_cb
            set(dist_min_slide(1), 'value', dist_min(1))
            set(dist_min_slide(2), 'value', dist_min(2))
            set(dist_max_slide(1), 'value', dist_max(1))
            set(dist_max_slide(2), 'value', dist_max(2))
            set(dist_min_edit(1), 'string', sprintf('%3.0f', dist_min(1)))
            set(dist_min_edit(2), 'string', sprintf('%3.0f', dist_min(2)))
            set(dist_max_edit(1), 'string', sprintf('%3.0f', dist_max(1)))
            set(dist_max_edit(2), 'string', sprintf('%3.0f', dist_max(2)))
        end
    end

%% Toggle gridlines

    function toggle_grid1(source, eventdata)
        [curr_gui, curr_ax] = deal(1);
        toggle_grid
    end

    function toggle_grid2(source, eventdata)
        [curr_gui, curr_ax] = deal(2);
        toggle_grid
    end

    function toggle_grid3(source, eventdata)
        [curr_gui, curr_ax] = deal(2, 3);
        toggle_grid
    end

    function toggle_grid
        if get(grid_check(curr_ax), 'value')
            axes(ax(curr_ax))
            grid on
        else
            grid off
        end
    end

%% Narrow color axis to +/- 2 standard deviations of current mean value

    function narrow_cb1(~, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 2, 1, 2);
        narrow_cb
    end

    function narrow_cb2(~, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 3, 2, 1);
        narrow_cb
    end

    function narrow_cb(source, eventdata)
        if (get(cbfix_check2(curr_rad), 'value') && data_done(curr_rad))
            axes(ax(curr_ax))
            tmp4            = zeros(2);
            tmp4(1, :)      = interp1(pk{curr_rad}.dist_lin(ind_decim{2, curr_rad}), 1:num_decim(2, curr_rad), [dist_min(curr_rad) dist_max(curr_rad)], 'nearest', 'extrap');
            switch disp_type{curr_ax}
                case 'amp.'
                    tmp4(2, :) = interp1(elev{curr_rad}, 1:num_sample(curr_rad), [elev_min(curr_gui) elev_max(curr_gui)], 'nearest', 'extrap');
                    tmp4(2, :) = flipud(tmp4(2, :));
                    tmp4    = amp_mean{curr_rad}(tmp4(2, 1):tmp4(2, 2), tmp4(1, 1):tmp4(1, 2));
                case 'flat'
                    tmp4(2, :) = interp1(depth{curr_rad}, 1:length(depth{curr_rad}), [depth_min depth_max], 'nearest', 'extrap');
                    tmp4    = amp_flat{curr_rad}(tmp4(2, 1):tmp4(2, 2), tmp4(1, 1):tmp4(1, 2));
            end
            [tmp5(1), tmp5(2)] ...
                            = deal(nanmean(tmp4(~isinf(tmp4))), nanstd(tmp4(~isinf(tmp4))));
            if any(isnan(tmp5))
                return
            end
            tmp4            = zeros(1, 2);
            if ((tmp5(1) - (2 * tmp5(2))) < db_min_ref(curr_rad))
                tmp4(1)     = db_min_ref(curr_rad);
            else
                tmp4(1)     = tmp5(1) - (2 * tmp5(2));
            end
            if ((tmp5(1) + (2 * tmp5(2))) > db_max_ref(curr_rad))
                tmp4(2)     = db_max_ref(curr_rad);
            else
                tmp4(2)     = tmp5(1) + (2 * tmp5(2));
            end
            [db_min(curr_rad), db_max(curr_rad)] ...
                            = deal(tmp4(1), tmp4(2));
            if (db_min(curr_rad) < get(cb_min_slide(curr_rad), 'min'))
                set(cb_min_slide(curr_rad), 'value', get(cb_min_slide(curr_rad), 'min'))
            else
                set(cb_min_slide(curr_rad), 'value', db_min(curr_rad))
            end
            if (db_max(curr_rad) > get(cb_max_slide(curr_rad), 'max'))
                set(cb_max_slide(curr_rad), 'value', get(cb_max_slide(curr_rad), 'max'))
            else
                set(cb_max_slide(curr_rad), 'value', db_max(curr_rad))
            end
            set(cb_min_edit(curr_rad), 'string', sprintf('%3.0f', db_min(curr_rad)))
            set(cb_max_edit(curr_rad), 'string', sprintf('%3.0f', db_max(curr_rad)))
            caxis([db_min(curr_rad) db_max(curr_rad)])
        end
    end

%% Switch display type

    function disp_radio1(~, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 2, 1, 2);
        disp_type{curr_ax}  = get(eventdata.NewValue, 'string');
        disp_radio
    end

    function disp_radio2(~, eventdata)
        [curr_gui, curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 3, 2, 1);
        disp_type{curr_ax}  = get(eventdata.NewValue, 'string');
        disp_radio
    end

    function disp_radio(source, eventdata)
        switch disp_type{curr_ax}
            case 'amp.'
                plot_db
            case 'flat'
                if (flat_done(curr_rad) && data_done(curr_rad))
                    plot_flat
                else
                    disp_type{curr_ax} = 'amp.';
                    set(disp_group(curr_rad), 'selectedobject', disp_check(curr_rad, 1));
                end
        end
    end

%% Switch viewing dimension

    function choose_dim(~, eventdata)
        [curr_gui, curr_ax] = deal(1);
        curr_dim            = get(eventdata.NewValue, 'string');
        axes(ax(1))
        switch curr_dim
            case '2D'
                [curr_az3, curr_el3] ...
                            = view;
                if any([curr_az2 curr_el2])
                    view(curr_az2, curr_el2)
                else
                    view(2)
                end
            case '3D'
                [curr_az2, curr_el2] ...
                            = view;
                view(curr_az3, curr_el3)
        end
    end

%% Keyboard shortcuts for various functions

    function keypress1(~, eventdata)
        [curr_gui, curr_ax] = deal(1);
        switch eventdata.Key
            case '1'
                [curr_rad, curr_rad_alt] ...
                            = deal(1, 2);
                if get(pk_check(1, 1), 'value')
                    set(pk_check(1, 1), 'value', 0)
                else
                    set(pk_check(1, 1), 'value', 1)
                end
                show_pk
            case '2'
                [curr_rad, curr_rad_alt] ...
                            = deal(2, 1);
                if get(pk_check(1, 2), 'value')
                    set(pk_check(1, 2), 'value', 0)
                else
                    set(pk_check(1, 2), 'value', 1)
                end
                show_pk
            case '3'
                if get(int_check(1), 'value')
                    set(int_check(1), 'value', 0)
                else
                    set(int_check(1), 'value', 1)
                end
                show_int1
            case '4'
                if get(core_check(1), 'value')
                    set(core_check(1), 'value', 0)
                else
                    set(core_check(1), 'value', 1)
                end
                show_core1
            case 'c'
                load_core
            case 'e'
                reset_xyz
            case 'g'
                if get(grid_check(1), 'value')
                    set(grid_check(1), 'value', 0)
                else
                    set(grid_check(1), 'value', 1)
                end
                toggle_grid
            case 'l'
                pk_last
            case 'n'
                pk_next
            case 'p'
                load_pk1
            case 'r'
                load_int
            case 't'
                misctest
            case 'v'
                pk_save
            case 'x'
                if get(xfix_check, 'value')
                    set(xfix_check, 'value', 0)
                else
                    set(xfix_check, 'value', 1)
                end
            case 'y'
                if get(yfix_check, 'value')
                    set(yfix_check, 'value', 0)
                else
                    set(yfix_check, 'value', 1)
                end
            case 'z'
                if get(zfix_check(1), 'value')
                    set(zfix_check(1), 'value', 0)
                else
                    set(zfix_check(1), 'value', 1)
                end
            case 'spacebar'
                axes(ax(1))
                switch curr_dim
                    case '2D'
                        set(dim_group, 'selectedobject', dim_check(2))
                        curr_dim ...
                            = '3D';
                        [curr_az2, curr_el2] ...
                            = view;
                        view(curr_az3, curr_el3)
                    case '3D'
                        set(dim_group, 'selectedobject', dim_check(1))
                        curr_dim ...
                            = '2D';
                        [curr_az3, curr_el3] ...
                            = view;
                        if any([curr_az2 curr_el2])
                            view(curr_az2, curr_el2)
                        else
                            view(2)
                        end
                end
        end
    end

    function keypress2(~, eventdata)
        curr_gui            = 2;
        switch eventdata.Key
            case '1'
                [curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 1, 2);
                axes(ax(curr_ax))
            case '2'
                [curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(3, 2, 1);
                axes(ax(curr_ax))
            case '3'
                [curr_rad, curr_rad_alt] ...
                            = deal(1, 2);
                if get(pk_check(curr_gui, curr_rad), 'value')
                    set(pk_check(curr_gui, curr_rad), 'value', 0)
                else
                    set(pk_check(curr_gui, curr_rad), 'value', 1)
                end
                show_pk
            case '4'
                [curr_rad, curr_rad_alt] ...
                            = deal(2, 1);
                if get(pk_check(curr_gui, curr_rad), 'value')
                    set(pk_check(curr_gui, curr_rad), 'value', 0)
                else
                    set(pk_check(curr_gui, curr_rad), 'value', 1)
                end
                show_pk
            case '5'
                [curr_rad, curr_rad_alt] ...
                            = deal(1, 2);
                if get(data_check(curr_rad), 'value')
                    set(data_check(curr_rad), 'value', 0)
                else
                    set(data_check(curr_rad), 'value', 1)
                end
                show_data
            case '6'
                [curr_rad, curr_rad_alt] ...
                            = deal(2, 1);
                if get(data_check(curr_rad), 'value')
                    set(data_check(curr_rad), 'value', 0)
                else
                    set(data_check(curr_rad), 'value', 1)
                end
                show_data
            case '7'
                [curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(2, 1, 2);
                if get(int_check(curr_ax), 'value')
                    set(int_check(curr_ax), 'value', 0)
                else
                    set(int_check(curr_ax), 'value', 1)
                end
                show_int
            case '8'
                [curr_ax, curr_rad, curr_rad_alt] ...
                            = deal(3, 2, 1);
                if get(int_check(curr_ax), 'value')
                    set(int_check(curr_ax), 'value', 0)
                else
                    set(int_check(curr_ax), 'value', 1)
                end
                show_int
            case '9'
                if get(core_check(curr_ax), 'value')
                    set(core_check(curr_ax), 'value', 0)
                else
                    set(core_check(curr_ax), 'value', 1)
                end
                show_core
            case 'a'
                pk_last
            case 'b'
                if get(cbfix_check2(curr_rad), 'value')
                    set(cbfix_check2(curr_rad), 'value', 0)
                else
                    set(cbfix_check2(curr_rad), 'value', 1)
                end
                narrow_cb
            case 'd'
                if get(distfix_check(curr_rad), 'value')
                    set(distfix_check(curr_rad), 'value', 0)
                else
                    set(distfix_check(curr_rad), 'value', 1)
                end
            case 'e'
                reset_xz
            case 'f'
                flatten
            case 'g'
                if get(grid_check(curr_ax), 'value')
                    set(grid_check(curr_ax), 'value', 0)
                else
                    set(grid_check(curr_ax), 'value', 1)
                end
                toggle_grid
            case 'l'
                load_data
            case 'm'
                pk_match
            case 'n'
                pk_next
            case 'x'
                choose_pk2
            case 'v'
                if get(nearest_check, 'value')
                    set(nearest_check, 'value', 0)
                else
                    set(nearest_check, 'value', 1)
                end
            case 'w'
                if (get(cmap_list, 'value') == 1)
                    set(cmap_list, 'value', 2)
                else
                    set(cmap_list, 'value', 1)
                end
                change_cmap
            case 'y'
                if get(zfix_check(curr_ax), 'value')
                    set(zfix_check(curr_ax), 'value', 0)
                else
                    set(zfix_check(curr_ax), 'value', 1)
                end
            case 'z'
                choose_pk1
            case 'downarrow'
                if (curr_ax == 1)
                    curr_ax = 2;
                end
                switch disp_type{curr_ax}
                    case 'amp.'                        
                        tmp1        = elev_max(curr_gui) - elev_min(curr_gui);
                        tmp2        = elev_min(curr_gui) - (0.25 * tmp1);
                        if (flat_done(curr_rad) && data_done(curr_rad))
                            tmp3    = [elev_min(curr_gui) elev_max(curr_gui)];
                        end
                        if (tmp2 < elev_min_ref)
                            elev_min(curr_gui) = elev_min_ref;
                        else
                            elev_min(curr_gui) = tmp2;
                        end
                        elev_max(curr_gui)    = elev_min(curr_gui) + tmp1;
                        if (elev_max(curr_gui) > elev_max_ref)
                            elev_max(curr_gui) = elev_max_ref;
                        end
                        if (flat_done(curr_rad) && data_done(curr_rad))
                            depth_min       = depth_min - (elev_max(curr_gui) - tmp3(2));
                            depth_max       = depth_max - (elev_min(curr_gui)  - tmp3(1));
                            if (depth_min < depth_min_ref)
                                depth_min   = depth_min_ref;
                            end
                            if (depth_max > depth_max_ref)
                                depth_max   = depth_max_ref;
                            end
                        end
                        set(z_min_edit(curr_ax), 'string', sprintf('%4.0f', elev_min(curr_gui)))
                        set(z_max_edit(curr_ax), 'string', sprintf('%4.0f', elev_max(curr_gui)))
                        if (elev_min(curr_gui) < get(z_min_slide(curr_ax), 'min'))
                            set(z_min_slide(curr_ax), 'value', get(z_min_slide(curr_ax), 'min'))
                        else
                            set(z_min_slide(curr_ax), 'value', elev_min(curr_gui))
                        end
                        if (elev_max(curr_gui) > get(z_max_slide(curr_ax), 'max'))
                            set(z_max_slide(curr_ax), 'value', get(z_max_slide(curr_ax), 'max'))
                        else
                            set(z_max_slide(curr_ax), 'value', elev_max(curr_gui))
                        end
                    case 'flat'
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
                        elev_min(curr_gui) = elev_min(curr_gui) - (depth_max - tmp3(2));
                        elev_max(curr_gui) = elev_max(curr_gui) - (depth_min - tmp3(1));
                        if (elev_min(curr_gui) < elev_min_ref)
                            elev_min(curr_gui)   = elev_min_ref;
                        end
                        if (elev_max(curr_gui) > elev_max_ref)
                            elev_max(curr_gui) = elev_max_ref;
                        end
                        set(z_min_edit(curr_ax), 'string', sprintf('%4.0f', depth_max))
                        set(z_max_edit(curr_ax), 'string', sprintf('%4.0f', depth_min))
                        if ((depth_max_ref - (depth_max - depth_min_ref)) < get(z_min_slide(curr_ax), 'min'))
                            set(z_min_slide(curr_ax), 'value', get(z_min_slide(curr_ax), 'min'))
                        elseif ((depth_max_ref - (depth_max - depth_min_ref)) > get(z_min_slide(curr_ax), 'max'))
                            set(z_min_slide(curr_ax), 'value', get(z_min_slide(curr_ax), 'max'))
                        else
                            set(z_min_slide(curr_ax), 'value', (depth_max_ref - (depth_max - depth_min_ref)))
                        end
                        if ((depth_max_ref - (depth_min - depth_min_ref)) < get(z_max_slide(curr_ax), 'min'))
                            set(z_max_slide(curr_ax), 'value', get(z_max_slide(curr_ax), 'min'))
                        elseif ((depth_max_ref - (depth_min - depth_min_ref)) > get(z_max_slide(curr_ax), 'max'))
                            set(z_max_slide(curr_ax), 'value', get(z_max_slide(curr_ax), 'max'))
                        else
                            set(z_max_slide(curr_ax), 'value', (depth_max_ref - (depth_min - depth_min_ref)))
                        end
                end
                update_z_range
                [curr_rad, curr_rad_alt] ...
                            = deal(curr_rad_alt, curr_rad);
                switch curr_ax
                    case 2
                        curr_ax = 3;
                    case 3
                        curr_ax = 2;
                end
                update_z_range
                [curr_rad, curr_rad_alt] ...
                            = deal(curr_rad_alt, curr_rad);
                switch curr_ax
                    case 2
                        curr_ax = 3;
                    case 3
                        curr_ax = 2;
                end
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
                set(dist_min_edit(curr_rad), 'string', sprintf('%3.1f', dist_min(curr_rad)))
                set(dist_max_edit(curr_rad), 'string', sprintf('%3.1f', dist_max(curr_rad)))
                if (dist_min(curr_rad) < get(dist_min_slide(curr_rad), 'min'))
                    set(dist_min_slide(curr_rad), 'value', get(dist_min_slide(curr_rad), 'min'))
                else
                    set(dist_min_slide(curr_rad), 'value', dist_min(curr_rad))
                end
                if (dist_max(curr_rad) > get(dist_max_slide(curr_rad), 'max'))
                    set(dist_max_slide(curr_rad), 'value', get(dist_max_slide(curr_rad), 'max'))
                else
                    set(dist_max_slide(curr_rad), 'value', dist_max(curr_rad))
                end
                update_dist_range
            case 'rightarrow'
                tmp1        = dist_max(curr_rad) - dist_min(curr_rad);
                tmp2        = dist_max(curr_rad) + (0.25 * tmp1);
                if (tmp2 > dist_max_ref(curr_rad));
                    dist_max(curr_rad) = dist_max_ref(curr_rad);
                else
                    dist_max(curr_rad) = tmp2;
                end
                dist_min(curr_rad) = dist_max(curr_rad) - tmp1;
                if (dist_min(curr_rad) < dist_min_ref(curr_rad))
                    dist_min(curr_rad) = dist_min_ref(curr_rad);
                end
                set(dist_min_edit(curr_rad), 'string', sprintf('%3.1f', dist_min(curr_rad)))
                set(dist_max_edit(curr_rad), 'string', sprintf('%3.1f', dist_max(curr_rad)))
                if (dist_min(curr_rad) < get(dist_min_slide(curr_rad), 'min'))
                    set(dist_min_slide(curr_rad), 'value', get(dist_min_slide(curr_rad), 'min'))
                else
                    set(dist_min_slide(curr_rad), 'value', dist_min(curr_rad))
                end
                if (dist_max(curr_rad) > get(dist_max_slide(curr_rad), 'max'))
                    set(dist_max_slide(curr_rad), 'value', get(dist_max_slide(curr_rad), 'max'))
                else
                    set(dist_max_slide(curr_rad), 'value', dist_max(curr_rad))
                end
                update_dist_range
            case 'uparrow'
                if (curr_ax == 1)
                    curr_ax = 2;
                end
                switch disp_type{curr_ax}
                    case 'amp.'
                        tmp1        = elev_max(curr_gui) - elev_min(curr_gui);
                        tmp2        = elev_max(curr_gui) + (0.25 * tmp1);
                        if (flat_done(curr_rad) && data_done(curr_rad))
                            tmp3    = [elev_min(curr_gui) elev_max(curr_gui)];
                        end
                        if (tmp2 > elev_max_ref)
                            elev_max(curr_gui) = elev_max_ref;
                        else
                            elev_max(curr_gui) = tmp2;
                        end
                        elev_min(curr_gui) = elev_max(curr_gui) - tmp1;
                        if (elev_min(curr_gui) < elev_min_ref)
                            elev_min(curr_gui) = elev_min_ref;
                        end
                        if (flat_done(curr_rad) && data_done(curr_rad))
                            depth_min       = depth_min - (elev_max(curr_gui) - tmp3(2));
                            depth_max       = depth_max - (elev_min(curr_gui) - tmp3(1));
                            if (depth_min < depth_min_ref)
                                depth_min   = depth_min_ref;
                            end
                            if (depth_max > depth_max_ref)
                                depth_max   = depth_max_ref;
                            end
                        end
                        set(z_min_edit(curr_ax), 'string', sprintf('%4.0f', elev_min(curr_gui)))
                        set(z_max_edit(curr_ax), 'string', sprintf('%4.0f', elev_max(curr_gui)))
                        if (elev_min(curr_gui) < get(z_min_slide(curr_ax), 'min'))
                            set(z_min_slide(curr_ax), 'value', get(z_min_slide(curr_ax), 'min'))
                        else
                            set(z_min_slide(curr_ax), 'value', elev_min(curr_gui))
                        end
                        if (elev_max(curr_gui) > get(z_max_slide(curr_ax), 'max'))
                            set(z_max_slide(curr_ax), 'value', get(z_max_slide(curr_ax), 'max'))
                        else
                            set(z_max_slide(curr_ax), 'value', elev_max(curr_gui))
                        end
                    case 'flat'
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
                        elev_min(curr_gui) = elev_min(curr_gui) - (depth_max - tmp3(2));
                        elev_max(curr_gui) = elev_max(curr_gui) - (depth_min - tmp3(1));
                        if (elev_min(curr_gui) < elev_min_ref)
                            elev_min(curr_gui) = elev_min_ref;
                        end
                        if (elev_max(curr_gui) > elev_max_ref)
                            elev_max(curr_gui) = elev_max_ref;
                        end
                        set(z_min_edit(curr_ax), 'string', sprintf('%4.0f', depth_max))
                        set(z_max_edit(curr_ax), 'string', sprintf('%4.0f', depth_min))
                        if ((depth_max_ref - (depth_max - depth_min_ref)) < get(z_min_slide(curr_ax), 'min'))
                            set(z_min_slide(curr_ax), 'value', get(z_min_slide(curr_ax), 'min'))
                        elseif ((depth_max_ref - (depth_max - depth_min_ref)) > get(y_min_slide, 'max'))
                            set(z_min_slide(curr_ax), 'value', get(z_min_slide(curr_ax), 'max'))
                        else
                            set(z_min_slide(curr_ax), 'value', (depth_max_ref - (depth_max - depth_min_ref)))
                        end
                        if ((depth_max_ref - (depth_min - depth_min_ref)) < get(z_max_slide(curr_ax), 'min'))
                            set(z_max_slide(curr_ax), 'value', get(z_max_slide(curr_ax), 'min'))
                        elseif ((depth_max_ref - (depth_min - depth_min_ref)) > get(z_max_slide(curr_ax), 'max'))
                            set(z_max_slide(curr_ax), 'value', get(z_max_slide(curr_ax), 'max'))
                        else
                            set(z_max_slide(curr_ax), 'value', (depth_max_ref - (depth_min - depth_min_ref)))
                        end
                end
                update_z_range
                [curr_rad, curr_rad_alt] ...
                            = deal(curr_rad_alt, curr_rad);
                switch curr_ax
                    case 2
                        curr_ax = 3;
                    case 3
                        curr_ax = 2;
                end
                update_z_range
                [curr_rad, curr_rad_alt] ...
                            = deal(curr_rad_alt, curr_rad);
                switch curr_ax
                    case 2
                        curr_ax = 3;
                    case 3
                        curr_ax = 2;
                end
            case 'space'
                switch disp_type{curr_ax}
                    case 'amp.'
                        if (flat_done(curr_rad) && data_done(curr_rad))
                            set(disp_group(curr_rad), 'selectedobject', disp_check(curr_rad, 2))
                            disp_type{curr_ax} = 'flat';
                            plot_flat
                        end
                    case 'flat'
                        set(disp_group(curr_rad), 'selectedobject', disp_check(curr_rad, 1))
                        disp_type{curr_ax} = 'amp.';
                        plot_db
                end
        end
    end

%% Test something

    function misctest(source, eventdata)
        
    end

%%
end