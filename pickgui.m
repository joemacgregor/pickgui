function pickgui
% PICKGUI Interactive radar-layer picker.
%
%   PICKGUI loads a GUI for tracing layers in CReSIS radar data, either the
%   original data files or those that have been pre-processed by
%   RADBLOCKPROC. This GUI includes semi-automatic tracing of layers and
%   layer prediction using the horizontal phase gradient and/or ARESP, if
%   the layer slopes for those methods have been added to the blocks by
%   PHASEINTERP or ARESP, respectively. Layers from separate data blocks
%   can be merged later using MERGEGUI and matched between transects using
%   FENCEGUI.
%
%   Refer to manual for operation (pickgui_man.pdf).
%
%   PICKGUI requires that the function SMOOTH_LOWESS be available within
%   the user's path. If original CReSIS data blocks are to be loaded, then
%   it also requires that either the function LL2PS be available within the
%   user's path or the Mapping Toolbox be licensed and available. If the
%   Parallel Computing Toolbox is licensed and available, then several
%   calculations related to data flattening will be parallelized.
%
% Joe MacGregor (UTIG), Mark Fahnestock (UAF-GI)
% Last updated: 10/01/14

if ~exist('smooth_lowess', 'file')
    error('pickgui:smoothlowess', 'Function SMOOTH_LOWESS is not available within this user''s path.')
end

%% Intialize variables

[pk, pk_ref]                = deal(struct);
pk.layer                    = struct;
[pk.num_layer, pk.num_man, pk.ind_trim_start, pk.ind_y_man, pk.num_keep_phase, pk.num_keep_aresp, pk.num_phase, pk.num_aresp, pk.ind_keep_phase, pk.ind_keep_aresp, pk.ind_x_start_phase, pk.ind_x_start_aresp, pk.ind_y_start_phase, pk.ind_y_start_aresp, pk.ind_y_phase_max, pk.ind_y_aresp_max] ...
                            = deal(0);
pk.predict_or_pk            = 'predict';

% two-way traveltime defaults
[twtt_min_ref, twtt_max_ref]= deal(0, 1);
[twtt_min, twtt_max]        = deal(twtt_min_ref, twtt_max_ref);

% distance defaults
[dist_min_ref, dist_max_ref]= deal(0, 1);
[dist_min, dist_max]        = deal(dist_min_ref, dist_max_ref);

% dB default
[db_min_ref, db_max_ref]    = deal(-130, 0);
[db_min, db_max]            = deal(-100, -20);

% phase difference default
[phase_diff_min_ref, phase_diff_max_ref] ...
                            = deal(-pi, pi);
[phase_diff_min, phase_diff_max] ...
                            = deal(phase_diff_min_ref, phase_diff_max_ref);

% ARESP default
[aresp_min_ref, aresp_max_ref] ...
                            = deal(-20, 20);
[aresp_min, aresp_max]      = deal(aresp_min_ref, aresp_max_ref);

% default values for several parameters
speed_vacuum                = 299792458;
permitt_ice                 = 3.15;
speed_ice                   = speed_vacuum / sqrt(permitt_ice);
decim                       = 5; % decimate radargram for display
length_chunk                = 10; % chunk length in km
int_track                   = 10; % number of indices (vertical) to separate phase- and ARESP-tracked layers
decim_flat                  = 5;
pk.num_win                  = 1; % +/- number of vertical indices in window within which to search for peak/trough in flattened layers
pk.length_smooth            = 1; % km length over which layers will be smoothed
pk.freq                     = 195e6; % radar center frequency
pk.twtt_match               = 0.01e-6; % traveltime range about which to search for matching layers
ord_poly                    = 3; % order of polynomial fit

[num_win_ref, length_smooth_ref, freq_ref, twtt_match_ref] ...
                            = deal(pk.num_win, pk.length_smooth, pk.freq, pk.twtt_match);
disp_type                   = 'twtt';
cmaps                       = {'bone' 'jet'}';
ref_start_or_end            = 'start';

if license('checkout', 'distrib_computing_toolbox')
    pool_check              = gcp('nocreate');
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
[aresp_avail, aresp_done, bed_avail, depth_avail, do_surfbed, flat_done, keep_phase_done, keep_aresp_done, load_done, load_flat, match_done, phase_avail, phase_done, pk_done, ref_done, smooth_done, surf_avail, trim_done] ...
                            = deal(false);
[amp_depth, amp_flat_mean, amp_mean, block, button, curr_chunk, curr_layer, dist_chunk, ii, ind_bed, ind_bed_flat, ind_decim, ind_decim_flat, ind_surf, ind_surf_flat, ind_decim_flat_old, ind_x_pk, ind_y_aresp, ind_y_curr, ind_y_flat, ind_y_mat, ind_y_phase, ind_y_pk, jj, kk, num_chunk, ...
 num_decim, num_decim_flat, num_sample_trim, p_aresp, p_arespdepth, p_bed, p_beddepth, p_bedflat, p_data, p_man, p_mandepth, p_phase, p_phasedepth, p_pk, p_pkdepth, p_pkflat, p_ref, p_refdepth, p_pksmooth, p_pksmoothdepth, p_pksmoothflat, p_startphase, p_startphasedepth, p_startaresp, ...
 p_startarespdepth, p_surf, p_surfflat, pkfig, rad_sample, tmp1, tmp2, tmp3, tmp4, tmp5] ...
                            = deal(0);
[ind_y_flat_mean, ind_y_flat_smooth] ...
                            = deal([]);
amp_flat                    = NaN;
[file_data, file_pk, file_ref, path_data, path_pk] ...
                            = deal('');

%% Draw the GUI

set(0, 'DefaultFigureWindowStyle', 'docked')
if ispc % windows switch
    pkgui                   = figure('toolbar', 'figure', 'name', 'PICKGUI', 'position', [1920 940 1 1], 'menubar', 'none', 'keypressfcn', @keypress, 'windowscrollwheelfcn', @wheel_zoom, 'windowbuttondownfcn', @mouse_click);
    ax_radar                = subplot('position', [0.065 0.06 1.42 0.81]);
    size_font               = 14;
    width_slide             = 0.01;
else
    pkgui                   = figure('toolbar', 'figure', 'name', 'PICKGUI', 'position', [1864 1100 1 1], 'menubar', 'none', 'keypressfcn', @keypress, 'windowscrollwheelfcn', @wheel_zoom, 'windowbuttondownfcn', @mouse_click);
    ax_radar                = subplot('position', [0.065 0.06 0.86 0.81]);
    size_font               = 18;
    width_slide             = 0.02;
end

hold on
colormap(bone)
caxis([db_min db_max])
axis ij tight
set(gca, 'fontsize', size_font, 'layer', 'top')
xlabel('Distance (km)')
ylabel('Traveltime ({\mu}s)')
colorbar('fontsize', size_font)
box on
% pan/zoom callbacks
h_pan                       = pan;
set(h_pan, 'actionpostcallback', @panzoom)
h_zoom                      = zoom;
set(h_zoom, 'actionpostcallback', @panzoom)

% sliders
twtt_min_slide              = uicontrol(pkgui, 'style', 'slider', 'units', 'normalized', 'position', [0.005 0.50 width_slide 0.32], 'callback', @slide_twtt_min, 'min', 0, 'max', 1, 'value', twtt_max_ref, 'sliderstep', [0.01 0.1]);
twtt_max_slide              = uicontrol(pkgui, 'style', 'slider', 'units', 'normalized', 'position', [0.005 0.07 width_slide 0.32], 'callback', @slide_twtt_max, 'min', 0, 'max', 1, 'value', twtt_min_ref, 'sliderstep', [0.01 0.1]);
cb_min_slide                = uicontrol(pkgui, 'style', 'slider', 'units', 'normalized', 'position', [0.96 0.07 width_slide 0.32], 'callback', @slide_db_min, 'min', -150, 'max', 0, 'value', db_min_ref, 'sliderstep', [0.01 0.1]);
cb_max_slide                = uicontrol(pkgui, 'style', 'slider', 'units', 'normalized', 'position', [0.96 0.50 width_slide 0.32], 'callback', @slide_db_max, 'min', -150, 'max', 0, 'value', db_max_ref, 'sliderstep', [0.01 0.1]);
dist_min_slide              = uicontrol(pkgui, 'style', 'slider', 'units', 'normalized', 'position', [0.12 0.005 0.27 0.02], 'callback', @slide_dist_min, 'min', 0, 'max', 1, 'value', dist_min_ref, 'sliderstep', [0.01 0.1]);
dist_max_slide              = uicontrol(pkgui, 'style', 'slider', 'units', 'normalized', 'position', [0.64 0.005 0.27 0.02], 'callback', @slide_dist_max, 'min', 0, 'max', 1, 'value', dist_max_ref, 'sliderstep', [0.01 0.1]);

% slider values
twtt_min_edit               = annotation('textbox', [0.005 0.82 0.04 0.03], 'string', num2str(twtt_min_ref), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
twtt_max_edit               = annotation('textbox', [0.005 0.39 0.04 0.03], 'string', num2str(twtt_max_ref), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
cb_min_edit                 = annotation('textbox', [0.965 0.39 0.04 0.03], 'string', num2str(db_min_ref), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
cb_max_edit                 = annotation('textbox', [0.9665 0.82 0.04 0.03], 'string', num2str(db_max_ref), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
dist_min_edit               = annotation('textbox', [0.08 0.005 0.04 0.03], 'string', num2str(dist_min_ref), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
dist_max_edit               = annotation('textbox', [0.59 0.005 0.04 0.03], 'string', num2str(dist_max_ref), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');

% push buttons
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Load data', 'units', 'normalized', 'position', [0.005 0.965 0.06 0.03], 'callback', @load_data, 'fontsize', size_font, 'foregroundcolor', 'b')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Trim y', 'units', 'normalized', 'position', [0.07 0.965 0.05 0.03], 'callback', @trim_y, 'fontsize', size_font, 'foregroundcolor', 'b')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Load picks', 'units', 'normalized', 'position', [0.30 0.965 0.06 0.03], 'callback', @load_pk, 'fontsize', size_font, 'foregroundcolor', 'b')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Load ref.', 'units', 'normalized', 'position', [0.30 0.925 0.06 0.03], 'callback', @load_ref, 'fontsize', size_font, 'foregroundcolor', 'b')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Pop figure', 'units', 'normalized', 'position', [0.30 0.885 0.06 0.03], 'callback', @pop_fig, 'fontsize', size_font, 'foregroundcolor', 'g')
phase_push                  = uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Phase', 'units', 'normalized', 'position', [0.37 0.965 0.04 0.03], 'callback', @track_phase, 'fontsize', size_font, 'foregroundcolor', 'm', 'visible', 'off');
keep_phase_push             = uicontrol(pkgui, 'style', 'pushbutton', 'string', 'keep', 'units', 'normalized', 'position', [0.42 0.965 0.04 0.03], 'callback', @pk_keep_phase, 'fontsize', size_font, 'foregroundcolor', 'm', 'visible', 'off');
aresp_push                  = uicontrol(pkgui, 'style', 'pushbutton', 'string', 'ARESP', 'units', 'normalized', 'position', [0.37 0.925 0.04 0.03], 'callback', @track_aresp, 'fontsize', size_font, 'foregroundcolor', 'm', 'visible', 'off');
keep_aresp_push             = uicontrol(pkgui, 'style', 'pushbutton', 'string', 'keep', 'units', 'normalized', 'position', [0.42 0.925 0.04 0.03], 'callback', @pk_keep_aresp, 'fontsize', size_font, 'foregroundcolor', 'm', 'visible', 'off');
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Manual', 'units', 'normalized', 'position', [0.37 0.885 0.04 0.03], 'callback', @track_man, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Flatten', 'units', 'normalized', 'position', [0.42 0.885 0.04 0.03], 'callback', @flatten, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Manually', 'units', 'normalized', 'position', [0.50 0.965 0.05 0.03], 'callback', @pk_man, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Semi-automatically', 'units', 'normalized', 'position', [0.47 0.925 0.08 0.03], 'callback', @pk_auto, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Delete', 'units', 'normalized', 'position', [0.77 0.885 0.05 0.03], 'callback', @pk_del, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Adjust', 'units', 'normalized', 'position', [0.77 0.925 0.05 0.03], 'callback', @pk_adj, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Merge', 'units', 'normalized', 'position', [0.82 0.925 0.04 0.03], 'callback', @pk_merge, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Focus', 'units', 'normalized', 'position', [0.72 0.925 0.05 0.03], 'callback', @pk_focus, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Surf./bed', 'units', 'normalized', 'position', [0.72 0.885 0.035 0.03], 'callback', @pk_surfbed, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Next', 'units', 'normalized', 'position', [0.68 0.925 0.035 0.03], 'callback', @pk_next, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Last', 'units', 'normalized', 'position', [0.64 0.925 0.035 0.03], 'callback', @pk_last, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Test', 'units', 'normalized', 'position', [0.82 0.885 0.04 0.03], 'callback', @misctest, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Match', 'units', 'normalized', 'position', [0.925 0.925 0.04 0.03], 'callback', @pk_match, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Save', 'units', 'normalized', 'position', [0.965 0.925 0.03 0.03], 'callback', @pk_save, 'fontsize', size_font, 'foregroundcolor', 'g')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Reset x/y', 'units', 'normalized', 'position', [0.945 0.885 0.05 0.03], 'callback', @reset_xy, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Reset', 'units', 'normalized', 'position', [0.005 0.46 0.03 0.03], 'callback', @reset_twtt_min, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Reset', 'units', 'normalized', 'position', [0.005 0.03 0.03 0.03], 'callback', @reset_twtt_max, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Reset', 'units', 'normalized', 'position', [0.40 0.005 0.03 0.03], 'callback', @reset_dist_min, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Reset', 'units', 'normalized', 'position', [0.92 0.005 0.03 0.03], 'callback', @reset_dist_max, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Reset', 'units', 'normalized', 'position', [0.965 0.03 0.03 0.03], 'callback', @reset_db_min, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Reset', 'units', 'normalized', 'position', [0.965 0.46 0.03 0.03], 'callback', @reset_db_max, 'fontsize', size_font, 'foregroundcolor', 'r')

% fixed text annotations
a(1)                        = annotation('textbox', [0.13 0.965 0.04 0.03], 'string', 'N_{decimate}', 'fontsize', size_font, 'color', 'b', 'edgecolor', 'none');
a(2)                        = annotation('textbox', [0.21 0.965 0.03 0.03], 'string', 'Chunk', 'fontsize', size_font, 'color', 'b', 'edgecolor', 'none');
a(3)                        = annotation('textbox', [0.21 0.925 0.03 0.03], 'string', 'L_{chunk}', 'fontsize', size_font, 'color', 'b', 'edgecolor', 'none');
a(4)                        = annotation('textbox', [0.21 0.885 0.03 0.03], 'string', 'f_{center}', 'fontsize', size_font, 'color', 'b', 'edgecolor', 'none');
a(5)                        = annotation('textbox', [0.965 0.42 0.03 0.03], 'string', 'dB_{min}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(6)                        = annotation('textbox', [0.965 0.85 0.03 0.03], 'string', 'dB_{max}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(7)                        = annotation('textbox', [0.56 0.965 0.06 0.03], 'string', 'N_{flat mean}', 'fontsize', size_font, 'color', 'm', 'edgecolor', 'none');
a(8)                        = annotation('textbox', [0.56 0.925 0.03 0.03], 'string', 'N_{win}', 'fontsize', size_font, 'color', 'm', 'edgecolor', 'none');
a(9)                        = annotation('textbox', [0.56 0.885 0.03 0.03], 'string', 'L_{smooth}', 'fontsize', size_font, 'color', 'm', 'edgecolor', 'none');
a(10)                       = annotation('textbox', [0.8575 0.885 0.03 0.03], 'string', 'Grid', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(11)                       = annotation('textbox', [0.005 0.85 0.03 0.03], 'string', 't_{min}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(12)                       = annotation('textbox', [0.005 0.42 0.03 0.03], 'string', 't_{max}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(13)                       = annotation('textbox', [0.04 0.005 0.03 0.03], 'string', 'dist_{min}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(14)                       = annotation('textbox', [0.55 0.005 0.03 0.03], 'string', 'dist_{max}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(15)                       = annotation('textbox', [0.025 0.85 0.03 0.03], 'string', 'fix', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(16)                       = annotation('textbox', [0.005 0.005 0.03 0.03], 'string', 'fix', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(17)                       = annotation('textbox', [0.95 0.005 0.03 0.03], 'string', 'fix 1', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(18)                       = annotation('textbox', [0.98 0.005 0.03 0.03], 'string', '2', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(19)                       = annotation('textbox', [0.86 0.925 0.05 0.03], 'string', '\Delta t', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(20)                       = annotation('textbox', [0.47 0.965 0.03 0.03], 'string', 'Pick', 'fontsize', size_font, 'color', 'm', 'edgecolor', 'none');
a(21)                       = annotation('textbox', [0.47 0.88 0.10 0.03], 'string', 'Smoothed layers', 'fontsize', size_font, 'color', 'm', 'edgecolor', 'none');
a(22)                       = annotation('textbox', [0.64 0.88 0.04 0.03], 'string', 'Layer', 'fontsize', size_font, 'color', 'm', 'edgecolor', 'none');
if ~ispc
    set(a, 'fontweight', 'bold')
end

% variable text annotations
file_box                    = annotation('textbox', [0.005 0.925 0.20 0.03], 'string', '', 'color', 'k', 'fontsize', size_font, 'backgroundcolor', 'w', 'edgecolor', 'k', 'interpreter', 'none');
status_box                  = annotation('textbox', [0.64 0.965 0.35 0.03], 'string', '', 'color', 'k', 'fontsize', size_font, 'backgroundcolor', 'w', 'edgecolor', 'k', 'interpreter', 'none');
cbl                         = annotation('textbox', [0.93 0.03 0.03 0.03], 'string', '', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
if ~ispc
    set(cbl, 'fontweight', 'bold')
end

% value boxes
decim_edit                  = uicontrol(pkgui, 'style', 'edit', 'string', num2str(decim), 'units', 'normalized', 'position', [0.175 0.965 0.03 0.03], 'fontsize', size_font, 'foregroundcolor', 'k', 'backgroundcolor', 'w', 'callback', @adj_decim);
length_chunk_edit           = uicontrol(pkgui, 'style', 'edit', 'string', num2str(length_chunk), 'units', 'normalized', 'position', [0.25 0.925 0.04 0.03], 'fontsize', size_font, 'foregroundcolor', 'k', 'backgroundcolor', 'w', 'callback', @adj_length_chunk);
freq_edit                   = uicontrol(pkgui, 'style', 'edit', 'string', num2str(1e-6 * pk.freq), 'units', 'normalized', 'position', [0.25 0.885 0.04 0.03], 'fontsize', size_font, 'foregroundcolor', 'k', 'backgroundcolor', 'w', 'callback', @adj_freq);
decim_flat_edit             = uicontrol(pkgui, 'style', 'edit', 'string', num2str(decim_flat), 'units', 'normalized', 'position', [0.595 0.965 0.04 0.03], 'fontsize', size_font, 'foregroundcolor', 'k', 'backgroundcolor', 'w', 'callback', @adj_decim_flat);
num_win_edit                = uicontrol(pkgui, 'style', 'edit', 'string', num2str(pk.num_win), 'units', 'normalized', 'position', [0.595 0.925 0.04 0.03], 'fontsize', size_font, 'foregroundcolor', 'k', 'backgroundcolor', 'w', 'callback', @adj_num_win);
length_smooth_edit          = uicontrol(pkgui, 'style', 'edit', 'string', num2str(pk.length_smooth), 'units', 'normalized', 'position', [0.595 0.885 0.04 0.03], 'fontsize', size_font, 'foregroundcolor', 'k', 'backgroundcolor', 'w', 'callback', @adj_length_smooth);
twtt_match_edit             = uicontrol(pkgui, 'style', 'edit', 'string', num2str(1e6 * pk.twtt_match), 'units', 'normalized', 'position', [0.88 0.925 0.04 0.03], 'fontsize', size_font, 'foregroundcolor', 'k', 'backgroundcolor', 'w', 'callback', @adj_twtt_match);

% menus
chunk_list                  = uicontrol(pkgui, 'style', 'popupmenu', 'string', 'N/A', 'value', 1, 'units', 'normalized', 'position', [0.25 0.955 0.045 0.04], 'fontsize', size_font, 'foregroundcolor', 'k', 'callback', @plot_chunk);
cmap_list                   = uicontrol(pkgui, 'style', 'popupmenu', 'string', cmaps, 'value', 1, 'units', 'normalized', 'position', [0.89 0.865 0.05 0.05], 'callback', @change_cmap, 'fontsize', size_font);
layer_list                  = uicontrol(pkgui, 'style', 'popupmenu', 'string', 'N/A', 'value', 1, 'units', 'normalized', 'position', [0.665 0.875 0.05 0.04], 'fontsize', size_font, 'foregroundcolor', 'k', 'callback', @pk_select);

% check boxes
twttfix_check               = uicontrol(pkgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.04 0.85 0.01 0.02], 'fontsize', size_font, 'value', 0, 'backgroundcolor', get(pkgui, 'color'));
distfix_check               = uicontrol(pkgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.02 0.005 0.01 0.02], 'fontsize', size_font, 'value', 0, 'backgroundcolor', get(pkgui, 'color'));
cbfix_check1                = uicontrol(pkgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.97 0.005 0.01 0.02], 'fontsize', size_font, 'value', 0, 'backgroundcolor', get(pkgui, 'color'));
cbfix_check2                = uicontrol(pkgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.985 0.005 0.01 0.02], 'callback', @narrow_cb, 'fontsize', size_font, 'value', 1, 'backgroundcolor', get(pkgui, 'color'));
ref_check                   = uicontrol(pkgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.36 0.92 0.01 0.02], 'callback', @show_ref, 'fontsize', size_font, 'value', 0, 'backgroundcolor', get(pkgui, 'color'));
man_check                   = uicontrol(pkgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.41 0.88 0.01 0.02], 'callback', @show_man, 'fontsize', size_font, 'value', 0, 'backgroundcolor', get(pkgui, 'color'));
phase_check                 = uicontrol(pkgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.41 0.96 0.01 0.02], 'callback', @show_phase, 'fontsize', size_font, 'value', 0, 'backgroundcolor', get(pkgui, 'color'), 'visible', 'off');
aresp_check                 = uicontrol(pkgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.41 0.92 0.01 0.02], 'callback', @show_aresp, 'fontsize', size_font, 'value', 0, 'backgroundcolor', get(pkgui, 'color'), 'visible', 'off');
pk_check                    = uicontrol(pkgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.55 0.92 0.01 0.02], 'callback', @show_pk, 'fontsize', size_font, 'value', 0, 'backgroundcolor', get(pkgui, 'color'));
smooth_check                = uicontrol(pkgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.55 0.88 0.01 0.02], 'callback', @show_smooth, 'fontsize', size_font, 'value', 0, 'backgroundcolor', get(pkgui, 'color'));
grid_check                  = uicontrol(pkgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.88 0.885 0.01 0.02], 'callback', @toggle_grid, 'fontsize', size_font, 'value', 0, 'backgroundcolor', get(pkgui, 'color'));
surfbed_check               = uicontrol(pkgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.76 0.88 0.01 0.02], 'callback', @show_surfbed, 'fontsize', size_font, 'value', 1, 'backgroundcolor', get(pkgui, 'color'));

% display buttons
disp_group                  = uibuttongroup('position', [0.005 0.885 0.20 0.03], 'selectionchangefcn', @disp_radio);
uicontrol(pkgui, 'style', 'text', 'parent', disp_group, 'units', 'normalized', 'position', [0 0.6 0.9 0.3], 'fontsize', size_font)
disp_check(1)               = uicontrol(pkgui, 'style', 'radio', 'string', 'twtt', 'units', 'normalized', 'position', [0.01 0.1 0.2 0.8], 'parent', disp_group, 'fontsize', size_font, 'handlevisibility', 'off');
disp_check(2)               = uicontrol(pkgui, 'style', 'radio', 'string', '~depth', 'units', 'normalized', 'position', [0.15 0.1 0.25 0.8], 'parent', disp_group, 'fontsize', size_font, 'handlevisibility', 'off', 'visible', 'off');
disp_check(3)               = uicontrol(pkgui, 'style', 'radio', 'string', 'phase', 'units', 'normalized', 'position', [0.4 0.1 0.2 0.8], 'parent', disp_group, 'fontsize', size_font, 'handlevisibility', 'off', 'visible', 'off');
disp_check(4)               = uicontrol(pkgui, 'style', 'radio', 'string', 'ARESP', 'units', 'normalized', 'position', [0.61 0.1 0.2 0.8], 'parent', disp_group, 'fontsize', size_font, 'handlevisibility', 'off', 'visible', 'off');
disp_check(5)               = uicontrol(pkgui, 'style', 'radio', 'string', 'flat', 'units', 'normalized', 'position', [0.85 0.1 0.15 0.8], 'parent', disp_group, 'fontsize', size_font, 'handlevisibility', 'off', 'visible', 'off');
set(disp_group, 'selectedobject', disp_check(1))

%% Clear plots

    function clear_plots(source, eventdata)
        if (any(p_aresp) & any(ishandle(p_aresp))) %#ok<AND2>
            delete(p_aresp(logical(p_aresp) & ishandle(p_aresp)))
        end
        if (any(p_man(:)) && any(ishandle(p_man(:))))
            delete(p_man(logical(p_man(:)) & ishandle(p_man(:))))
        end
        if (any(p_mandepth(:)) && any(ishandle(p_mandepth(:))))
            delete(p_mandepth(logical(p_mandepth(:)) & ishandle(p_mandepth(:))))
        end
        if (any(p_phase) & any(ishandle(p_phase))) %#ok<AND2>
            delete(p_phase(logical(p_phase) & ishandle(p_phase)))
        end
        if (any(p_pk) && any(ishandle(p_pk)))
            delete(p_pk(logical(p_pk) & ishandle(p_pk)))
        end
        if (any(p_pkdepth) && any(ishandle(p_pkdepth)))
            delete(p_pkdepth(logical(p_pkdepth) & ishandle(p_pkdepth)))
        end
        if (any(p_pkflat) && any(ishandle(p_pkflat)))
            delete(p_pkflat(logical(p_pkflat) & ishandle(p_pkflat)))
        end
        if (any(p_pksmooth) && any(ishandle(p_pksmooth)))
            delete(p_pksmooth(logical(p_pksmooth) & ishandle(p_pksmooth)))
        end
        if (any(p_pksmoothdepth) && any(ishandle(p_pksmoothdepth)))
            delete(p_pksmoothdepth(logical(p_pksmoothdepth) & ishandle(p_pksmoothdepth)))
        end
        if (any(p_pksmoothflat) && any(ishandle(p_pksmoothflat)))
            delete(p_pksmoothflat(logical(p_pksmoothflat) & ishandle(p_pksmoothflat)))
        end
        if (logical(p_startphase) && ishandle(p_startphase))
            delete(p_startphase)
        end
        if (logical(p_startphasedepth) && ishandle(p_startphasedepth))
            delete(p_startphasedepth)
        end
        if (logical(p_startaresp) && ishandle(p_startaresp))
            delete(p_startaresp)
        end
        if (logical(p_startarespdepth) && ishandle(p_startarespdepth))
            delete(p_startarespdepth)
        end
        pause(0.1)
        set([aresp_check phase_check man_check pk_check smooth_check surfbed_check], 'value', 0)
        set(disp_group, 'selectedobject', disp_check(1))
        set(disp_check(2:5), 'visible', 'off')
        set(layer_list, 'string', 'N/A', 'value', 1)
        set(pkgui, 'keypressfcn', @keypress, 'windowbuttondownfcn', @mouse_click)
    end

%% Clear data and picks

    function clear_data(source, eventdata)
        [pk, pk_ref]        = deal(struct);
        pk.layer            = struct;
        [pk.num_layer, pk.ind_trim_start, pk.num_man, pk.ind_y_man, pk.num_keep_phase, pk.num_keep_aresp, pk.num_aresp, pk.num_phase, pk.ind_keep_phase, pk.ind_keep_aresp, ...
         pk.ind_x_start_phase, pk.ind_x_start_aresp, pk.ind_y_start_phase, pk.ind_y_start_aresp, pk.ind_y_phase_max, pk.ind_y_aresp_max] ...
                            = deal(0);
        pk.predict_or_pk    = 'predict';
        [aresp_avail, aresp_done, bed_avail, depth_avail, do_surfbed, pk_done, flat_done, keep_phase_done, keep_aresp_done, load_done, load_flat, match_done, phase_avail, phase_done, ref_done, smooth_done, surf_avail, trim_done] ...
                            = deal(false);
        [amp_depth, amp_flat_mean, amp_mean, button, curr_chunk, curr_layer, dist_chunk, ii, ind_bed, ind_bed_flat, ind_decim, ind_surf, ind_surf_flat, ind_decim_flat_old, ind_x_pk, ind_y_aresp, ind_y_curr, ind_y_flat, ind_y_mat, ind_y_phase, ind_y_pk, jj, num_chunk, num_decim, ...
         num_decim_flat, num_sample_trim, p_aresp, p_bed, p_bedflat, p_data, p_man, p_mandepth, p_phase, p_pk, p_pkdepth, p_pkflat, p_pksmooth, p_pksmoothdepth, p_pksmoothflat, p_ref, p_refdepth, p_startaresp, p_startarespdepth, p_startphase, p_startphasedepth, p_surf, p_surfflat, pkfig, ...
         rad_sample, tmp1, tmp2, tmp3, tmp4, tmp5] ...
                            = deal(0);
        [ind_y_flat_mean, ind_y_flat_smooth] ...
                            = deal([]);
        amp_flat            = NaN;
        [file_ref, file_pk, ref_start_or_end] ...
                            = deal('');
        [pk.num_win, pk.length_smooth, pk.freq, pk.twtt_match] ...
                            = deal(num_win_ref, length_smooth_ref, freq_ref, twtt_match_ref);
        set(decim_flat_edit, 'string', num2str(decim_flat))
        set(num_win_edit, 'string', num2str(pk.num_win))
        set(length_smooth_edit, 'string', num2str(pk.length_smooth))
        set(freq_edit, 'string', num2str(1e-6 * pk.freq))
        set(twtt_match_edit, 'string', num2str(1e6 * pk.twtt_match))
    end

%% Load radar data

    function load_data(source, eventdata) %#ok<*INUSD>
        
        tmp1                = file_data;
        
        % see if user just wants to move on to next data block
        if ~isempty(file_data)
            try %#ok<TRYNC>
                tmp2        = str2double(file_data((end - 5):(end - 4))) + 1;
                if (tmp2 < 10)
                    tmp2    = ['0' num2str(tmp2)];
                else
                    tmp2    = num2str(tmp2);
                end
                if exist([path_data file_data(1:(end - 6)) tmp2 '.mat'], 'file')
                    set(status_box, 'string', 'Load next block in transect? Y: yes; otherwise: no...')
                    waitforbuttonpress
                    if strcmpi(get(pkgui, 'currentcharacter'), 'Y')
                        file_data ...
                            = [file_data(1:(end - 6)) tmp2 '.mat'];
                    end
                end
            end
        end
        
        % Dialog box to choose radar data file to load
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
        
        if isempty(file_data)
            file_data       = tmp1;
            set(status_box, 'string', 'No data loaded.')
            return
        end
        
        % attempt to load data
        set(status_box, 'string', 'Loading data...')
        pause(0.1)
        tmp1                = load([path_data file_data]);
        
        try
            
            block           = tmp1.block;
            tmp1            = 0;
            set(file_box, 'string', file_data(1:(end - 4)))
            
        catch % if not block structure than possibly an original CReSIS data file
            
            if isfield(tmp1, 'Data')
                block       = struct;
                try
                    [block.amp, block.lat, block.lon, block.num_trace, block.twtt, block.elev_air, block.dt, block.num_sample, block.file_in, block.ind_overlap] ...
                            = deal(double(tmp1.Data), tmp1.Latitude, tmp1.Longitude, length(tmp1.Latitude), tmp1.Time, tmp1.Elevation, (tmp1.Time(2) - tmp1.Time(1)), length(tmp1.Time), file_data(1:(end - 4)), NaN(1, 2));
                catch
                    set(status_box, 'string', 'Selected file does not contain expected variables.')
                    return
                end
                try
                    block.twtt_surf ...
                            = tmp1.Surface;
                catch
                    block.twtt_surf ...
                            = NaN(1, block.num_trace);
                    set(status_box, 'string', 'Selected file does not contain the Surface variable. Setting bed pick to NaN.')
                end
                try
                    block.twtt_bed ...
                            = tmp1.Bottom;
                catch
                    block.twtt_bed ...
                            = NaN(1, block.num_trace);
                    set(status_box, 'string', 'Selected file does not contain the Bottom variable. Setting bed pick to NaN.')
                end
                block.amp   = single(block.amp);
                if isrow(block.twtt)
                    block.twtt ...
                            = block.twtt';
                end
                block.dist  = 1e-3 .* cumsum([0 distance([block.lat(1:(end - 1))' block.lon(1:(end - 1))'], [block.lat(2:end)' block.lon(2:end)'], wgs84Ellipsoid)']);
                block.dist_lin ...
                            = interp1([1 block.num_trace], block.dist([1 end]), 1:block.num_trace);
                tmp1        = 0;
                set(file_box, 'string', file_data(6:(end - 4)))
                
            else
                set(status_box, 'string', 'Selected file does not contain expected variables.')
                return
            end
        end
        
        if load_done
            set(ref_check, 'value', 0)
            if (logical(p_bed) && ishandle(p_bed))
                delete(p_bed)
            end
            if (logical(p_beddepth) && ishandle(p_beddepth))
                delete(p_beddepth)
            end
            if (logical(p_bedflat) && ishandle(p_bedflat))
                delete(p_bedflat)
            end
            if (logical(p_data) && ishandle(p_data))
                delete(p_data)
            end
            if (any(p_ref) && any(ishandle(p_ref)))
                delete(p_ref(logical(p_ref) & ishandle(p_ref)))
            end
            if (any(p_refdepth) && any(ishandle(p_refdepth)))
                delete(p_refdepth(logical(p_refdepth) & ishandle(p_refdepth)))
            end
            if (logical(p_surf) && ishandle(p_surf))
                delete(p_surf)
            end
            if (logical(p_surfflat) && ishandle(p_surfflat))
                delete(p_surfflat)
            end
            set([disp_check(2:5) aresp_check aresp_push keep_aresp_push keep_phase_push phase_check phase_push], 'visible', 'off')
            clear_plots
            clear_data
            set(twttfix_check, 'value', 0)
            load_done       = false;
            set(disp_group, 'selectedobject', disp_check(1))
            disp_type       = 'twtt';
            pause(0.1)
        end
        
        if ~isfield(block, 'phase_diff_filt') % no phase to work with
            phase_avail     = false;
        else
            set(phase_push, 'visible', 'on')
            set(phase_check, 'visible', 'on')
            set(keep_phase_push, 'visible', 'on')
            set(disp_check(3), 'visible', 'on')
            [phase_diff_min_ref, phase_diff_min] ...
                            = deal(nanmin(block.phase_diff_filt(:)));
            [phase_diff_max_ref, phase_diff_max] ...
                            = deal(nanmax(block.phase_diff_filt(:)));
            phase_avail     = true;
        end
        
        if ~isfield(block, 'slope_aresp') % no aresp
            aresp_avail     = false;
        else
            set(aresp_push, 'visible', 'on')
            set(aresp_check, 'visible', 'on')
            set(keep_aresp_push, 'visible', 'on')
            set(disp_check(4), 'visible', 'on')
            [aresp_min_ref, aresp_min] ...
                            = deal(atand(nanmin(block.slope_aresp(:))));
            [aresp_max_ref, aresp_max] ...
                            = deal(atand(nanmax(block.slope_aresp(:))));
            aresp_avail     = true;
        end
        
        % decimation vector
        if (decim > 1)
            ind_decim       = (1 + ceil(decim / 2)):decim:(block.num_trace - ceil(decim / 2));
        else
            ind_decim       = 1:block.num_trace;
        end
        num_decim           = length(ind_decim);
        
        if (decim_flat > 1)
            ind_decim_flat  = ceil((decim_flat / 2) + 1):decim_flat:(block.num_trace - ceil(decim_flat / 2));
        else
            ind_decim_flat  = 1:block.num_trace;
        end
        num_decim_flat      = length(ind_decim_flat);
        
        % convert to dB
        block.amp(block.amp == 0) ...
                            = NaN;
        block.amp(isinf(block.amp)) ...
                            = NaN;
        block.amp           = 10 .* log10(abs(block.amp));
        num_sample_trim     = block.num_sample;
        if (length(block.twtt) ~= block.num_sample) % fix for data that were preprocessed with radblockproc before this same fix was implemented there
            block.twtt      = interp1(1:length(block.twtt), block.twtt, 1:block.num_sample, 'linear', 'extrap');
        end
        if isrow(block.twtt)
            block.twtt      = block.twtt';
        end
        if (decim > 1)
            amp_mean        = NaN(num_sample_trim, num_decim, 'single');
            tmp1            = floor(decim / 2);
            for ii = 1:num_decim
                amp_mean(:, ii) ...
                            = nanmean(block.amp(:, (ind_decim(ii) - tmp1):(ind_decim(ii) + tmp1)), 2);
            end
        else
            amp_mean        = single(block.amp);
        end
        
        % make chunks
        adj_length_chunk
        
        % assign traveltime and distance reference values/sliders based on data
        [twtt_min_ref, twtt_max_ref, twtt_min, twtt_max, db_min_ref, db_max_ref, db_min, db_max, dist_min_ref, dist_max_ref, dist_min, dist_max] ...
                            = deal(block.twtt(1), block.twtt(end), block.twtt(1), block.twtt(end), min(block.amp(:)), max(block.amp(:)), min(block.amp(:)), max(block.amp(:)), block.dist_lin(1), block.dist_lin(end), block.dist_lin(1), block.dist_lin(end));
        set(twtt_min_slide, 'min', (1e6 * twtt_min_ref), 'max', (1e6 * twtt_max_ref), 'value', (1e6 * twtt_max_ref))
        set(twtt_max_slide, 'min', (1e6 * twtt_min_ref), 'max', (1e6 * twtt_max_ref), 'value', (1e6 * twtt_min_ref))
        set(dist_min_slide, 'min', dist_min_ref, 'max', dist_max_ref, 'value', dist_min_ref)
        set(dist_max_slide, 'min', dist_min_ref, 'max', dist_max_ref, 'value', dist_max_ref)
        set(twtt_min_edit, 'string', sprintf('%3.1f', (1e6 * twtt_min_ref)))
        set(twtt_max_edit, 'string', sprintf('%3.1f', (1e6 * twtt_max_ref)))
        set(dist_min_edit, 'string', sprintf('%3.1f', dist_min_ref))
        set(dist_max_edit, 'string', sprintf('%3.1f', dist_max_ref))
        set(cb_min_slide, 'min', db_min_ref, 'max', db_max_ref, 'value', db_min)
        set(cb_max_slide, 'min', db_min_ref, 'max', db_max_ref, 'value', db_max)
        update_dist_range
        update_twtt_range
        
        [ind_surf, ind_bed] = deal(NaN(1, block.num_trace));
        if isfield(block, 'twtt_surf')
            if any(~isnan(block.twtt_surf) & any(~isinf(block.twtt_surf)))
                ind_surf(~isnan(block.twtt_surf) & ~isinf(block.twtt_surf)) ...
                            = interp1(block.twtt, 1:num_sample_trim, block.twtt_surf(~isnan(block.twtt_surf) & ~isinf(block.twtt_surf)), 'nearest', 'extrap');
                if any(isnan(ind_surf))
                    ind_surf(isnan(ind_surf)) ...
                            = round(interp1(find(~isnan(ind_surf)), ind_surf(~isnan(ind_surf)), find(isnan(ind_surf)), 'linear', 'extrap'));
                end
                ind_surf(ind_surf <= 0) ...
                            = 1;
                surf_avail  = true;
                p_surf      = plot(block.dist_lin(ind_decim), (1e6 .* block.twtt_surf(ind_decim)), 'm.', 'markersize', 12, 'visible', 'off');
            end
        end
        if isfield(block, 'twtt_bed')
            ind_bed         = NaN(1, block.num_trace);
            if any(~isnan(block.twtt_bed) & any(~isinf(block.twtt_bed)))
                ind_bed(~isnan(block.twtt_bed) & ~isinf(block.twtt_bed)) ...
                            = interp1(block.twtt, 1:num_sample_trim, block.twtt_bed(~isnan(block.twtt_bed) & ~isinf(block.twtt_bed)), 'nearest', 'extrap');
                ind_bed(ind_bed > num_sample_trim) ...
                            = num_sample_trim;
                bed_avail   = true;
                set(surfbed_check, 'value', 1)
                p_bed       = plot(block.dist_lin(ind_decim), (1e6 .* block.twtt_bed(ind_decim)), 'm.', 'markersize', 12, 'visible', 'off');
            end
        end
        
        if surf_avail
            if ~isempty(find(~isnan(ind_surf(ind_decim)), 1))
                amp_depth   = NaN(num_sample_trim, num_decim, 'single');
                for ii = find(~isnan(ind_surf(ind_decim)))
                    amp_depth(1:(num_sample_trim - ind_surf(ind_decim(ii)) + 1), ii) ...
                            = amp_mean(ind_surf(ind_decim(ii)):num_sample_trim, ii); % shift data up to surface
                end
                set(disp_check(2), 'visible', 'on')
                depth_avail = true;
                if bed_avail
                    tmp1    = ind_decim(~isnan(ind_bed(ind_decim)) & ~isnan(ind_surf(ind_decim)));
                    p_beddepth ...
                            = plot(block.dist_lin(tmp1), (1e6 .* (block.twtt_bed(tmp1) - block.twtt_surf(tmp1) + block.twtt(1))), 'm.', 'markersize', 12, 'visible', 'off');
                end
                set(surfbed_check, 'value', 1)
            end
        end
        
        set(layer_list, 'string', {'surface' 'bed'}, 'value', 1)
        pk_select
        
        % plot data
        set(disp_group, 'selectedobject', disp_check(1))
        load_done           = true;
        plot_twtt
        set(status_box, 'string', 'Block loaded.')
    end

%% Trim excess data (e.g., before surface reflection and after deepest bed reflection)

    function trim_y(source, eventdata)
        
        if ~load_done
            set(status_box, 'string', 'No data to trim.')
            return
        end
        if ((pk.num_keep_phase || pk.num_keep_aresp || pk.num_man) && ~pk_done)
            set(status_box, 'string', 'Can''t trim once predicted layers have been made.')
            return
        end
        if (trim_done && flat_done)
            set(status_box, 'string', 'Data that has already been flattened should not be trimmed twice.')
            return
        end
        if ~strcmp(disp_type, 'twtt')
            set(disp_group, 'selectedobject', disp_check(1))
            disp_type       = 'twtt';
            plot_twtt
        end
        
        set(status_box, 'string', 'Trimming data...')
        pause(0.1)
        
        % trim most of the data
        tmp1                = (interp1(block.twtt, 1:num_sample_trim, twtt_min, 'nearest', 'extrap'):interp1(block.twtt, 1:num_sample_trim, twtt_max, 'nearest', 'extrap'))';
        if trim_done
            tmp2            = pk.ind_trim_start + tmp1(1);
        end
        pk.ind_trim_start   = tmp1(1);
        
        if surf_avail
            ind_surf        = ind_surf - pk.ind_trim_start;
        end
        if bed_avail
            ind_bed         = ind_bed - pk.ind_trim_start;
        end
        block.amp           = block.amp(tmp1, :);
        amp_mean            = amp_mean(tmp1, :);
        if depth_avail
            tmp2            = amp_depth;
            amp_depth       = NaN(size(amp_mean), 'single');
            for ii = find(~isnan(ind_surf(ind_decim)))
                amp_depth(1:(length(tmp1) - ind_surf(ind_decim(ii)) + 1), ii) ...
                            = tmp2(1:(length(tmp1) - ind_surf(ind_decim(ii)) + 1), ii);
            end
        end
        
        if phase_avail
            block.phase_diff_filt ...
                            = block.phase_diff_filt(tmp1, :);
        end
        if aresp_avail
            block.slope_aresp ...
                            = block.slope_aresp(tmp1, :);
        end
        num_sample_trim     = length(tmp1);
        
        % if layers already exist, then trim them too
        if pk_done
            for ii = 1:pk.num_layer
                pk.layer(ii).ind_y ...
                            = pk.layer(ii).ind_y - pk.ind_trim_start;
            end
            if any(smooth_done)
                for ii = find(smooth_done)
                    pk.layer(ii).ind_y_smooth ...
                            = pk.layer(ii).ind_y_smooth - pk.ind_trim_start;
                end
            end
        end
        
        % save traveltime for the end
        block.twtt          = block.twtt(tmp1);
        [twtt_min_ref, twtt_max_ref] ...
                            = deal(block.twtt(1), block.twtt(end));
        set(twtt_min_slide, 'min', (1e6 * twtt_min_ref), 'max', (1e6 * twtt_max_ref), 'value', (1e6 * twtt_max_ref))
        set(twtt_max_slide, 'min', (1e6 * twtt_min_ref), 'max', (1e6 * twtt_max_ref), 'value', (1e6 * twtt_min_ref))
        
        if (bed_avail && depth_avail)
            if (logical(p_beddepth) && ishandle(p_beddepth))
                delete(p_beddepth)
            end
            tmp2            = ind_decim(~isnan(ind_bed(ind_decim)) & ~isnan(ind_surf(ind_decim)));
            p_beddepth      = plot(block.dist_lin(tmp2), (1e6 .* (block.twtt_bed(tmp2) - block.twtt_surf(tmp2) + block.twtt(1))), 'm.', 'markersize', 12, 'visible', 'off');
        end        
        
        if trim_done
            pk.ind_trim_start ...
                            = tmp2;
        else
            trim_done       = true;
        end
        
        update_twtt_range
        set(status_box, 'string', ['Data trimmed off before ' num2str((1e6 * block.twtt(1)), '%2.1f') ' us and after ' num2str((1e6 * block.twtt(end)), '%2.1f') ' us.'])
    end

%% Load existing layer picks for this block

    function load_pk(source, eventdata)
        
        if ~load_done
            set(status_box, 'string', 'Load data before picks.')
            return
        end
        if trim_done
            set(status_box, 'string', 'Cannot load picks if data have already been trimmed. Reload data.')
            return
        end
        if ~strcmp(disp_type, 'twtt')
            set(disp_group, 'selectedobject', disp_check(1))
            disp_type       = 'twtt';
            plot_twtt
        end
        
        tmp1                = file_data(1:11);
        tmp2                = path_data(end - 1);
        if (isnan(str2double(tmp2)) || ~isreal(str2double(tmp2))) % check for a/b/c/etc in file_pk_short
            if ispc
                tmp3        = '..\..\..\pk\';
            else
                tmp3        = '../../../pk/';
            end
        else
            if ispc
                tmp3        = '..\..\pk\';
            else
                tmp3        = '../../pk/';
            end
        end
        
        % look for picks file in the expected location
        if ispc
            if (~isempty(path_data) && exist([path_data tmp3 file_data(1:11)], 'dir'))
                path_pk     = [path_data(1:strfind(path_data, '\block')) 'pk\' file_data(1:11) '\'];
            elseif (~isempty(path_data) && exist([path_data tmp3 file_data(1:11) '\' tmp2(end)], 'dir'))
                path_pk     = [path_data(1:strfind(path_data, '\block')) 'pk\' file_data(1:11) '\' tmp2(end) '\'];
            end
        else
            if (~isempty(path_data) && exist([path_data tmp3 file_data(1:11)], 'dir'))
                path_pk     = [path_data(1:strfind(path_data, '/block')) 'pk/' file_data(1:11) '/'];
            elseif (~isempty(path_data) && exist([path_data tmp3 file_data(1:11) '/' tmp2(end)], 'dir'))
                path_pk     = [path_data(1:strfind(path_data, '/block')) 'pk/' file_data(1:11) '/' tmp2(end) '/'];
            end
        end
        
        if (~isempty(path_pk) && exist([path_pk file_data(1:(end - 4)) '_pk.mat'], 'file'))
            file_pk         = [file_data(1:(end - 4)) '_pk.mat'];
        else % Dialog box to choose picks file to load
            [tmp1, tmp2]    = deal(file_pk, path_pk);
            if ~isempty(path_pk)
                [file_pk, path_pk] = uigetfile('*.mat', 'Load picks:', path_pk);
            elseif ~isempty(path_data)
                [file_pk, path_pk] = uigetfile('*.mat', 'Load picks:', path_data);
            else
                [file_pk, path_pk] = uigetfile('*.mat', 'Load picks:');
            end
            if ~ischar(file_pk)
                [file_pk, path_pk] = deal('', tmp2);
            end
        end
        
        if isempty(file_pk)
            file_pk         = tmp1;
            set(status_box, 'string', 'No picks loaded.')
            return
        end
        
        % check against data file
        try %#ok<TRYNC>
            if ~strcmp(file_data(1:(end - 4)), file_pk(1:(end - 7)))
                set(status_box, 'string', ['Selected picks file (' file_pk(1:(end - 4)) ') might not match data file. Continue loading? Y: yes; otherwise: no...'])
                waitforbuttonpress
                if ~strcmpi(get(pkgui, 'currentcharacter'), 'Y')
                    set(status_box, 'string', 'Loading of picks file cancelled.')
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
            set(status_box, 'string', 'Selected file does not contained a pk structure.')
            return
        end
        
        if ~pk.num_layer
            set(status_box, 'string', 'Picks file has no picks.')
            return
        end
        
        clear_plots
        
        [aresp_done, flat_done, keep_aresp_done, keep_phase_done, load_flat, match_done, phase_done, pk_done] ...
                            = deal(false);
        
        if isfield(pk, 'poly_flat')
            load_flat       = true;
            ord_poly        = size(pk.poly_flat, 1) - 1;
        else % sort layers by index
            tmp1            = NaN(pk.num_layer, block.num_trace);
            tmp2            = NaN(1, pk.num_layer);
            for ii = 1:pk.num_layer
                pk.layer(ii).ind_y ...
                            = interp1(block.twtt, 1:num_sample_trim, pk.layer(ii).twtt, 'nearest', 'extrap');
                pk.layer(ii).ind_y_smooth ...
                            = interp1(block.twtt, 1:num_sample_trim, pk.layer(ii).twtt_smooth, 'nearest', 'extrap');
                tmp1(ii, :) = pk.layer(ii).twtt;
                tmp2(ii)    = nanmean(pk.layer(ii).ind_y);
            end
            [~, tmp3]       = sort(tmp2);
            pk.layer        = pk.layer(tmp3);
            if (length(find(all(~isnan(tmp1')))) > 3)
                pk.predict_or_pk ...
                            = 'pk';
            else
                pk.predict_or_pk ...
                            = 'predict';
            end
        end
        
        if ~isfield(pk, 'predict_or_pk')
            pk.predict_or_pk= 'predict';
        end
        
        set(freq_edit, 'string', num2str(1e-6 * pk.freq))
        set(num_win_edit, 'string', num2str(pk.num_win))
        set(length_smooth_edit, 'string', num2str(pk.length_smooth))
        set(twtt_match_edit, 'string', num2str(1e6 * pk.twtt_match))
        
        set(layer_list, 'string', [num2cell(1:pk.num_layer) 'surface' 'bed'], 'value', 1)
        
        % trim data as it was before
        [twtt_min, twtt_min_ref] ...
                            = deal(pk.twtt_min_ref);
        if surf_avail
            if (twtt_min_ref > min(block.twtt_surf(~isinf(block.twtt_surf))))
                [twtt_min, twtt_min_ref] ...
                            = deal(min(block.twtt_surf(~isinf(block.twtt_surf))));
            end
        end
        [twtt_max, twtt_max_ref] ...
                            = deal(pk.twtt_max_ref);
        if bed_avail
            if (twtt_max_ref < max(block.twtt_bed(~isinf(block.twtt_bed))))
                [twtt_max, twtt_max_ref] ...
                            = deal(max(block.twtt_bed(~isinf(block.twtt_bed))));
            end
        end
        set(twtt_min_edit, 'string', sprintf('%3.1f', (1e6 * twtt_min_ref)))
        set(twtt_max_edit, 'string', sprintf('%3.1f', (1e6 * twtt_max_ref)))
        pk_done             = true;
        smooth_done         = true(1, pk.num_layer);
        pause(0.1)
        set(status_box, 'string', 'Continuing pick loading...')
        
        % change old variable names
        if isfield(pk, 'num_layer_keep')
            [pk.ind_keep_phase, pk.num_keep_phase] ...
                            = deal(pk.ind_layer_keep, pk.num_layer_keep);
            pk              = rmfield(pk, {'ind_layer_keep' 'num_layer_keep'});
            if isfield(pk, 'ind_x_start')
                pk.ind_x_start_phase ...
                            = pk.ind_x_start;
                pk          = rmfield(pk, 'ind_x_start');
            else
                pk.ind_x_start_phase ...
                            = 1;
            end
            if isfield(pk, 'ind_y_start')
                pk.ind_y_start_phase ...
                            = pk.ind_y_start;
                pk          = rmfield(pk, 'ind_y_start');
            end
        end
        
        % load and plot phase tracking
        if (phase_avail && pk.num_keep_phase)
            prop_phase
            p_phase         = zeros(1, pk.num_phase);
            for ii = 1:pk.num_phase
                p_phase(ii) = plot(block.dist_lin(ind_decim), (1e6 .* block.twtt(round(ind_y_phase(ii, ind_decim)))), 'b', 'linewidth', 1, 'visible', 'off'); % plot phase-tracked layers
            end
            p_startphase    = plot((ones(1, pk.num_phase) .* block.dist_lin(pk.ind_x_start_phase)), (1e6 .* block.twtt(pk.ind_y_start_phase)), 'm.', 'markersize', 12, 'visible', 'off');
            if depth_avail
                p_phasedepth= zeros(1, pk.num_phase);
                for ii = 1:pk.num_phase
                        tmp1= ind_decim(~isnan(ind_y_phase(ii, ind_decim)) & ~isnan(ind_surf(ind_decim)));
                        tmp2= round(ind_y_phase(ii, tmp1) - ind_surf(tmp1) + 1);
                        tmp2(tmp2 <= 0) ...
                            = 1;
                    p_phasedepth(ii) ...
                            = plot(block.dist_lin(tmp1), (1e6 .* block.twtt(tmp2)), 'b', 'linewidth', 1, 'visible', 'off');
                end
                p_startphasedepth ...
                            = plot((ones(1, pk.num_phase) .* block.dist_lin(pk.ind_x_start_phase)), (1e6 .* (block.twtt(pk.ind_y_start_phase) - (1e6 .* (block.twtt_surf(pk.ind_x_start_phase) - block.twtt(1))))), 'm.', 'markersize', 12, 'visible', 'off');
            end
            phase_done      = true;
            if pk.num_keep_phase
                for ii = 1:pk.num_keep_phase
                    set([p_phase(pk.ind_keep_phase(ii)) p_phasedepth(pk.ind_keep_phase(ii))], 'linewidth', 2, 'color', 'w') % increase linewidth and change color
                end
                keep_phase_done ...
                            = true;
            end
        else
            pk.num_keep_phase ...
                            = 0;
        end
        
        % do same for ARESP, change old variable names
        if isfield(pk, 'num_layer_keep_aresp')
            [pk.ind_keep_aresp, pk.num_keep_aresp] ...
                            = deal(pk.ind_layer_keep_aresp, pk.num_layer_keep_aresp);
            pk              = rmfield(pk, {'ind_layer_keep_aresp' 'num_layer_keep_aresp'});
        end
        
        if (aresp_avail && isfield(pk, 'num_keep_aresp'))
            if pk.num_keep_aresp
                prop_aresp
                p_aresp     = zeros(1, pk.num_aresp);
                for ii = 1:pk.num_aresp
                    p_aresp(ii) ...
                            = plot(block.dist_lin(ind_decim), (1e6 .* block.twtt(round(ind_y_aresp(ii, ind_decim)))), 'c', 'linewidth', 1, 'visible', 'off'); % plot ARESP-tracked layers
                end
                p_startaresp= plot((ones(1, pk.num_aresp) .* block.dist_lin(pk.ind_x_start_aresp)), (1e6 .* block.twtt(pk.ind_y_start_aresp)), 'm.', 'markersize', 12, 'visible', 'off'); % plot the starting y indices
                if depth_avail
                    p_arespdepth ...
                            = zeros(1, pk.num_aresp);
                    for ii = 1:pk.num_aresp
                        tmp1= ind_decim(~isnan(ind_y_aresp(ii, ind_decim)) & ~isnan(ind_surf(ind_decim)));
                        tmp2= round(ind_y_aresp(ii, tmp1) - ind_surf(tmp1) + 1);
                        tmp2(tmp2 <= 0) ...
                            = 1;
                        p_arespdepth(ii) ...
                            = plot(block.dist_lin(tmp1), (1e6 .* block.twtt(tmp2)), 'c', 'linewidth', 1, 'visible', 'off');
                    end
                    p_startarespdepth ...
                            = plot((ones(1, pk.num_aresp) .* block.dist_lin(pk.ind_x_start_aresp)), (1e6 .* (block.twtt(pk.ind_y_start_aresp) - (1e6 .* (block.twtt_surf(pk.ind_x_start_aresp) - block.twtt(1))))), 'm.', 'markersize', 12, 'visible', 'off');
                end
                aresp_done  = true;
                if pk.num_keep_aresp
                    for ii = 1:pk.num_keep_aresp
                        set([p_aresp(pk.ind_keep_aresp(ii)) p_arespdepth(pk.ind_keep_aresp(ii))], 'linewidth', 2, 'color', 'w')
                    end
                    keep_aresp_done ...
                            = true;
                end
            end
        else
            pk.num_keep_aresp ...
                            = 0;
        end
        
        % do same for manual layers
        if pk.num_man
            if any((pk.ind_y_man(:) < 1) | (pk.ind_y_man(:) > num_sample_trim)) % trimming not done when manual picks were done (no longer possible)
                set(status_box, 'string', 'Cannot display manually picked layers due to trimming issue...')
                [pk.num_man, pk.ind_y_man] ...
                            = deal(0);
                pause(0.1)
            else
                p_man       = zeros(pk.num_man, 2);
                if depth_avail
                    p_mandepth ...
                            = zeros(pk.num_man, 2);
                end
                for ii = 1:pk.num_man
                    p_man(ii, 1) ...
                            = plot(block.dist_lin(ind_decim), (1e6 .* block.twtt(round(pk.ind_y_man(ii, ind_decim)))), 'linewidth', 2, 'color', [0.85 0.85 0.85], 'visible', 'off'); % max-picking spline
                    p_man(ii, 2) ...
                            = plot(block.dist_lin(ind_decim), (1e6 .* block.twtt(round(pk.ind_y_man(ii, ind_decim)))), 'linewidth', 2, 'color', [0.85 0.85 0.85], 'visible', 'off');
                    if depth_avail
                        tmp1= ind_decim(~isnan(pk.ind_y_man(ii, ind_decim)) & ~isnan(ind_surf(ind_decim)));
                        p_mandepth(ii, 1) ...
                            = plot(block.dist_lin(tmp1), (1e6 .* block.twtt(round(pk.ind_y_man(ii, tmp1) - ind_surf(tmp1) + 1))), 'linewidth', 2, 'color', [0.85 0.85 0.85], 'visible', 'off');
                        p_mandepth(ii, 2) ...
                            = plot(block.dist_lin(tmp1), (1e6 .* block.twtt(round(pk.ind_y_man(ii, tmp1) - ind_surf(tmp1) + 1))), 'linewidth', 2, 'color', [0.85 0.85 0.85], 'visible', 'off');
                    end
                end
            end
        end
        
        trim_y
        
        % remove old fields
        if isfield(pk.layer, 'ind_y_flat_mean')
            pk.layer        = rmfield(pk.layer, 'ind_y_flat_mean');
        end
        if isfield(pk.layer, 'ind_y_flat_smooth')
            pk.layer        = rmfield(pk.layer, 'ind_y_flat_smooth');
        end
        if isfield(pk.layer, 'type')
            pk.layer        = rmfield(pk.layer, 'type');
        end
        if isfield(pk, 'ind_x_mean')
            pk              = rmfield(pk, 'ind_x_mean');
        end
        if isfield(pk, 'keep_or_flat')
            pk              = rmfield(pk, 'keep_or_flat');
        end
        if isfield(pk, 'num_ind_mean')
            pk              = rmfield(pk, 'num_ind_mean');
        end
        
        % plot layers in twtt/~depth
        [p_pk, p_pkdepth, p_pkflat, p_pksmooth, p_pksmoothdepth, p_pksmoothflat] ...
                            = deal(zeros(1, pk.num_layer));
        [ind_y_flat_mean, ind_y_flat_smooth] ...
                            = deal(NaN(pk.num_layer, num_decim_flat));
        tmp3                = [];
        for ii = 1:pk.num_layer
            if ~isempty(find(~isnan(pk.layer(ii).ind_y(ind_decim)), 1))
                try
                    p_pk(ii)= plot(block.dist_lin(ind_decim(~isnan(pk.layer(ii).ind_y(ind_decim)))), (1e6 .* block.twtt(round(pk.layer(ii).ind_y(ind_decim(~isnan(pk.layer(ii).ind_y(ind_decim))))))), 'r.', 'markersize', 12, 'visible', 'off');
                    if depth_avail
                        tmp1= ind_decim(~isnan(pk.layer(ii).ind_y(ind_decim)) & ~isnan(ind_surf(ind_decim)));
                        tmp2= pk.layer(ii).ind_y(tmp1) - ind_surf(tmp1) + 1;
                        tmp2((tmp2 < 1) | (tmp2 > num_sample_trim)) ...
                            = NaN;
                        p_pkdepth(ii) ...
                            = plot(block.dist_lin(tmp1(~isnan(tmp2))), (1e6 .* block.twtt(round(tmp2(~isnan(tmp2))))), 'r.', 'markersize', 12, 'visible', 'off');
                    end
                catch
                    tmp3    = [tmp3 ii]; %#ok<AGROW>
                    continue
                end
            else
                tmp3        = [tmp3 ii]; %#ok<AGROW>
                continue
            end
            if ~isempty(find(~isnan(pk.layer(ii).ind_y_smooth(ind_decim)), 1))
                try
                    tmp1    = ind_decim(~isnan(pk.layer(ii).ind_y_smooth(ind_decim)));
                    p_pksmooth(ii) ...
                            = plot(block.dist_lin(tmp1), (1e6 .* block.twtt(round(pk.layer(ii).ind_y_smooth(tmp1)))), 'g.', 'markersize', 12, 'visible', 'off');
                    if depth_avail
                        tmp1= ind_decim(~isnan(pk.layer(ii).ind_y_smooth(ind_decim)) & ~isnan(ind_surf(ind_decim)));
                        tmp2= pk.layer(ii).ind_y_smooth(tmp1) - ind_surf(tmp1) + 1;
                        tmp2((tmp2 < 1) | (tmp2 > num_sample_trim)) ...
                            = NaN;
                        p_pksmoothdepth(ii) ...
                            = plot(block.dist_lin(tmp1(~isnan(tmp2))), (1e6 .* block.twtt(round(tmp2(~isnan(tmp2))))), 'g.', 'markersize', 12, 'visible', 'off');
                    end
                catch
                    tmp3    = [tmp3 ii]; %#ok<AGROW>
                end
            else
                tmp3        = [tmp3 ii]; %#ok<AGROW>
            end
        end
        
        % get rid of empty layers
        if ~isempty(tmp3)
            pk.num_layer
            tmp2            = setdiff(1:pk.num_layer, tmp3);
            for ii = tmp3
                if (logical(p_pk(ii)) && ishandle(p_pk(ii)))
                    delete(p_pk(ii))
                end
                if (logical(p_pkdepth(ii)) && ishandle(p_pkdepth(ii)))
                    delete(p_pkdepth(ii))
                end
                if (logical(p_pksmooth(ii)) && ishandle(p_pksmooth(ii)))
                    delete(p_pksmooth(ii))
                end
                if (logical(p_pksmoothdepth(ii)) && ishandle(p_pksmoothdepth(ii)))
                    delete(p_pksmoothdepth(ii))
                end
            end
            if (length(tmp2) > 1)
                curr_layer  = interp1(1:length(tmp2), tmp2, curr_layer, 'nearest', 'extrap');
            else
                curr_layer  = 1;
            end
            [pk.num_layer, pk.layer, smooth_done, p_pk, p_pkdepth, p_pkflat, p_pksmooth, p_pksmoothdepth, p_pksmoothflat] ...
                            = deal(length(tmp2), pk.layer(tmp2), smooth_done(tmp2), p_pk(tmp2), p_pkdepth(tmp2), p_pkflat(tmp2), p_pksmooth(tmp2), p_pksmoothdepth(tmp2), p_pksmoothflat(tmp2));            
            set(layer_list, 'string', [num2cell(1:pk.num_layer) 'surface' 'bed'], 'value', curr_layer)
        end
        
        if load_flat % reflatten if possible
            
            pause(0.1)
            set(status_box, 'string', 'Flattening data...')
            
            % load polynomial flattening
            ind_y_mat       = single((1:size(block.amp, 1))');
            ind_y_mat       = ind_y_mat(:, ones(1, block.num_trace)); % matrix of y indices
            switch ord_poly
                case 2
                    ind_y_flat ...
                            = ((ind_y_mat .^ 2) .* pk.poly_flat(ones(num_sample_trim, 1), :)) + (ind_y_mat .* (pk.poly_flat((2 .* ones(num_sample_trim, 1)), :))) + pk.poly_flat((3 .* ones(num_sample_trim, 1)), :);
                case 3
                    ind_y_flat ...
                            = ((ind_y_mat .^ 3) .* pk.poly_flat(ones(num_sample_trim, 1), :)) + ((ind_y_mat .^ 2) .* pk.poly_flat((2 .* ones(num_sample_trim, 1)), :)) + ...
                              (ind_y_mat .* (pk.poly_flat((3 .* ones(num_sample_trim, 1)), :))) + pk.poly_flat((4 .* ones(num_sample_trim, 1)), :);
            end
            ind_y_mat       = 0;
            ind_y_flat(ind_y_flat < 1) ...
                            = 1;
            ind_y_flat(ind_y_flat > num_sample_trim) ...
                            = num_sample_trim;
            amp_flat        = NaN(size(block.amp, 1), block.num_trace, 'single');
            tmp2            = find(~sum(isnan(ind_y_flat)));
            
            if parallel_check
                tmp1        = block.amp(:, tmp2);
                tmp3        = ind_y_flat(:, tmp2);
                tmp4        = amp_flat(:, tmp2);
                parfor ii = 1:length(tmp2)
                    tmp4(:, ii) ...
                            = interp1(tmp1(:, ii), tmp3(:, ii));
                end
                amp_flat(:, tmp2) ...
                            = tmp4;
                tmp1        = 0;
            else
                for ii = 1:length(tmp2)
                    amp_flat(:, tmp2(ii)) ...
                            = interp1(block.amp(:, tmp2(ii)), ind_y_flat(:, tmp2(ii)));
                end
            end
            flat_done       = true;
            
            % flatten layers
            warning('off', 'MATLAB:interp1:NaNinY')
            tmp1            = NaN(pk.num_layer, num_decim_flat);
            for ii = 1:pk.num_layer
                tmp1(ii, :) = pk.layer(ii).ind_y(ind_decim_flat);
            end
            for ii = 1:num_decim_flat
                if sum(isnan(ind_y_flat(:, ind_decim_flat(ii))))
                    continue
                end
                [~, tmp3]   = unique(ind_y_flat(:, ind_decim_flat(ii)));
                if any(~isnan(tmp1(:, ii)))
                    ind_y_flat_mean(~isnan(tmp1(:, ii)), ii) ...
                            = interp1(ind_y_flat(tmp3, ind_decim_flat(ii)), tmp3, tmp1(~isnan(tmp1(:, ii)), ii), 'nearest', 'extrap');
                end
            end
            
            tmp1            = NaN(pk.num_layer, num_decim_flat);
            for ii = 1:pk.num_layer
                tmp1(ii, :) = pk.layer(ii).ind_y_smooth(ind_decim_flat);
            end
            for ii = 1:num_decim_flat
                if sum(isnan(ind_y_flat(:, ind_decim_flat(ii))))
                    continue
                end
                [~, tmp3]   = unique(ind_y_flat(:, ind_decim_flat(ii)));
                if any(~isnan(tmp1(:, ii)))
                    ind_y_flat_smooth(~isnan(tmp1(:, ii)), ii) ...
                            = interp1(ind_y_flat(tmp3, ind_decim_flat(ii)), tmp3, tmp1(~isnan(tmp1(:, ii)), ii), 'nearest', 'extrap');
                end
            end
            warning('on', 'MATLAB:interp1:NaNinY')
            
            % flatten surface pick
            if surf_avail
                ind_surf_flat ...
                            = NaN(1, block.num_trace);
                for ii = find(~sum(isnan(ind_y_flat)))
                    [~, tmp2] ...
                            = unique(ind_y_flat(:, ii));
                    if (length(tmp2) > 1)
                        ind_surf_flat(ii) ...
                            = interp1(ind_y_flat(tmp2, ii), tmp2, ind_surf(ii), 'nearest', 'extrap');
                    end
                end
                ind_surf_flat ...
                            = round(ind_surf_flat);
                ind_surf_flat((ind_surf_flat < 1) | (ind_surf_flat > num_sample_trim)) ...
                            = NaN;
                if ~all(isnan(ind_surf_flat(ind_decim)))
                    tmp1    = ind_surf_flat(ind_decim_flat);
                    p_surfflat ...
                            = plot(block.dist_lin(ind_decim_flat(~isnan(tmp1))), (1e6 .* block.twtt(tmp1(~isnan(tmp1)))), 'm.', 'markersize', 12, 'visible', 'off');
                end
            end
            
            % flatten bed pick
            if bed_avail
                ind_bed_flat= NaN(1, block.num_trace);
                for ii = find(~sum(isnan(ind_y_flat)))
                    [~, tmp2] ...
                            = unique(ind_y_flat(:, ii));
                    if (length(tmp2) > 1)
                        ind_bed_flat(ii) ...
                            = interp1(ind_y_flat(tmp2, ii), tmp2, ind_bed(ii), 'nearest', 'extrap');
                    end
                end
                ind_bed_flat= round(ind_bed_flat);
                ind_bed_flat((ind_bed_flat < 1) | (ind_bed_flat > num_sample_trim)) ...
                            = NaN;
                if ~all(isnan(ind_bed_flat(ind_decim)))
                    tmp1    = ind_bed_flat(ind_decim_flat);
                    p_bedflat ...
                            = plot(block.dist_lin(ind_decim_flat(~isnan(tmp1))), (1e6 .* block.twtt(tmp1(~isnan(tmp1)))), 'm.', 'markersize', 12, 'visible', 'off');
                end
            end
            
            mean_flat
            pause(0.1)
            set(disp_check(5), 'visible', 'on')
            match_done      = true;
        end
        
        set([pk_check smooth_check], 'value', 1)
        if depth_avail
            set(disp_check(2), 'visible', 'on')
        end
        if phase_avail
            set(disp_check(3), 'visible', 'on')
        end
        if aresp_avail
            set(disp_check(4), 'visible', 'on')
        end
        show_phase
        show_aresp
        show_ref
        show_man
        show_pk
        show_smooth
        if load_flat
            set(status_box, 'string', ['Picks loaded from ' file_pk(1:(end - 4)) '.'])
        else
            set(status_box, 'string', ['Picks loaded from ' file_pk(1:(end - 4)) ' (no flattening).'])
        end
    end

%% Load reference layer picks from previous/next block

    function load_ref(source, eventdata)
        
        if ~load_done
            set(status_box, 'string', 'Data not loaded yet.')
            return
        end
        
        if ispc
            tmp1            = '..\..\pk\';
        else
            tmp1            = '../../pk/';
        end
        
        % check if reference layers can be quickly located
        try %#ok<TRYNC>
            if (~isempty(path_data) && exist([path_data tmp1 file_data(1:11)], 'dir'))
                if ispc
                    if isempty(path_pk)
                        path_pk = [path_data(1:strfind(path_data, '\block')) 'pk\' file_data(1:11) '\'];
                    elseif ~strcmp(path_pk(1:(end - 2)), [path_data(1:strfind(path_data, '\block')) 'pk\' file_data(1:11) '\'])
                        path_pk = [path_data(1:strfind(path_data, '\block')) 'pk\' file_data(1:11) '\'];
                    end
                else
                    if isempty(path_pk)
                        path_pk = [path_data(1:strfind(path_data, '/block')) 'pk/' file_data(1:11) '/'];
                    elseif ~strcmp(path_pk(1:(end - 2)), [path_data(1:strfind(path_data, '/block')) 'pk/' file_data(1:11) '/'])
                        path_pk = [path_data(1:strfind(path_data, '/block')) 'pk/' file_data(1:11) '/'];
                    end
                end
            end
        end
        
        if isempty(file_pk)
            file_pk         = [file_data(1:(end - 4)) '_pk.mat'];
        end
        
        ref_start_or_end       = '';
        
        % check for existing pick files before and after current block
        if (~isempty(path_pk) && ~isempty(file_pk))
            try %#ok<TRYNC>
                tmp1        = str2double(file_pk((end - 8):(end - 7))) - 1;
                if (tmp1 < 10)
                    tmp1    = ['0' num2str(tmp1)];
                else
                    tmp1    = num2str(tmp1);
                end
                tmp2        = str2double(file_pk((end - 8):(end - 7))) + 1;
                if (tmp2 < 10)
                    tmp2    = ['0' num2str(tmp2)];
                else
                    tmp2    = num2str(tmp2);
                end
                if (exist([path_pk file_pk(1:(end - 9)) tmp1 '_pk.mat'], 'file') && ~exist([path_pk file_pk(1:(end - 9)) tmp2 '_pk.mat'], 'file'))
                    file_ref= [file_pk(1:(end - 9)) tmp1 '_pk.mat'];
                    ref_start_or_end ...
                            = 'start';
                elseif (~exist([path_pk file_pk(1:(end - 9)) tmp1 '_pk.mat'], 'file') && exist([path_pk file_pk(1:(end - 9)) tmp2 '_pk.mat'], 'file'))
                    set(status_box, 'string', 'Only right-hand-side reference picks available. Load these? Y: yes; otherwise: no...')
                    waitforbuttonpress
                    if strcmpi(get(pkgui, 'currentcharacter'), 'Y')
                        file_ref ...
                            = [file_pk(1:(end - 9)) tmp2 '_pk.mat'];
                        ref_start_or_end ...
                            = 'end';
                    else
                        [file_ref, ref_start_or_end] ...
                            = deal('');
                    end
                elseif (exist([path_pk file_pk(1:(end - 9)) tmp1 '_pk.mat'], 'file') && exist([path_pk file_pk(1:(end - 9)) tmp2 '_pk.mat'], 'file'))
                    set(status_box, 'string', 'Left (L) or right (R) reference picks?')
                    waitforbuttonpress
                    if strcmpi(get(pkgui, 'currentcharacter'), 'L')
                        file_ref= [file_pk(1:(end - 9)) tmp1 '_pk.mat'];
                        ref_start_or_end ...
                            = 'start';
                    elseif strcmpi(get(pkgui, 'currentcharacter'), 'R')
                        file_ref= [file_pk(1:(end - 9)) tmp2 '_pk.mat'];
                        ref_start_or_end ...
                            = 'end';
                    end
                end
            end
        end
        
        % Dialog box to choose reference picks file to load
        if isempty(file_ref)
            tmp1            = path_pk;
            if ~isempty(path_pk)
                [file_ref, path_pk] = uigetfile('*.mat', 'Load picks:', path_pk);
            elseif ~isempty(path_data)
                [file_ref, path_pk] = uigetfile('*.mat', 'Load picks:', path_data);
            else
                [file_ref, path_pk] = uigetfile('*.mat', 'Load picks:');
            end
            if ~ischar(file_ref)
                file_ref    = '';
                path_pk     = tmp1;
            end
        end
        
        if isempty(file_ref)
            set(status_box, 'string', 'No reference picks loaded.')
            return
        end
        
        if strcmp(file_ref, file_pk)
            set(status_box, 'string', 'Reference picks file should not be same as current block''s pick file.')
            return
        end
        
        % load reference picks file
        tmp1                = load([path_pk file_ref]);
        try
            pk_ref          = tmp1.pk;
            tmp1            = 0;
        catch
            set(status_box, 'string', 'Selected file does not contain a pk structure.')
            return
        end
        
        if (any(p_ref) && any(ishandle(p_ref)))
            delete(p_ref(logical(p_ref) & ishandle(p_ref)))
        end
        if (any(p_refdepth) && any(ishandle(p_refdepth)))
            delete(p_refdepth(logical(p_refdepth) & ishandle(p_refdepth)))
        end
        
        ref_done            = false;
        set(status_box, 'string', 'Loading reference picks...')
        
        % plot reference picks
        p_ref               = zeros(1, pk_ref.num_layer);
        
        if ~isempty(ref_start_or_end)
            switch ref_start_or_end
                case 'start'
                    for ii = 1:pk_ref.num_layer
                        p_ref(ii) = plot(block.dist_lin(1:decim:block.ind_overlap(1)), (1e6 .* pk_ref.layer(ii).twtt_smooth(pk_ref.ind_overlap(2):decim:end)), 'y.', 'markersize', 12, 'visible', 'off');
                    end
                case 'end'
                    for ii = 1:pk_ref.num_layer
                        p_ref(ii) = plot(block.dist_lin(block.ind_overlap(2):decim:end), (1e6 .* pk_ref.layer(ii).twtt_smooth(1:decim:pk_ref.ind_overlap(1))), 'y.', 'markersize', 12, 'visible', 'off');
                    end
            end
        elseif ~isnan(pk_ref.ind_overlap(2))
            if (abs(pk_ref.dist(pk_ref.ind_overlap(2)) - block.dist(1)) < 0.01)
                ref_start_or_end ...
                            = 'start'; % reference picks from start of transect
                for ii = 1:pk_ref.num_layer
                    p_ref(ii) = plot(block.dist_lin(1:decim:block.ind_overlap(1)), (1e6 .* pk_ref.layer(ii).twtt_smooth(pk_ref.ind_overlap(2):decim:end)), 'y.', 'markersize', 12, 'visible', 'off');
                end
            end
        elseif ~isnan(pk_ref.ind_overlap(1))
            if (abs(pk_ref.dist(pk_ref.ind_overlap(1)) - block.dist(end)) < 0.01)
                ref_start_or_end ...
                            = 'end';
                for ii = 1:pk_ref.num_layer
                    p_ref(ii) = plot(block.dist_lin(block.ind_overlap(2):decim:end), (1e6 .* pk_ref.layer(ii).twtt_smooth(1:decim:pk_ref.ind_overlap(1))), 'y.', 'markersize', 12, 'visible', 'off');
                end
            end
        else
            set(status_box, 'string', 'This picks file does not overlap with the current block.')
            pk_ref          = struct;
            return
        end
        
        if depth_avail % plot reference picks in ~depth
            p_refdepth      = zeros(1, pk_ref.num_layer);
            switch ref_start_or_end
                case 'start'
                    for ii = 1:pk_ref.num_layer
                        p_refdepth(ii) = plot(block.dist_lin(1:decim:block.ind_overlap(1)), (1e6 .* (pk_ref.layer(ii).twtt_smooth(pk_ref.ind_overlap(2):decim:end) - pk_ref.twtt_surf(pk_ref.ind_overlap(2):decim:end) + block.twtt(1))), 'y.', 'markersize', 12, 'visible', 'off');
                    end
                case 'end'
                    for ii = 1:pk_ref.num_layer
                        p_refdepth(ii) = plot(block.dist_lin(block.ind_overlap(2):decim:end), (1e6 .* (pk_ref.layer(ii).twtt_smooth(1:decim:pk_ref.ind_overlap(1)) - pk_ref.twtt_surf(1:decim:pk_ref.ind_overlap(1)) + block.twtt(1))), 'y.', 'markersize', 12, 'visible', 'off');
                    end
            end
        end
        
        set(ref_check, 'value', 1)
        ref_done        = true;
        show_ref
        set(status_box, 'string', ['Reference picks loaded from ' file_ref(1:(end - 4)) '.'])
    end

%% Track layers using the horizontal phase gradient

    function track_phase(source, eventdata)

        if ~load_done
            set(status_box, 'string', 'Load data first.')
            return
        end
        if ~phase_avail
            set(status_box, 'string', 'Horizontal phase gradient not available.')
            return
        end
        if ~any(strcmp(disp_type, {'twtt' '~depth'}))
            set(status_box, 'string', 'ARESP must be traced in twtt or ~depth.')
            return
        end
        
        if (logical(p_startphase) && ishandle(p_startphase))
            delete(p_startphase)
        end
        if (logical(p_startphasedepth) && ishandle(p_startphasedepth))
            delete(p_startphasedepth)
        end
        if (any(p_phase) && any(ishandle(p_phase)))
            delete(p_phase(logical(p_phase) & ishandle(p_phase)))
        end
        if (any(p_phasedepth) && any(ishandle(p_phasedepth)))
            delete(p_phasedepth(logical(p_phasedepth) & ishandle(p_phasedepth)))
        end
        
        if keep_phase_done
            [pk.num_keep_phase, pk.ind_keep_phase] ...
                            = deal(0);
        end
        [phase_done, keep_phase_done] ...
                            = deal(false);
        
        set(phase_check, 'value', 0)
        set(status_box, 'string', 'Choose trace to propagate layers from...(Q: cancel)')
        
        set(pkgui, 'keypressfcn', [], 'windowbuttondownfcn', [])
        
        while true
            
            [ind_x_pk, ind_y_pk, button] ...
                            = ginput(1); % get x/y indices to start tracing from (e.g., thickest ice)
            
            if strcmpi(char(button), 'Q')
                
                set(status_box, 'string', 'Cancelled tracking phase.')
                set(pkgui, 'keypressfcn', @keypress, 'windowbuttondownfcn', @mouse_click)
                return
                
            elseif (button == 1) % left-click
                
                pause(0.1)
                set(status_box, 'string', 'Tracking phase...')
                
                pk.ind_x_start_phase ...
                            = interp1(block.dist_lin(ind_decim), ind_decim, ind_x_pk, 'nearest', 'extrap'); % convert to index
                
                switch disp_type
                    case 'twtt'
                        pk.ind_y_phase_max ...
                            = interp1(block.twtt, 1:num_sample_trim, (1e-6 * ind_y_pk), 'nearest', 'extrap'); % maximum index to pick, not much in data below here for layer tracing; used in flattening section
                    case '~depth'
                        pk.ind_y_phase_max ...
                            = interp1(block.twtt, 1:num_sample_trim, (1e-6 * (ind_y_pk + (1e6 * (block.twtt_surf(pk.ind_x_start_phase) - block.twtt(1))))), 'nearest', 'extrap'); % maximum index to pick, not much in data below here for layer tracing; used in flattening section
                end
                
                if (block.twtt_surf(pk.ind_x_start_phase) > twtt_min)
                    pk.ind_y_start_phase ...
                            = (interp1(block.twtt, 1:num_sample_trim, block.twtt_surf(pk.ind_x_start_phase), 'nearest', 'extrap') + int_track):int_track:pk.ind_y_phase_max; % y indices from which to propagate phase-tracked layers
                else
                    pk.ind_y_start_phase ...
                            = (1 + int_track):int_track:pk.ind_y_phase_max; % y indices from which to propagate first attempt at phase-tracked layers
                end
                pk.num_phase= length(pk.ind_y_start_phase); % number of y indices to test in first attempt
                
                prop_phase % propagate phase in a separate sub-function
                
                axes(ax_radar) %#ok<*LAXES>
                p_phase     = zeros(1, pk.num_phase);
                for ii = 1:pk.num_phase
                    p_phase(ii) = plot(block.dist_lin(ind_decim(~isnan(ind_y_phase(ii, ind_decim)))), (1e6 .* block.twtt(round(ind_y_phase(ii, ind_decim(~isnan(ind_y_phase(ii, ind_decim))))))), 'b', 'linewidth', 1, 'visible', 'off'); % plot the phase-tracked layers
                end
                p_startphase= plot((ones(1, pk.num_phase) .* block.dist_lin(pk.ind_x_start_phase)), (1e6 .* block.twtt(pk.ind_y_start_phase)), 'm.', 'markersize', 12, 'visible', 'off'); % plot the starting y indices
                if depth_avail
                    p_phasedepth ...
                            = zeros(1, pk.num_phase);
                    for ii = 1:pk.num_phase
                        tmp1= ind_decim(~isnan(ind_y_phase(ii, ind_decim)) & ~isnan(ind_surf(ind_decim)));
                        tmp2= round(ind_y_phase(ii, tmp1) - ind_surf(tmp1) + 1);
                        tmp2(tmp2 <= 0) ...
                            = 1;
                        tmp2(tmp2 > num_sample_trim) ...
                            = num_sample_trim;                        
                        p_phasedepth(ii) ...
                            = plot(block.dist_lin(tmp1), (1e6 .* block.twtt(tmp2)), 'b', 'linewidth', 1, 'visible', 'off'); % plot the phase-tracked layers
                    end
                    tmp1    = pk.ind_y_start_phase - ind_surf(pk.ind_x_start_phase) + 1;                    
                    tmp1    = tmp1((tmp1 > 0) & (tmp1 < num_sample_trim));                    
                    p_startphasedepth ...
                            = plot((ones(1, length(tmp1)) .* block.dist_lin(pk.ind_x_start_phase)), (1e6 .* block.twtt(tmp1)), 'm.', 'markersize', 12, 'visible', 'off');
                end
                set(phase_check, 'value', 1)
                phase_done  = true;
                show_phase
                set(status_box, 'string', ['Traced ' num2str(pk.num_phase) ' layers starting from ' num2str(block.dist(pk.ind_x_start_phase), '%3.1f') ' km.'])
                break
            end
        end
        
        set(pkgui, 'keypressfcn', @keypress, 'windowbuttondownfcn', @mouse_click)
    end

%% Propagate horizontal phase gradient

    function prop_phase(source, eventdata)
        
        rad_sample          = (pk.freq * block.dt) / (2 * pi); % radians traveled per time sampling interval
        tmp1                = block.phase_diff_filt ./ rad_sample; % convert into index space translated
        
        ind_y_phase         = NaN(pk.num_phase, block.num_trace, 'single'); % y indices of first-attempt layer tracing by phase propagation
        ind_y_phase(:, pk.ind_x_start_phase) ...
                            = pk.ind_y_start_phase'; % assign starting y index at starting x index
        
        % loop through each starter y index at pk.ind_x_start_phase and propagate from that x index using the filtered phase gradient
        for ii = 1:pk.num_phase %#ok<*FXUP>
            for jj = (pk.ind_x_start_phase - 1):-1:1 % loop heading left of pk.ind_x_start_phase
                ind_y_phase(ii, jj)     = ind_y_phase(ii, (jj + 1)) + tmp1(round(ind_y_phase(ii, (jj + 1))), (jj + 1));
                if (ind_y_phase(ii, jj) < 1)
                    ind_y_phase(ii, jj) = 1;
                elseif (ind_y_phase(ii, jj) > num_sample_trim)
                    ind_y_phase(ii, jj) = num_sample_trim;
                end
            end
            for jj = (pk.ind_x_start_phase + 1):block.num_trace % loop heading right of pk.ind_x_start_phase
                ind_y_phase(ii, jj)     = ind_y_phase(ii, (jj - 1)) - tmp1(round(ind_y_phase(ii, (jj - 1))), (jj - 1));
                if (ind_y_phase(ii, jj) < 1)
                    ind_y_phase(ii, jj) = 1;
                elseif (ind_y_phase(ii, jj) > num_sample_trim)
                    ind_y_phase(ii, jj) = num_sample_trim;
                end
            end
        end
        
        tmp1                = 0;
        
        % remove phase-tracked layers that NaN'd out
        tmp1                = find(sum(isnan(ind_y_phase), 2));
        if ~isempty(tmp1)
            ind_y_phase     = ind_y_phase(setdiff(1:pk.num_phase, tmp1), :);
            pk.ind_y_start_phase ...
                            = pk.ind_y_start_phase(setdiff(1:pk.num_phase, tmp1));
            pk.num_phase    = pk.num_phase - length(tmp1);
        end
    end

%% Pick phase-propagated layers to keep for flattening

    function pk_keep_phase(source, eventdata)
        
        if ~phase_done
            set(status_box, 'string', 'No phase-tracked layers yet.')
            return
        end
        if ~any(strcmp(disp_type, {'twtt' '~depth'}))
            set(status_box, 'string', 'ARESP layers can only be kept in twtt or ~depth.')
            return
        end
        if ~get(phase_check, 'value')
            set(phase_check, 'value', 1)
            show_phase
            keep_phase_done = false;
        end
        
        set(status_box, 'string', 'Pick best phase-traced layers for flattening...')
        pause(0.1)
        set(status_box, 'string', 'Left-click: select; return/enter: keep; D: unkeep; Q: stop...')
        
        axes(ax_radar)
        ii                  = 0;
        
        set(pkgui, 'keypressfcn', [], 'windowbuttondownfcn', [])
        
        while true
            
            [ind_x_pk, ind_y_pk, button] ...
                            = ginput(1); % point near layer to keep
            if ~ii % dummy start
                ii          = 1;
                continue
            end
            
            if ~any(strcmpi(char(button), {'Q' 'D'}))
                if (button ~= 1)
                    continue
                end
            end
            
            if strcmpi(char(button), 'Q') % quit/stop picking keepers
                
                break
                
            elseif strcmpi(char(button), 'D') % delete a kept layer
                
                if pk.num_keep_phase
                    ind_x_pk= interp1(block.dist_lin(ind_decim), ind_decim, ind_x_pk, 'nearest', 'extrap');
                    switch disp_type
                        case 'twtt'
                            ind_y_pk = interp1(block.twtt, 1:num_sample_trim, (1e-6 * ind_y_pk), 'nearest', 'extrap');
                        case '~depth'
                            ind_y_pk = interp1(block.twtt, 1:num_sample_trim, (1e-6 * (ind_y_pk + (1e6 * (block.twtt_surf(ind_x_pk) - block.twtt(1))))), 'nearest', 'extrap');
                    end
                    tmp1    = pk.ind_keep_phase(~isnan(ind_y_phase(pk.ind_keep_phase, ind_x_pk)));
                    if (length(tmp1) > 1)
                        tmp1= interp1(ind_y_phase(tmp1, ind_x_pk), tmp1, ind_y_pk, 'nearest', 'extrap');
                    end
                    set(p_phase(tmp1), 'color', 'b', 'linewidth', 1)
                    if depth_avail
                        set(p_phasedepth(tmp1), 'color', 'b', 'linewidth', 1)
                    end
                    if (pk.num_keep_phase > 1)
                        tmp2=interp1(pk.ind_keep_phase, 1:pk.num_keep_phase, tmp1, 'nearest', 'extrap');
                    else
                        tmp2= 1;
                    end
                    pk.ind_keep_phase ...
                            = pk.ind_keep_phase(setdiff(1:pk.num_keep_phase, tmp2));
                    pk.num_keep_phase ...
                            = pk.num_keep_phase - 1;
                    set(status_box, 'string', 'Layer unkept. Pick again...')
                end
                
            else
                
                ind_x_pk    = interp1(block.dist_lin(ind_decim), ind_decim, ind_x_pk, 'nearest', 'extrap'); % convert x index to non-decimated range
                tmp2        = find(~isnan(ind_y_phase(:, ind_x_pk)));
                tmp1        = ind_y_phase(tmp2, ind_x_pk);
                [tmp1, tmp3]= unique(tmp1);
                if (length(tmp1) > 1)
                    switch disp_type
                        case 'twtt'
                            ind_y_pk = interp1(tmp1, tmp2(tmp3), interp1(block.twtt, 1:num_sample_trim, (1e-6 * ind_y_pk), 'nearest', 'extrap'), 'nearest', 'extrap'); % index of nearest layer
                        case '~depth'
                            ind_y_pk = interp1(tmp1, tmp2(tmp3), interp1(block.twtt, 1:num_sample_trim, (1e-6 * (ind_y_pk + (1e6 * (block.twtt_surf(ind_x_pk) - block.twtt(1))))), 'nearest', 'extrap'), 'nearest', 'extrap');
                    end
                else
                    ind_y_pk= tmp2;
                end
                if any(ind_y_pk == pk.ind_keep_phase)
                    set(status_box, 'string', 'Chosen layer already kept. Pick again...')
                    continue
                end
                set(p_phase(ind_y_pk), 'color', 'w') % change layer color to white
                if depth_avail
                    set(p_phasedepth(ind_y_pk), 'color', 'w')
                end
                waitforbuttonpress
                
                if (double(get(pkgui, 'currentcharacter')) == 13)
                    set(p_phase(ind_y_pk), 'linewidth', 2) % increase linewidth
                    if depth_avail
                        set(p_phasedepth(ind_y_pk), 'linewidth', 2) % increase linewidth
                    end                        
                    pk.ind_keep_phase(pk.num_keep_phase + 1) = ind_y_pk; % preserve layer number
                    pk.num_keep_phase = pk.num_keep_phase + 1;
                else
                    set(p_phase(ind_y_pk), 'color', 'b') % return to previous color
                    if depth_avail
                        set(p_phasedepth(ind_y_pk), 'color', 'b')
                    end
                    set(status_box, 'string', 'Then pick again...')
                end
            end
        end
        
        if pk.num_keep_phase
            pk.ind_keep_phase ...
                            = sort(pk.ind_keep_phase); % sort kept layer numbers
            set(status_box, 'string', ['Picked ' num2str(pk.num_keep_phase) ' phase-tracked layers to keep.'])
            keep_phase_done = true;
        else
            set(status_box, 'string', 'No layers kept.')
        end
        
        set(pkgui, 'keypressfcn', @keypress, 'windowbuttondownfcn', @mouse_click)
    end

%% Track layers using ARESP

    function track_aresp(source, eventdata)

        if ~load_done
            set(status_box, 'string', 'Load data first.')
            return
        end
        if ~aresp_avail
            set(status_box, 'string', 'ARESP not available.')
            return
        end
        if ~any(strcmp(disp_type, {'twtt' '~depth'}))
            set(status_box, 'string', 'ARESP must be traced in twtt or ~depth.')
            return
        end
        
        if (logical(p_startaresp) && ishandle(p_startaresp))
            delete(p_startaresp)
        end
        if (logical(p_startarespdepth) && ishandle(p_startarespdepth))
            delete(p_startarespdepth)
        end
        if (any(p_aresp) && any(ishandle(p_aresp)))
            delete(p_aresp(logical(p_aresp) & ishandle(p_aresp)))
        end
        if (any(p_arespdepth) && any(ishandle(p_arespdepth)))
            delete(p_arespdepth(logical(p_arespdepth) & ishandle(p_arespdepth)))
        end
        
        if keep_aresp_done
            [pk.num_keep_aresp, pk.ind_keep_aresp] ...
                            = deal(0);
        end
        [aresp_done, keep_aresp_done] ...
                            = deal(false);
        
        set(aresp_check, 'value', 0)
        set(status_box, 'string', 'Choose trace to propagate layers from...(Q: cancel)')
        
        set(pkgui, 'keypressfcn', [], 'windowbuttondownfcn', [])
        axes(ax_radar)
        
        while true
            
            [ind_x_pk, ind_y_pk, button] ...
                            = ginput(1); % get x index to start tracing from (e.g., thickest ice)
            
            if strcmpi(char(button), 'Q')
                
                set(status_box, 'string', 'Cancelled tracking ARESP.')
                set(pkgui, 'keypressfcn', @keypress, 'windowbuttondownfcn', @mouse_click)
                return
                
            elseif (button == 1) % left-click
                
                pause(0.1)
                set(status_box, 'string', 'Tracking ARESP slopes...')
                
                pk.ind_x_start_aresp ...
                            = interp1(block.dist_lin(ind_decim), ind_decim, ind_x_pk, 'nearest', 'extrap'); % convert to index
                switch disp_type
                    case 'twtt'
                        pk.ind_y_aresp_max ...
                            = interp1(block.twtt, 1:num_sample_trim, (1e-6 * ind_y_pk), 'nearest', 'extrap'); % maximum index to pick, not much in data below here for layer tracing
                    case '~depth'
                        pk.ind_y_aresp_max ...
                            = interp1(block.twtt, 1:num_sample_trim, (1e-6 * (ind_y_pk + (1e6 * (block.twtt_surf(pk.ind_x_start_aresp) - block.twtt(1))))), 'nearest', 'extrap');
                end
                if (block.twtt_surf(pk.ind_x_start_aresp) > twtt_min)
                    pk.ind_y_start_aresp ...
                            = (interp1(block.twtt, 1:num_sample_trim, block.twtt_surf(pk.ind_x_start_aresp), 'nearest', 'extrap') + int_track):int_track:pk.ind_y_aresp_max; % y indices to propagate ARESP layers from
                else
                    pk.ind_y_start_aresp ...
                            = (1 + int_track):int_track:pk.ind_y_aresp_max; % y indices from which to propagate first attempt at ARESP-tracked layers
                end
                pk.ind_y_start_aresp ...
                            = pk.ind_y_start_aresp(~isnan(block.slope_aresp(pk.ind_y_start_aresp, pk.ind_x_start_aresp)));
                pk.num_aresp= length(pk.ind_y_start_aresp); % number of y indices to test in first attempt
                
                prop_aresp % propagate ARESP in a separate sub-function
                
                p_aresp     = zeros(1, pk.num_aresp);
                for ii = 1:pk.num_aresp
                    p_aresp(ii) ...
                            = plot(block.dist_lin(ind_decim(~isnan(ind_y_aresp(ii, ind_decim)))), (1e6 .* block.twtt(round(ind_y_aresp(ii, ind_decim(~isnan(ind_y_aresp(ii, ind_decim))))))), 'c', 'linewidth', 1, 'visible', 'off');
                end
                p_startaresp= plot((ones(1, pk.num_aresp) .* block.dist_lin(pk.ind_x_start_aresp)), (1e6 .* block.twtt(pk.ind_y_start_aresp)), 'r.', 'markersize', 12, 'visible', 'off'); % plot the starting y indices
                if depth_avail
                    p_arespdepth ...
                            = zeros(1, pk.num_aresp);
                    for ii = 1:pk.num_aresp
                        tmp1= ind_decim(~isnan(ind_y_aresp(ii, ind_decim)) & ~isnan(ind_surf(ind_decim)));
                        tmp2= round(ind_y_aresp(ii, tmp1) - ind_surf(tmp1) + 1);
                        tmp2(tmp2 <= 0) ...
                            = 1;
                        tmp2(tmp2 > num_sample_trim) ...
                            = num_sample_trim;
                        p_arespdepth(ii) ...
                            = plot(block.dist_lin(tmp1), (1e6 .* block.twtt(tmp2)), 'c', 'linewidth', 1, 'visible', 'off');
                    end
                    tmp1    = pk.ind_y_start_aresp - ind_surf(pk.ind_x_start_aresp) + 1;
                    tmp1    = tmp1((tmp1 > 0) & (tmp1 < num_sample_trim));
                    p_startarespdepth ...
                            = plot((ones(1, length(tmp1)) .* block.dist_lin(pk.ind_x_start_aresp)), (1e6 .* block.twtt(tmp1)), 'r.', 'markersize', 12, 'visible', 'off'); % plot the starting y indices
                end
                set(aresp_check, 'value', 1)
                aresp_done  = true;
                show_aresp
                set(status_box, 'string', ['Traced ' num2str(pk.num_aresp) ' layers starting from ' num2str(block.dist(pk.ind_x_start_aresp), '%3.1f') ' km.'])
                break
            end
        end
        
        set(pkgui, 'keypressfcn', @keypress, 'windowbuttondownfcn', @mouse_click)
    end

%% Propagate ARESP

    function prop_aresp(source, eventdata)
        
        ind_y_aresp         = NaN(pk.num_aresp, block.num_trace, 'single'); % y indices of first-attempt layer tracing by ARESP
        ind_y_aresp(:, pk.ind_x_start_aresp) ...
                            = pk.ind_y_start_aresp'; % assign starting y index at starting x index
        
        % loop through each starter y index at pk.ind_x_start and propagate from that x index using the ARESP slope
        for ii = 1:pk.num_aresp %#ok<*FXUP>
            for jj = (pk.ind_x_start_aresp - 1):-1:1 % loop heading left of pk.ind_x_start
                ind_y_aresp(ii, jj)     = ind_y_aresp(ii, (jj + 1)) + block.slope_aresp(round(ind_y_aresp(ii, (jj + 1))), (jj + 1));
                if (ind_y_aresp(ii, jj) < 1)
                    ind_y_aresp(ii, jj) = 1;
                elseif (ind_y_aresp(ii, jj) > num_sample_trim)
                    ind_y_aresp(ii, jj) = num_sample_trim;
                elseif isnan(ind_y_aresp(ii, jj))
                    break
                end
            end
            for jj = (pk.ind_x_start_aresp + 1):block.num_trace % loop heading right of pk.ind_x_start
                ind_y_aresp(ii, jj)     = ind_y_aresp(ii, (jj - 1)) - block.slope_aresp(round(ind_y_aresp(ii, (jj - 1))), (jj - 1));
                if (ind_y_aresp(ii, jj) < 1)
                    ind_y_aresp(ii, jj) = 1;
                elseif (ind_y_aresp(ii, jj) > num_sample_trim)
                    ind_y_aresp(ii, jj) = num_sample_trim;
                elseif isnan(ind_y_aresp(ii, jj))
                    break
                end
            end
        end
        
        % remove ARESP layers that NaN'd out
        tmp1                = find(sum(isnan(ind_y_aresp), 2));
        if ~isempty(tmp1)
            ind_y_aresp     = ind_y_aresp(setdiff(1:pk.num_aresp, tmp1), :);
            pk.ind_y_start_aresp ...
                            = pk.ind_y_start_aresp(setdiff(1:pk.num_aresp, tmp1));
            pk.num_aresp    = pk.num_aresp - length(tmp1);
        end
    end

%% Pick ARESP-propagated layers to keep for flattening

    function pk_keep_aresp(source, eventdata)
        
        if ~aresp_done
            set(status_box, 'string', 'No ARESP-tracked layers yet.')
            return
        end
        if ~any(strcmp(disp_type, {'twtt' '~depth'}))
            set(status_box, 'string', 'ARESP layers can only be kept in twtt or ~depth.')
            return
        end
        if ~get(aresp_check, 'value')
            set(aresp_check, 'value', 1)
            show_aresp
            keep_aresp_done = false;
        end
        
        set(status_box, 'string', 'Pick best ARESP-traced layers for flattening...')
        pause(0.1)
        set(status_box, 'string', 'Left-click: select; return/enter: keep; D: unkeep; Q: stop...')
        
        axes(ax_radar)
        ii                  = 0;
        
        set(pkgui, 'keypressfcn', [], 'windowbuttondownfcn', [])
        
        while true
            
            [ind_x_pk, ind_y_pk, button] ...
                            = ginput(1); % point near layer to keep
            if ~ii % dummy start
                ii          = 1;
                continue
            end
            
            if ~any(strcmpi(char(button), {'Q' 'D'}))
                if (button ~= 1)
                    continue
                end
            end
            
            if strcmpi(char(button), 'Q') % quit/stop picking keepers
                
                break
                
            elseif strcmpi(char(button), 'D') % delete a keeper
                
                if pk.num_keep_aresp
                    ind_x_pk= interp1(block.dist_lin(ind_decim), ind_decim, ind_x_pk, 'nearest', 'extrap');
                    switch disp_type
                        case 'twtt'
                            ind_y_pk = interp1(block.twtt, 1:num_sample_trim, (1e-6 * ind_y_pk), 'nearest', 'extrap');
                        case '~depth'
                            ind_y_pk = interp1(block.twtt, 1:num_sample_trim, (1e-6 * (ind_y_pk + (1e6 * (block.twtt_surf(ind_x_pk) - block.twtt(1))))), 'nearest', 'extrap');
                    end
                    tmp1    = pk.ind_keep_aresp(~isnan(ind_y_aresp(pk.ind_keep_aresp, ind_x_pk)));
                    if (length(tmp1) > 1)
                        tmp1= interp1(ind_y_aresp(tmp1, ind_x_pk), tmp1, ind_y_pk, 'nearest', 'extrap');
                    end
                    set(p_aresp(tmp1), 'color', 'c', 'linewidth', 1)
                    if depth_avail
                        set(p_arespdepth(tmp1), 'color', 'c', 'linewidth', 1)
                    end
                    if (pk.num_keep_aresp > 1)
                        tmp2=interp1(pk.ind_keep_aresp, 1:pk.num_keep_aresp, tmp1, 'nearest', 'extrap');
                    else
                        tmp2= 1;
                    end
                    pk.ind_keep_aresp ...
                            = pk.ind_keep_aresp(setdiff(1:pk.num_keep_aresp, tmp2));
                    pk.num_keep_aresp ...
                            = pk.num_keep_aresp - 1;
                    set(status_box, 'string', 'Layer unkept. Pick again...')
                end
                
            else
                
                ind_x_pk    = interp1(block.dist_lin(ind_decim), ind_decim, ind_x_pk, 'nearest', 'extrap'); % convert x index to non-decimated range
                tmp2        = find(~isnan(ind_y_aresp(:, ind_x_pk)));
                tmp1        = ind_y_aresp(tmp2, ind_x_pk);
                [tmp1, tmp3]= unique(tmp1);
                if (length(tmp1) > 1)
                    switch disp_type
                        case 'twtt'
                            ind_y_pk = interp1(tmp1, tmp2(tmp3), interp1(block.twtt, 1:num_sample_trim, (1e-6 * ind_y_pk), 'nearest', 'extrap'), 'nearest', 'extrap'); % index of nearest layer
                        case '~depth'
                            ind_y_pk = interp1(tmp1, tmp2(tmp3), interp1(block.twtt, 1:num_sample_trim, (1e-6 * (ind_y_pk + (1e6 * (block.twtt_surf(ind_x_pk) - block.twtt(1))))), 'nearest', 'extrap'), 'nearest', 'extrap');
                    end
                else
                    ind_y_pk= tmp2;
                end
                if any(ind_y_pk == pk.ind_keep_aresp)
                    set(status_box, 'string', 'Chosen layer already kept. Pick again...')
                    continue
                end
                set(p_aresp(ind_y_pk), 'color', 'w') % change layer color to white
                if depth_avail
                    set(p_arespdepth(ind_y_pk), 'color', 'w')
                end
                
                waitforbuttonpress
                
                if (double(get(pkgui, 'currentcharacter')) == 13)
                    set(p_aresp(ind_y_pk), 'linewidth', 2) % increase linewidth
                    if depth_avail
                        set(p_arespdepth(ind_y_pk), 'linewidth', 2)
                    end
                    pk.ind_keep_aresp(pk.num_keep_aresp + 1) ...
                            = ind_y_pk; % preserve layer number
                    pk.num_keep_aresp ...
                            = pk.num_keep_aresp + 1;
                else
                    set(p_aresp(ind_y_pk), 'color', 'c') % return to previous color
                    if depth_avail
                        set(p_arespdepth(ind_y_pk), 'color', 'c')
                    end
                    set(status_box, 'string', 'Then pick again...')
                end
            end
        end
        
        if pk.num_keep_aresp
            pk.ind_keep_aresp ...
                            = sort(pk.ind_keep_aresp); % sort kept layer numbers
            set(status_box, 'string', ['Picked ' num2str(pk.num_keep_aresp) ' ARESP-tracked layers to keep.'])
            keep_aresp_done = true;
        else
            set(status_box, 'string', 'No layers kept.')
        end
        
        set(pkgui, 'keypressfcn', @keypress, 'windowbuttondownfcn', @mouse_click)
    end

%% Manually pick a layer to use for flattening

    function track_man(source, eventdata)
        
        if ~any(strcmp(disp_type, {'twtt' '~depth'}))
            set(status_box, 'string', 'Layers can only be tracked in twtt or ~depth.')
            return
        end
        if ~trim_done
            set(status_box, 'string', 'Cannot pick manual layers for flattening if trimming not done.')
            return
        end
        
        axes(ax_radar)
        ii                  = 0;
        [ind_x_pk, ind_y_pk, button] ...
                            = deal([]);
        if all(~p_man)
            p_man           = zeros(0, 2);
        end
        if depth_avail
            if all(~p_mandepth)
                p_mandepth  = zeros(0, 2);
            end
        end
        
        tmp4                = 0;
        
        set(pkgui, 'keypressfcn', [], 'windowbuttondownfcn', [])
        
        while true
            
            set(status_box, 'string', 'Left-click: pick; U: undo pick; D: delete layer; return: done; Q: quit...')
            [tmp1, tmp2, button] ...
                            = ginput(1); % pick manually, press enter when done
            if ~ii % dummy initialization
                ii          = 1;
                continue
            end
            
            if (button == 1)
                
                [ind_x_pk(ii), ind_y_pk(ii)] ...
                            = deal(tmp1, tmp2);
                ind_x_pk(ii)= interp1(block.dist_lin(ind_decim), ind_decim, ind_x_pk(ii), 'nearest', 'extrap'); % raw picks must be indices, not dimensionalized vectors (horizontal)
                switch disp_type
                    case 'twtt'
                        ind_y_pk(ii) ...
                            = interp1(block.twtt, 1:num_sample_trim, (1e-6 .* ind_y_pk(ii)), 'nearest', 'extrap');
                    case '~depth'
                        ind_y_pk(ii) ...
                            = interp1(block.twtt, 1:num_sample_trim, (1e-6 .* (ind_y_pk(ii) + (1e6 * (block.twtt_surf(ind_x_pk(ii)) - block.twtt(1))))), 'nearest', 'extrap');
                end
                if (ii > 1)
                    delete(p_man((pk.num_man + 1), 1))
                    if depth_avail
                        delete(p_mandepth((pk.num_man + 1), 1))
                    end
                end
                p_man((pk.num_man + 1), 1) ...
                            = plot(block.dist_lin(ind_x_pk), (1e6 .* block.twtt(ind_y_pk)), 'wx', 'markersize', 12, 'visible', 'off'); % original picks
                if depth_avail
                    p_mandepth((pk.num_man + 1), 1) ...
                            = plot(block.dist_lin(ind_x_pk), (1e6 .* block.twtt(ind_y_pk - ind_surf(ind_x_pk) + 1)), 'wx', 'markersize', 12, 'visible', 'off');
                end
                switch disp_type
                    case 'twtt'
                        set(p_man(logical(p_man) & ishandle(p_man)), 'visible', 'on')
                    case '~depth'
                        set(p_mandepth(logical(p_mandepth) & ishandle(p_mandepth)), 'visible', 'on')
                end
                ii          = ii + 1;
                tmp4        = button;
                
            elseif (strcmpi(char(button), 'U') && (ii > 1))
                
                [ind_x_pk, ind_y_pk] ...
                            = deal(ind_x_pk(1:(end - 1)), ind_y_pk(1:(end - 1)));
                delete(p_man((pk.num_man + 1), 1))
                if depth_avail
                    delete(p_mandepth((pk.num_man + 1), 1))
                end
                if (ii > 2)
                    p_man((pk.num_man + 1), 1) ...
                            = plot(block.dist_lin(ind_x_pk), (1e6 .* block.twtt(ind_y_pk)), 'wx', 'markersize', 12, 'visible', 'off'); % original picks
                    if depth_avail
                        p_mandepth((pk.num_man + 1), 1) ...
                            = plot(block.dist_lin(ind_x_pk), (1e6 .* block.twtt(ind_y_pk - ind_surf(ind_x_pk) + 1)), 'wx', 'markersize', 12, 'visible', 'off');
                    end
                end
                switch disp_type
                    case 'twtt'
                        set(p_man(logical(p_man) & ishandle(p_man)), 'visible', 'on')
                    case '~depth'
                        set(p_mandepth(logical(p_man) & ishandle(p_man)), 'visible', 'on')
                end
                ii          = ii - 1;
                tmp4        = button;
                
            elseif strcmpi(char(button), 'D') % delete a manual layer
                
                if (tmp4 == 1)
                    set(status_box, 'string', 'Cannot delete a manual layer in the middle of picking one...')
                    pause(0.1)
                    continue
                end
                if pk.num_man
                    tmp1    = interp1(block.dist_lin(ind_decim), ind_decim, tmp1, 'nearest', 'extrap'); % raw picks must be indices, not dimensionalized vectors (horizontal)
                    switch disp_type
                        case 'twtt'
                            tmp2 = interp1(block.twtt, 1:num_sample_trim, (1e-6 .* tmp2), 'nearest', 'extrap');
                        case '~depth'
                            tmp2 = interp1(block.twtt, 1:num_sample_trim, (1e-6 .* (tmp2 + (1e6 * (block.twtt_surf(tmp1) - block.twtt(1))))), 'nearest', 'extrap');
                    end
                    tmp3    = NaN(1, pk.num_man);
                    for ii = 1:pk.num_man
                        tmp3(ii) ...
                            = pk.ind_y_man(ii, tmp1);
                    end
                    if (pk.num_man > 1)
                        tmp3= interp1(tmp3, 1:pk.num_man, tmp2, 'nearest', 'extrap');
                    else
                        tmp3= 1;
                    end
                    pk.ind_y_man ...
                            = pk.ind_y_man(setdiff(1:pk.num_man, tmp3), :);
                    delete(p_man(tmp3, :))
                    if depth_avail
                        delete(p_mandepth(tmp3, :))
                    end
                    p_man   = p_man(setdiff(1:pk.num_man, tmp3), :);
                    if depth_avail
                        p_mandepth ...
                            = p_mandepth(setdiff(1:pk.num_man, tmp3), :);
                    end
                    pk.num_man ...
                            = pk.num_man - 1;
                    set(status_box, 'string', ['Deleted manual layer #' num2str(tmp3) '.'])
                    pause(0.1)
                end
                
            elseif (double(get(pkgui, 'currentcharacter')) == 13)
                
                if (length(ind_x_pk) < 3)
                    if (any(ishandle(p_man(end, :))) && any(p_man(end, :)))
                        delete(p_man(end, (ishandle(p_man(end, :)) & logical(p_man(end, :)))))
                    end
                    if (any(ishandle(p_mandepth(end, :))) && any(p_mandepth(end, :)))
                        delete(p_mandepth(end, (ishandle(p_mandepth(end, :)) & logical(p_mandepth(end, :)))))
                    end
                    p_man   = p_man(1:(end - 1), :);
                    if depth_avail
                        p_mandepth ...
                            = p_mandepth(1:(end - 1), :);
                    end
                    set(status_box, 'string', 'Not enough picked points to make a manual layer. Start over.')
                    pause(0.1)
                    ii      = 1;
                    continue
                end
                pk.num_man  = pk.num_man + 1;
                
                if ~issorted(ind_x_pk) % resort picks
                    [ind_x_pk, tmp1] ...
                            = sort(ind_x_pk);
                    ind_y_pk= ind_y_pk(tmp1);
                end
                if (length(unique(ind_x_pk)) < length(ind_x_pk)) % don't keep picks that are accidentally at the same horizontal index, otherwise spline will fail
                    [ind_x_pk, tmp1] ...
                            = unique(ind_x_pk);
                    ind_y_pk= ind_y_pk(tmp1);
                end
                if (ind_x_pk(1) ~= 1) % a pick at the very beginning of the chunk is necessary for interpolation
                    ind_y_pk= [interp1(ind_x_pk, ind_y_pk, 1, 'linear', 'extrap') ind_y_pk]; %#ok<AGROW>
                    ind_x_pk= [1 ind_x_pk]; %#ok<AGROW>
                end
                if (ind_x_pk(end) ~= block.num_trace) % a pick at the very end is also necessary
                    ind_y_pk= [ind_y_pk interp1(ind_x_pk, ind_y_pk, block.num_trace, 'linear', 'extrap')]; %#ok<AGROW>
                    ind_x_pk= [ind_x_pk block.num_trace]; %#ok<AGROW>
                end
                
                switch disp_type
                    case 'twtt'
                        pk.ind_y_man(pk.num_man, 1:block.num_trace) ...
                            = spline(ind_x_pk, ind_y_pk, 1:block.num_trace); % interpolate spline through picks
                    case '~depth'
                        pk.ind_y_man(pk.num_man, 1:block.num_trace) ...
                            = spline(ind_x_pk, (ind_y_pk - ind_surf(ind_x_pk) + 1), 1:block.num_trace) + ind_surf - 1; % interpolate spline through picks
                end
                
                p_man(pk.num_man, 2) ...
                            = plot(block.dist_lin(ind_decim), (1e6 .* block.twtt(round(pk.ind_y_man(pk.num_man, ind_decim)))), 'linewidth', 2, 'color', [0.85 0.85 0.85], 'visible', 'off'); % max-picking spline
                if depth_avail
                    tmp1    = ind_decim(~isnan(pk.ind_y_man(pk.num_man, ind_decim)) & ~isnan(ind_surf(ind_decim)));
                    tmp2    = round(pk.ind_y_man(pk.num_man, tmp1) - ind_surf(tmp1) + 1);
                    tmp2((tmp2 < 1) | (tmp2 > num_sample_trim)) ...
                            = NaN;
                    p_mandepth(pk.num_man, 2) ...
                            = plot(block.dist_lin(tmp1(~isnan(tmp2))), (1e6 .* block.twtt(tmp2(~isnan(tmp2)))), 'linewidth', 2, 'color', [0.85 0.85 0.85], 'visible', 'off');
                end
                
                set(man_check, 'value', 1)
                show_man
                ii          = 1;
                [ind_x_pk, ind_y_pk] ...
                            = deal([]);
                set(status_box, 'string', 'Manual layer successfully picked. Now pick another...')
                pause(0.1)
                tmp4        = button;
                
            elseif strcmpi(char(button), 'Q')
                
                if (size(p_man, 1) > pk.num_man)
                    if (any(p_man(end, :)) && any(ishandle(p_man(end, :))))
                        delete(p_man(end, (logical(p_man(end, :)) & ishandle(p_man(end, :)))))
                    end
                    if (any(p_mandepth(end, :)) && any(ishandle(p_mandepth(end, :))))
                        delete(p_mandepth(end, (logical(p_mandepth(end, :)) & ishandle(p_mandepth(end, :)))))
                    end
                    p_man   = p_man(1:(end - 1), :);
                    if depth_avail
                        p_mandepth ...
                            = p_mandepth(1:(end - 1), :);
                    end
                end
                set(status_box, 'string', 'Ended manual picking.')
                set(pkgui, 'keypressfcn', @keypress, 'windowbuttondownfcn', @mouse_click)
                return
            end
        end
    end

%% Flatten data

    function flatten(source, eventdata)
        
        switch pk.predict_or_pk
            case 'predict'
                if ((pk.num_keep_phase + pk.num_keep_aresp + pk.num_man) < ord_poly)
                    set(status_box, 'string', 'Not enough predicted layers to flatten (need 2-3+ including surface).')
                    return
                end
            case 'pk'
                if ~all(smooth_done)
                    set(status_box, 'string', 'Not all layers have been smoothed prior to flattening.')
                    return
                end
                if (pk.num_layer < ord_poly)
                    set(status_box, 'string', 'Not enough layers to attempt flattening.')
                    return
                end
        end
        
        if (logical(p_surfflat) && ishandle(p_surfflat))
            delete(p_surfflat)
        end
        if (logical(p_bedflat) && ishandle(p_bedflat))
            delete(p_bedflat)
        end
        
        flat_done           = false;
        
        set(status_box, 'string', 'Starting polynomial fitting...')
        pause(0.1)
        
        pk.poly_flat        = NaN((ord_poly + 1), block.num_trace, 'single');
        % 2nd-order polynomial fits between kept layers at pk.ind_x_start (x) and kept layers at each trace (y)
        
        if parallel_check
            pctRunOnAll warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale')
            pctRunOnAll warning('off', 'MATLAB:polyfit:PolyNotUnique')
        else
            warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale')
            warning('off', 'MATLAB:polyfit:PolyNotUnique')
        end
        
        switch pk.predict_or_pk
            
            case 'predict' % concatenate predicted layers to use in flattening, then sort them all together
                
                if surf_avail
                    tmp1    = ind_surf;
                else
                    tmp1    = [];
                end
                if pk.num_keep_phase
                    tmp1    = [tmp1; ind_y_phase(pk.ind_keep_phase, :)];
                end
                if pk.num_keep_aresp
                    tmp1    = [tmp1; ind_y_aresp(pk.ind_keep_aresp, :)];
                end
                if pk.num_man
                    tmp1    = [tmp1; pk.ind_y_man];
                end
                
                [~, tmp2]   = sort(nanmean(tmp1, 2));
                tmp1        = tmp1(tmp2, :); % all layers sorted by mean value
                tmp2        = tmp1(:, 1);
                ind_x_pk    = 1;
                
            case 'pk' % flattening using picked layers
                
                ind_y_curr  = NaN(pk.num_layer, block.num_trace);
                for ii = 1:pk.num_layer
                    ind_y_curr(ii, :) ...
                            = pk.layer(ii).ind_y_smooth;
                end
                if surf_avail
                    ind_y_curr ...
                            = [ind_surf; ind_y_curr];
                end
                
                ind_x_pk    = sum(isnan(tmp1));
                ind_x_pk    = find((ind_x_pk == min(ind_x_pk)), 1); % first trace in the record that has the maximum number of layers, i.e., minimum number of NaNs
                ind_y_pk    = ind_y_curr(:, ind_x_pk);
                tmp1        = ind_y_curr(~isnan(ind_y_pk), :); % depth of all layers that are not nan at the reference trace
                tmp2        = tmp1(:, ind_x_pk); % depths of non-nan layers at reference trace
        end
        
        tmp3                = find(sum(~isnan(tmp1)) > ord_poly); % traces where it will be worth doing the polynomial
        
        % polyfit to layers
        if parallel_check
            tmp1            = tmp1(:, tmp3);
            tmp4            = pk.poly_flat(:, tmp3);
            parfor ii = 1:length(tmp3)
                tmp4(:, ii) = polyfit(tmp2(~isnan(tmp1(:, ii))), tmp1(~isnan(tmp1(:, ii)), ii), ord_poly)'; %#ok<PFBNS>
            end
            pk.poly_flat(:, tmp3) ...
                            = tmp4;
            tmp4            = 0;
        else
            for ii = 1:length(tmp3)
                pk.poly_flat(:, tmp3(ii)) ...
                            = polyfit(tmp2(~isnan(tmp1(:, tmp3(ii)))), tmp1(~isnan(tmp1(:, tmp3(ii))), tmp3(ii)), ord_poly)';
            end
        end
        tmp4                = 0;
        
        % smooth polynomials
        tmp1                = round((2 * pk.length_smooth) / nanmean(diff(block.dist_lin)));
        if parallel_check
            tmp4            = pk.poly_flat;
            parfor ii = 1:(ord_poly + 1)
                tmp4(ii, :) = smooth_lowess(tmp4(ii, :), tmp1);
            end
            pk.poly_flat    = tmp4;
            tmp4            = 0;
        else
            for ii = 1:(ord_poly + 1)
                pk.poly_flat(ii, :) ...
                            = smooth_lowess(pk.poly_flat(ii, :), tmp1);
            end
        end
        
        % iterate using other layers if any are available
        if (strcmp(pk.predict_or_pk, 'pk') && any(isnan(ind_y_pk)))
            
            % determine which layers overlap with original polyfit layers and order polyfit iteration based on the length of their overlap
            tmp1            = zeros(1, length(ind_y_pk));
            for ii = find(isnan(ind_y_pk))'
                tmp1(ii)    = length(find(~isnan(ind_y_curr(ii, tmp3))));
            end
            [tmp1, tmp2]    = sort(tmp1, 'descend');
            tmp5            = tmp2(logical(tmp1));
            set(status_box, 'string', ['Iterating flattening for ' num2str(length(tmp2(logical(tmp1)))) ' overlapping layers out of ' num2str(length(find(isnan(ind_y_pk)))) ' not fitted initially...'])
            
            pause(0.1)
            
            kk              = 0;
            
            for ii = tmp5
                
                kk                  = kk + 1;
                set(status_box, 'string', ['Adding layer #' num2str(ii) ' (' num2str(kk) ' / ' num2str(length(tmp5)) ') to flattening...'])
                pause(0.1)
                
                % calculate flattening matrix based on current polynomials
                tmp1                = find(~isnan(ind_y_curr(ii, :)));
                ind_y_mat           = single((1:num_sample_trim)');
                ind_y_mat           = ind_y_mat(:, ones(1, length(tmp1))); % matrix of y indices
                switch ord_poly
                    case 2
                        ind_y_flat  = ((ind_y_mat .^ 2) .* pk.poly_flat(ones(num_sample_trim, 1), tmp1)) + (ind_y_mat .* (pk.poly_flat((2 .* ones(num_sample_trim, 1)), tmp1))) + pk.poly_flat((3 .* ones(num_sample_trim, 1)), tmp1);
                    case 3
                        ind_y_flat  = ((ind_y_mat .^ 3) .* pk.poly_flat(ones(num_sample_trim, 1), tmp1)) + ((ind_y_mat .^ 2) .* pk.poly_flat((2 .* ones(num_sample_trim, 1)), tmp1)) + ...
                                      (ind_y_mat .* (pk.poly_flat((3 .* ones(num_sample_trim, 1)), tmp1))) + pk.poly_flat((4 .* ones(num_sample_trim, 1)), tmp1);
                end
                ind_y_mat           = 0;
                ind_y_flat(ind_y_flat < 1) ...
                                    = 1; % limit too-low indices
                ind_y_flat(ind_y_flat > num_sample_trim) ...
                                    = num_sample_trim; % limit too-high indices
                tmp3                = find(sum(~isnan(ind_y_flat)));
                
                % find best-fit value at original trace
                tmp4                = NaN(1, length(tmp1));
                tmp3                = tmp3(~isnan(ind_y_curr(ii, tmp1(tmp3)))); % indices where both flattening and new layer exist
                
                if isempty(tmp3)
                    set(status_box, 'string', ['Layer #' num2str(ii) ' has no overlap...'])
                    pause(0.1)
                    continue
                end
                
                % index of current layer at reference trace for all overlapping traces
                for jj = tmp3
                    [~, tmp2]       = unique(ind_y_flat(:, jj), 'last');
                    tmp2            = intersect((1 + find(diff(ind_y_flat(:, jj)) > 0)), tmp2);
                    if ~isempty(tmp2)
                        tmp4(jj)    = interp1(ind_y_flat(tmp2, jj), tmp2, ind_y_curr(ii, tmp1(jj)), 'linear', NaN);
                    end
                end
                
                ind_y_pk(ii)        = nanmean(tmp4); % best guess y index at reference trace 
                
                % extract best layers again, now including the new layer
                tmp4                = tmp1(tmp3);
                tmp1                = ind_y_curr(~isnan(ind_y_pk), tmp4);
                tmp2                = ind_y_pk(~isnan(ind_y_pk));
                tmp3                = find(sum(~isnan(tmp1)) > ord_poly);
                tmp4                = tmp4(sum(~isnan(tmp1)) > ord_poly);
                
                % new polynomials using additional "depth" for this layer
                if parallel_check
                    tmp1            = tmp1(:, tmp3);
                    tmp3            = pk.poly_flat(:, tmp3);
                    parfor jj = 1:length(tmp4)
                        tmp3(:, jj) = polyfit(tmp2(~isnan(tmp1(:, jj))), tmp1(~isnan(tmp1(:, jj)), jj), ord_poly)'; %#ok<PFBNS>
                    end
                    pk.poly_flat(:, tmp4) ...
                                    = tmp3;
                else
                    for jj = 1:length(tmp3)
                        pk.poly_flat(:, tmp4(jj)) ...
                                    = polyfit(tmp2(~isnan(tmp1(:, tmp3(jj)))), tmp1(~isnan(tmp1(:, tmp3(jj))), tmp3(jj)), ord_poly)';
                    end
                end
                
                % smooth polynomials again
                tmp2                = round(mean(pk.length_smooth) / nanmean(diff(block.dist_lin(ind_decim))));
                if parallel_check
                    tmp1            = pk.poly_flat;
                    parfor jj = 1:(ord_poly + 1)
                        tmp1(jj, :) = smooth_lowess(tmp1(jj, :), tmp2);
                    end
                    pk.poly_flat    = tmp1;
                    tmp1            = 0;
                else
                    for jj = 1:(ord_poly + 1)
                        pk.poly_flat(jj, :) ...
                                    = smooth_lowess(pk.poly_flat(jj, :), tmp2);
                    end
                end
            end
            if surf_avail
                ind_y_curr  = ind_y_curr(2:end, :);
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
        ind_y_mat           = single((1:num_sample_trim)');
        ind_y_mat           = ind_y_mat(:, ones(1, block.num_trace)); % matrix of y indices
        switch ord_poly
            case 2
                ind_y_flat  = ((ind_y_mat .^ 2) .* pk.poly_flat(ones(num_sample_trim, 1), :)) + (ind_y_mat .* (pk.poly_flat((2 .* ones(num_sample_trim, 1)), :))) + pk.poly_flat((3 .* ones(num_sample_trim, 1)), :);
            case 3
                ind_y_flat  = ((ind_y_mat .^ 3) .* pk.poly_flat(ones(num_sample_trim, 1), :)) + ((ind_y_mat .^ 2) .* pk.poly_flat((2 .* ones(num_sample_trim, 1)), :)) + ...
                              (ind_y_mat .* (pk.poly_flat((3 .* ones(num_sample_trim, 1)), :))) + pk.poly_flat((4 .* ones(num_sample_trim, 1)), :);
        end
        ind_y_mat           = 0;
        ind_y_flat(ind_y_flat < 1) ...
                            = 1;
        ind_y_flat(ind_y_flat > num_sample_trim) ...
                            = num_sample_trim;
        
        set(status_box, 'string', 'Done polynomial fitting. Now flattening radargram...')
        pause(0.1)
        
        % flattened radargram based on layer fits
        amp_flat            = NaN(size(block.amp, 1), block.num_trace, 'single');
        tmp2                = find(sum(~isnan(ind_y_flat)));
        if parallel_check
            pctRunOnAll warning('off', 'MATLAB:interp1:NaNinY')
            tmp1            = block.amp(:, tmp2);
            tmp3            = ind_y_flat(:, tmp2);
            tmp4            = amp_flat(:, tmp2);
            parfor ii = 1:length(tmp2)
                tmp4(:, ii) = interp1(tmp1(:, ii), tmp3(:, ii));
            end
            amp_flat(:, tmp2) ...
                            = tmp4;
            tmp1            = 0;
            pctRunOnAll warning('on', 'MATLAB:interp1:NaNinY')
        else
            warning('off', 'MATLAB:interp1:NaNinY')
            for ii = 1:length(tmp2)
                amp_flat(:, tmp2(ii)) ...
                            = interp1(block.amp(:, ii), ind_y_flat(:, ii));
            end
            warning('on', 'MATLAB:interp1:NaNinY')
        end
        
        set(status_box, 'string', 'Flattening layers...')
        pause(0.1)
                
        % flatten surface pick
        ind_surf_flat       = NaN(1, block.num_trace);
        if surf_avail
            for ii = tmp2
                [~, tmp1]   = unique(ind_y_flat(:, ii));
                if (length(tmp1) > 1)
                    ind_surf_flat(ii) ...
                            = interp1(ind_y_flat(tmp1, ii), tmp1, ind_surf(ii), 'nearest', 'extrap');
                end
            end
            ind_surf_flat   = round(ind_surf_flat);
            ind_surf_flat((ind_surf_flat < 1) | (ind_surf_flat > num_sample_trim)) ...
                            = NaN;
        end
        
        % flatten bed pick
        ind_bed_flat        = NaN(1, block.num_trace);
        if bed_avail
            for ii = tmp2
                [~, tmp1]   = unique(ind_y_flat(:, ii));
                if (length(tmp1) > 1)
                    ind_bed_flat(ii) ...
                            = interp1(ind_y_flat(tmp1, ii), tmp1, ind_bed(ii), 'nearest', 'extrap');
                end
            end
            ind_bed_flat    = round(ind_bed_flat);
            ind_bed_flat((ind_bed_flat < 1) | (ind_bed_flat > num_sample_trim)) ...
                            = NaN;
        end
        
        flat_done           = true;
        smooth_done         = false(1, pk.num_layer);
        
        % re-flatten layers if any exist
        if (pk_done && pk.num_layer)
            warning('off', 'MATLAB:interp1:NaNinY')
            ind_y_curr      = NaN(pk.num_layer, block.num_trace);
            for ii = 1:pk.num_layer
                ind_y_curr(ii, :) ...
                            = pk.layer(ii).ind_y;
            end
            [ind_y_flat_mean, ind_y_flat_smooth] ...
                            = deal(NaN(pk.num_layer, num_decim_flat));
            for ii = find(~sum(isnan(ind_y_flat(:, ind_decim_flat))))
                [~, tmp1]   = unique(ind_y_flat(:, ind_decim_flat(ii)));
                ind_y_flat_mean(~isnan(ind_y_curr(:, ind_decim_flat(ii))), ii) ...
                            = interp1(ind_y_flat(tmp1, ind_decim_flat(ii)), tmp1, ind_y_curr(~isnan(ind_y_curr(:, ind_decim_flat(ii))), ind_decim_flat(ii)), 'nearest', 'extrap');
            end
            warning('on', 'MATLAB:interp1:NaNinY')
            pk_smooth
        end
        
        set(disp_check(5), 'visible', 'on')
        mean_flat
        pk_select
        set(status_box, 'string', 'Flattened radargram and layers.')
    end

%% Horizontally average flattened data

    function mean_flat(source, eventdata)
        if ~flat_done
            set(status_box, 'string', 'No flattened data to average yet.')
            return
        end
        
        set(status_box, 'string', 'Horizontally averaging flattened data...')
        pause(0.1)
        
        ind_decim_flat_old  = ind_decim_flat;
        
        if (decim_flat > 1)
            ind_decim_flat  = ceil((decim_flat / 2) + 1):decim_flat:(block.num_trace - ceil(decim_flat / 2));
            num_decim_flat  = length(ind_decim_flat);
            amp_flat_mean   = NaN(num_sample_trim, num_decim_flat, 'single');
            tmp1            = floor(decim_flat / 2);
            for ii = 1:num_decim_flat
                amp_flat_mean(:, ii) ...
                            = nanmean(amp_flat(:, (ind_decim_flat(ii) - tmp1):(ind_decim_flat(ii) + tmp1)), 2);
            end
        elseif (decim_flat == 1)
            ind_decim_flat  = 1:block.num_trace;
            num_decim_flat  = block.num_trace;
            amp_flat_mean   = amp_flat;
        end
        if (logical(p_bedflat) && ishandle(p_bedflat))
            delete(p_bedflat)
        end
        if (logical(p_surfflat) && ishandle(p_surfflat))
            delete(p_surfflat)
        end
        if (any(p_pkflat) && any(ishandle(p_pkflat)))
            delete(p_pkflat(logical(p_pkflat) & ishandle(p_pkflat)))
        end
        if (any(p_pksmoothflat) && any(ishandle(p_pksmoothflat)))
            delete(p_pksmoothflat(logical(p_pksmoothflat) & ishandle(p_pksmoothflat)))
        end
        [p_pkflat, p_pksmoothflat] ...
                            = deal(zeros(1, pk.num_layer));
        [tmp1, tmp2]        = deal(ind_y_flat_mean, ind_y_flat_smooth);
        [ind_y_flat_mean, ind_y_flat_smooth] ...
                            = deal(NaN(pk.num_layer, num_decim_flat));
        if surf_avail
            if ~all(isnan(ind_surf_flat(ind_decim_flat)))
                tmp3        = ind_surf_flat(ind_decim_flat);
                p_surfflat  = plot(block.dist_lin(ind_decim_flat(~isnan(tmp3))), (1e6 .* block.twtt(tmp3(~isnan(tmp3)))), 'm.', 'markersize', 12, 'visible', 'off');
            end
        end
        if bed_avail
            if ~all(isnan(ind_bed_flat(ind_decim_flat)))
                tmp3        = ind_bed_flat(ind_decim_flat);
                p_bedflat   = plot(block.dist_lin(ind_decim_flat(~isnan(tmp3))), (1e6 .* block.twtt(tmp3(~isnan(tmp3)))), 'm.', 'markersize', 12, 'visible', 'off');
            end
        end
        warning('off', 'MATLAB:interp1:NaNinY')
        for ii = 1:pk.num_layer
            ind_y_flat_mean(ii, :) ...
                            = round(interp1(ind_decim_flat_old, tmp1(ii, :), ind_decim_flat, 'linear', 'extrap'));
            if all(isnan(ind_y_flat_mean(ii, :)))
                continue
            end
            p_pkflat(ii)    = plot(block.dist_lin(ind_decim_flat(~isnan(ind_y_flat_mean(ii, :)))), (1e6 .* block.twtt(ind_y_flat_mean(ii, ~isnan(ind_y_flat_mean(ii, :))))), 'r.', 'markersize', 12, 'visible', 'off');
            if smooth_done(ii)
                ind_y_flat_smooth(ii, :) ...
                            = round(interp1(ind_decim_flat_old, tmp2(ii, :), ind_decim_flat, 'linear', 'extrap'));
                if ~isempty(find(~isnan(ind_y_flat_smooth(ii, :)), 1))
                    p_pksmoothflat(ii) ...
                            = plot(block.dist_lin(ind_decim_flat(~isnan(ind_y_flat_smooth(ii, :)))), (1e6 .* block.twtt(ind_y_flat_smooth(ii, ~isnan(ind_y_flat_smooth(ii, :))))), 'g.', 'markersize', 12, 'visible', 'off');
                end
            end
        end
        warning('on', 'MATLAB:interp1:NaNinY')

        pk_smooth
       
        set(disp_group, 'selectedobject', disp_check(5))
        disp_type           = 'flat';
        plot_flat
        set(status_box, 'string', 'Horizontally averaged and flattened radargram.')
    end

%% Semi-automatic layer picking

    function pk_auto(source, eventdata)
        
        if ~any(strcmp(disp_type, {'twtt' '~depth' 'flat'}))
            set(status_box, 'string', 'Layers can only be traced in twtt, ~depth or flat.')
            return
        end
        
        axes(ax_radar)
        
        if ~load_done
            set(status_box, 'string', 'Data not yet loaded.')
            return
        end
        
        if ~pk.num_layer % initialize for no layers
            smooth_done     = logical([]);
            [p_pksmooth, p_pksmoothdepth, p_pksmoothflat, ind_y_flat_mean, ind_y_flat_smooth] ...
                            = deal([]);
        end
        
        % layer picking callbacks
        tmp3                = pk.num_layer + 1; % used to shorten real/flatten conversion loop after this while picking loop
        set(pkgui, 'keypressfcn', [], 'windowbuttondownfcn', [])
        
        while true
            
            set(status_box, 'string', 'Left-click: Pick; D: delete; U: undo; L: cut left; R: cut right; C: cut chunk; M: merge; Q: done...')
            
            % get pick and convert to indices
            [ind_x_pk, ind_y_pk, button] ...
                            = ginput(1);
            switch disp_type
                case {'twtt' '~depth'}
                    ind_x_pk= interp1(block.dist_lin(ind_decim), ind_decim, ind_x_pk, 'nearest', 'extrap');
                case 'flat'
                    ind_x_pk= interp1(block.dist_lin(ind_decim_flat), 1:num_decim_flat, ind_x_pk, 'nearest', 'extrap');
            end
            switch disp_type
                case {'twtt' 'flat'}
                    ind_y_pk= interp1(block.twtt, 1:num_sample_trim, (1e-6 * ind_y_pk), 'nearest', 'extrap');
                case '~depth'
                    ind_y_pk= interp1(block.twtt, 1:num_sample_trim, (1e-6 * (ind_y_pk + (1e6 * (block.twtt_surf(ind_x_pk) - block.twtt(1))))), 'nearest', 'extrap');
            end
            
            if (button == 1) % trace layer
                
                pk.num_layer= pk.num_layer + 1;
                curr_layer  = pk.num_layer;
                [smooth_done, p_pksmooth, p_pksmoothdepth, p_pksmoothflat, ind_y_flat_mean, ind_y_flat_smooth] ...
                            = deal([smooth_done false], [p_pksmooth 0], [p_pksmoothdepth 0], [p_pksmoothflat 0], [ind_y_flat_mean; NaN(1, num_decim_flat)], [ind_y_flat_smooth; NaN(1, num_decim_flat)]);
                pk_prop
                set(status_box, 'string', ['Layer #' num2str(curr_layer) ' picked.'])
                pause(0.5)
                
            elseif strcmpi(char(button), 'D') % delete layer
                
                if ~pk.num_layer
                    continue
                end
                tmp1        = NaN((pk.num_layer - tmp3 + 1), 1);
                for ii = 1:(pk.num_layer - tmp3 + 1)
                    switch disp_type
                        case 'twtt'
                            tmp1(ii) = pk.layer(tmp3 + ii - 1).ind_y(ind_x_pk); % y index at x index pick for each layer
                        case '~depth'
                            tmp1(ii) = pk.layer(tmp3 + ii - 1).ind_y(ind_x_pk) - ind_surf(ind_x_pk) + 1;
                        case 'flat'
                            tmp1(ii) = ind_y_flat_mean((tmp3 + ii - 1), ind_x_pk);
                    end
                end
                if (length(tmp1(~isnan(tmp1))) > 1)
                    tmp1= interp1(tmp1(~isnan(tmp1)), find(~isnan(tmp1)), ind_y_pk, 'nearest', 'extrap') + tmp3 - 1;
                elseif (length(tmp1(~isnan(tmp1))) == 1)
                    tmp1=find(~isnan(tmp1)) + tmp3 - 1;
                else
                    set(status_box, 'string', 'Cannot determine which layer to delete. Pick a more distinct x index.')
                    pause(0.5)
                    continue
                end
                switch disp_type
                    case 'twtt'
                        delete(p_pk(tmp1))
                    case '~depth'
                        delete(p_pkdepth(tmp1))
                    case 'flat'
                        delete(p_pkflat(tmp1))
                end
                tmp2        = setdiff(1:pk.num_layer, tmp1);
                switch disp_type
                    case 'twtt'
                        p_pk= p_pk(tmp2);
                    case '~depth'
                        p_pkdepth ...
                            = p_pkdepth(tmp2);
                    case 'flat'
                        p_pkflat ...
                            = p_pkflat(tmp2);
                end
                [pk.layer, smooth_done, p_pksmooth, p_pksmoothdepth, p_pksmoothflat, pk.num_layer, ind_y_flat_mean, ind_y_flat_smooth] ...
                            = deal(pk.layer(tmp2), smooth_done(tmp2), p_pksmooth(tmp2), p_pksmoothdepth(tmp2), p_pksmoothflat(tmp2), (pk.num_layer - 1), ind_y_flat_mean(tmp2, :), ind_y_flat_smooth(tmp2, :));
                set(status_box, 'string', ['Deleted layer #' num2str(tmp1) '.'])
                pause(0.5)
                
            elseif strcmpi(char(button), 'U') % undo
                
                if (pk.num_layer < tmp3)
                    continue
                end
                switch disp_type
                    case 'twtt'
                        delete(p_pk(end))
                        p_pk= p_pk(1:(end - 1));
                    case '~depth'
                        delete(p_pkdepth(end))
                        p_pkdepth ...
                            = p_pkdepth(1:(end - 1));
                    case 'flat'
                        delete(p_pkflat(end))
                        p_pkflat ...
                            = p_pkflat(1:(end - 1));
                end
                [pk.layer, smooth_done, p_pksmooth, p_pksmoothdepth, p_pksmoothflat, pk.num_layer, ind_y_flat_mean, ind_y_flat_smooth] ...
                            = deal(pk.layer(1:(end - 1)), smooth_done(1:(end - 1)), p_pksmooth(1:(end - 1)), p_pksmoothdepth(1:(end - 1)), p_pksmoothflat(1:(end - 1)), (pk.num_layer - 1), ind_y_flat_mean(1:(end - 1), :), ind_y_flat_smooth(1:(end - 1), :));
                set(status_box, 'string', 'Undid last layer.')
                pause(0.5)
                
            elseif (strcmpi(char(button), 'L') || strcmpi(char(button), 'R') || strcmpi(char(button), 'C'))  % delete portion of layer
                
                if (pk.num_layer < tmp3)
                    continue
                end
                tmp1        = NaN((pk.num_layer - tmp3 + 1), 1);
                try
                    for ii = 1:(pk.num_layer - tmp3 + 1) % y index at x index pick for each layer
                        switch disp_type
                            case {'twtt' '~depth'}
                                tmp1(ii) = pk.layer(tmp3 + ii - 1).ind_y(ind_x_pk);
                            case 'flat'
                                tmp1(ii) = ind_y_flat_mean((tmp3 + ii - 1), ind_x_pk);
                        end
                    end
                catch
                    set(status_box, 'string', 'Something went wrong compiling new layer y-indices. Try again.')
                    pause(0.5)
                    continue
                end
                tmp2        = find(~isnan(tmp1));
                if ((length(tmp2) > 1) && (length(unique(tmp1(tmp2))) == length(tmp2)))
                    tmp1    = interp1(tmp1(tmp2), tmp2, ind_y_pk, 'nearest', 'extrap') + tmp3 - 1;
                elseif (length(tmp2) == 1)
                    tmp1    = tmp3 - 1 + tmp2;
                else
                    set(status_box, 'string', 'Cannot determine which new layer to edit. Pick a more distinct x index.')
                    pause(0.5)
                    continue
                end
                if strcmpi(char(button), 'C')
                    set(status_box, 'string', 'Now choose right end of cut...');
                    [ind_x_pk(2), ~] ...
                            = ginput(1);
                    switch disp_type
                        case {'twtt' '~depth'}
                            ind_x_pk(2) = interp1(block.dist_lin(ind_decim), ind_decim, ind_x_pk(2), 'nearest', 'extrap');
                        case 'flat'
                            ind_x_pk(2) = interp1(block.dist_lin(ind_decim_flat), 1:num_decim_flat, ind_x_pk(2), 'nearest', 'extrap');
                    end
                end
                switch char(button)
                    case {'l' 'L'}
                        tmp2= 1:ind_x_pk;
                    case {'r' 'R'}
                        switch disp_type
                            case {'twtt' '~depth'}
                                tmp2 = ind_x_pk:block.num_trace;
                            case 'flat'
                                tmp2 = ind_x_pk:num_decim_flat;
                        end
                    case {'c' 'C'}
                        tmp2= ind_x_pk(1):ind_x_pk(2);
                end
                switch disp_type
                    case {'twtt' '~depth'}
                        pk.layer(tmp1).ind_y(tmp2) ...
                            = NaN;
                        if strcmp(disp_type, 'twtt')
                            delete(p_pk(tmp1))
                        else
                            delete(p_pkdepth(tmp1))
                        end
                    case 'flat'
                        ind_y_flat_mean(tmp1, tmp2) ...
                            = NaN;
                        delete(p_pkflat(tmp1))
                end
                switch disp_type
                    case {'twtt' '~depth'}
                        if all(isnan(pk.layer(tmp1).ind_y))
                            tmp2 = true;
                        else
                            tmp2 = false;
                        end
                    case 'flat'
                        if all(isnan(ind_y_flat_mean(tmp1, :)))
                            tmp2 = true;
                        else
                            tmp2 = false;
                        end
                end
                if tmp2
                    tmp4    = setdiff(1:pk.num_layer, tmp1);
                    switch disp_type
                        case 'twtt'
                            p_pk = p_pk(tmp4);
                        case '~depth'
                            p_pkdepth = p_pkdepth(tmp4);
                        case 'flat'
                            p_pkflat = p_pkflat(tmp4);
                    end
                    [pk.layer, smooth_done, p_pksmooth, p_pksmoothdepth, p_pksmoothflat, pk.num_layer, ind_y_flat_mean, ind_y_flat_smooth] ...
                            = deal(pk.layer(tmp4), smooth_done(tmp4), p_pksmooth(tmp4), p_pksmoothdepth(tmp4), p_pksmoothflat(tmp4), (pk.num_layer - 1), ind_y_flat_mean(tmp4, :), ind_y_flat_smooth(tmp4, :));
                    set(status_box, 'string', 'Deleted edited layer because it is now empty.')
                    pause(0.5)
                else
                    switch disp_type
                        case 'twtt'
                            p_pk(tmp1) ...
                                = plot(block.dist_lin(ind_decim(~isnan(pk.layer(tmp1).ind_y(ind_decim)))), (1e6 .* block.twtt(pk.layer(tmp1).ind_y(ind_decim(~isnan(pk.layer(tmp1).ind_y(ind_decim)))))), '.', 'color', [1 0.7 0.7], 'markersize', 12);
                        case '~depth'
                            tmp2 = ind_decim(~isnan(pk.layer(tmp1).ind_y(ind_decim)) & ~isnan(ind_surf(ind_decim)));
                            tmp4 = pk.layer(tmp1).ind_y(tmp2) - ind_surf(tmp2) + 1;
                            tmp4((tmp4 < 1) | (tmp4 > num_sample_trim)) ...
                                 = NaN;
                            p_pkdepth(tmp1) ...
                                = plot(block.dist_lin(tmp2(~isnan(tmp4))), (1e6 .* block.twtt(round(tmp4(~isnan(tmp4))))), '.', 'color', [1 0.7 0.7], 'markersize', 12);
                        case 'flat'
                            p_pkflat(tmp1) ...
                                = plot(block.dist_lin(ind_decim_flat(~isnan(ind_y_flat_mean(tmp1, :)))), (1e6 .* block.twtt(ind_y_flat_mean(tmp1, ~isnan(ind_y_flat_mean(tmp1, :))))), '.', 'color', [1 0.7 0.7], 'markersize', 12);
                    end
                    switch char(button)
                        case {'l' 'L'}
                            set(status_box, 'string', ['Trimmed left of layer #' num2str(tmp1) ' at ' num2str(ind_x_pk, '%3.1f') ' km.'])
                        case {'r' 'R'}
                            set(status_box, 'string', ['Trimmed right of layer #' num2str(tmp1) ' at ' num2str(ind_x_pk, '%3.1f') ' km.'])
                        case {'c' 'C'}
                            set(status_box, 'string', ['Trimmed layer #' num2str(tmp1) ' between ' num2str(ind_x_pk(1), '%3.1f') ' and ' num2str(ind_x_pk(2), '%3.1f') ' km.'])
                    end
                end
                
            elseif strcmpi(char(button), 'M') % merge two new layers
                
                if (pk.num_layer < (tmp3 + 1))
                    set(status_box, 'string', 'Not enough new layers to merge.')
                    continue
                end
                tmp1        = NaN((pk.num_layer - tmp3 + 1), 1);
                for ii = 1:(pk.num_layer - tmp3 + 1)
                    switch disp_type
                        case {'twtt' '~depth'}
                            tmp1(ii) = pk.layer(tmp3 + ii - 1).ind_y(ind_x_pk);
                        case 'flat'
                            tmp1(ii) = ind_y_flat_mean((tmp3 + ii - 1), ind_x_pk);
                    end
                end
                tmp2        = find(~isnan(tmp1));
                if ((length(tmp2) > 1) && (length(unique(tmp1(tmp2))) == length(tmp2)))
                    tmp1    = interp1(tmp1(tmp2), tmp2, ind_y_pk, 'nearest', 'extrap') + tmp3 - 1;
                elseif (length(tmp2) == 1)
                    tmp1    = tmp3 - 1 + tmp2;
                else
                    set(status_box, 'string', 'Cannot determine which new layer to edit. Pick a more distinct x index.')
                    pause(0.5)
                    continue
                end
                
                switch disp_type
                    case 'twtt'
                        set(p_pk(tmp1), 'color', 'y')
                    case '~depth'
                        set(p_pkdepth(tmp1), 'color', 'y')
                    case 'flat'
                        set(p_pkflat(tmp1), 'color', 'y')
                end
                
                set(status_box, 'string', 'Now choose layer to merge with...')
                [ind_x_pk, ind_y_pk] ...
                            = ginput(1);
                switch disp_type
                    case {'twtt' '~depth'}
                        ind_x_pk = interp1(block.dist_lin(ind_decim), ind_decim, ind_x_pk, 'nearest', 'extrap');
                    case 'flat'
                        ind_x_pk = interp1(block.dist_lin(ind_decim_flat), 1:num_decim_flat, ind_x_pk, 'nearest', 'extrap');
                end
                switch disp_type
                    case {'twtt' 'flat'}
                        ind_y_pk ...
                            = interp1(block.twtt, 1:num_sample_trim, (1e-6 * ind_y_pk), 'nearest', 'extrap');
                    case '~depth'
                        ind_y_pk ...
                            = interp1(block.twtt, 1:num_sample_trim, (1e-6 * (ind_y_pk + (1e6 * (block.twtt_surf(ind_x_pk) - block.twtt(1))))), 'nearest', 'extrap');
                end
                tmp2        = NaN((pk.num_layer - tmp3 + 1), 1);
                for ii = 1:(pk.num_layer - tmp3 + 1)
                    switch disp_type
                        case {'twtt' '~depth'}
                            tmp2(ii) ...
                                = pk.layer(tmp3 + ii - 1).ind_y(ind_x_pk);
                        case 'flat'
                            tmp2(ii) ...
                                = ind_y_flat_mean((tmp3 + ii - 1), ind_x_pk);
                    end
                end
                tmp4        = find(~isnan(tmp2));
                if ((length(tmp4) > 1) && (length(unique(tmp2(tmp4))) == length(tmp4)))
                    tmp2    = interp1(tmp2(tmp4), tmp4, ind_y_pk, 'nearest', 'extrap') + tmp3 - 1;
                elseif (length(tmp4) == 1)
                    tmp2    = tmp3 - 1 + tmp4;
                else
                    set(status_box, 'string', 'Cannot determine which layer to merge with. Pick a more distinct x index.')
                    switch disp_type
                        case 'twtt'
                            set(p_pk(tmp1), 'color', [1 0.7 0.7])
                        case '~depth'
                            set(p_pkdepth(tmp1), 'color', [1 0.7 0.7])
                        case 'flat'
                            set(p_pkflat(tmp1), 'color', [1 0.7 0.7])
                    end
                    pause(0.5)
                    continue
                end
                
                switch disp_type
                    case 'twtt'
                        set(p_pk(tmp2), 'color', 'y')
                    case '~depth'
                        set(p_pkdepth(tmp2), 'color', 'y')
                    case 'flat'
                        set(p_pkflat(tmp2), 'color', 'y')
                end
                pause(0.25)
                
                tmp4        = setdiff(1:pk.num_layer, tmp2);
                
                switch disp_type
                    case {'twtt' '~depth'}
                        pk.layer(tmp1).ind_y(isnan(pk.layer(tmp1).ind_y)) ...
                            = pk.layer(tmp2).ind_y(isnan(pk.layer(tmp1).ind_y));
                        if strcmp(disp_type, 'twtt')
                            delete(p_pk([tmp1 tmp2]))
                            p_pk(tmp1) ...
                                = plot(block.dist_lin(ind_decim(~isnan(pk.layer(tmp1).ind_y(ind_decim)))), (1e6 .* block.twtt(pk.layer(tmp1).ind_y(ind_decim(~isnan(pk.layer(tmp1).ind_y(ind_decim)))))), '.', 'color', [1 0.7 0.7], 'markersize', 12);
                            p_pk= p_pk(tmp4);
                        else
                            delete(p_pkdepth([tmp1 tmp2]))
                            tmp2= ind_decim(~isnan(pk.layer(tmp1).ind_y(ind_decim)) & ~isnan(ind_surf(ind_decim)));
                            tmp5= pk.layer(tmp1).ind_y(tmp2) - ind_surf(tmp2) + 1;
                            tmp5((tmp5 < 1) | (tmp5 > num_sample_trim)) ...
                                 = NaN;
                            p_pkdepth(tmp1) ...
                                = plot(block.dist_lin(tmp2(~isnan(tmp5))), (1e6 .* block.twtt(round(tmp5(~isnan(tmp5))))), '.', 'color', [1 0.7 0.7], 'markersize', 12);
                            p_pkdepth ...
                                = p_pkdepth(tmp4);
                        end
                    case 'flat'
                        ind_y_flat_mean(tmp1, isnan(ind_y_flat_mean(tmp1, :))) ...
                            = ind_y_flat_mean(tmp2, isnan(ind_y_flat_mean(tmp1, :)));
                        delete(p_pkflat([tmp1 tmp2]))
                        p_pkflat(tmp1) ...
                            = plot(block.dist_lin(ind_decim_flat(~isnan(ind_y_flat_mean(tmp1, :)))), (1e6 .* block.twtt(ind_y_flat_mean(tmp1, ~isnan(ind_y_flat_mean(tmp1, :))))), '.', 'color', [1 0.7 0.7], 'markersize', 12);
                        p_pkflat ...
                            = p_pkflat(tmp4);
                end
                
                [pk.layer, smooth_done, p_pksmooth, p_pksmoothdepth, p_pksmoothflat, pk.num_layer, ind_y_flat_mean, ind_y_flat_smooth] ...
                            = deal(pk.layer(tmp4), smooth_done(tmp4), p_pksmooth(tmp4), p_pksmoothdepth(tmp4), p_pksmoothflat(tmp4), (pk.num_layer - 1), ind_y_flat_mean(tmp4, :), ind_y_flat_smooth(tmp4, :));
                
                set(status_box, 'string', 'Layers merged.')
                pause(0.5)
                
            elseif strcmpi(char(button), 'E')
                
                reset_xy
                
            elseif strcmpi(char(button), 'W')
                
                pk.num_win  = pk.num_win + 1;
                set(num_win_edit, 'string', num2str(pk.num_win))
                set(status_box, 'string', ['Vertical search window adjusted to +/- ' num2str(pk.num_win) ' samples.'])
                pause(0.5)
            
            elseif strcmpi(char(button), 'S')
                
                if (pk.num_win > 1)
                    pk.num_win ...
                            = pk.num_win - 1;
                    set(num_win_edit, 'string', num2str(pk.num_win))
                    set(status_box, 'string', ['Vertical search window adjusted to +/- ' num2str(pk.num_win) ' samples.'])
                    pause(0.5)
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
                
            elseif strcmpi(char(button), 'Q') % done picking lines
                
                set(status_box, 'string', 'Done picking layers...')
                break
                
            end
        end
        
        set(pkgui, 'keypressfcn', @keypress, 'windowbuttondownfcn', @mouse_click)
        
        if ~pk.num_layer
            set(status_box, 'string', 'No layers picked/left.')
            return
        end
        
        switch disp_type
            case 'twtt'
                set(p_pk(logical(p_pk) & ishandle(p_pk)), 'color', 'r')
            case '~depth'
                set(p_pkdepth(logical(p_pkdepth) & ishandle(p_pkdepth)), 'color', 'r')
            case 'flat'
                set(p_pkflat(logical(p_pkflat) & ishandle(p_pkflat)), 'color', 'r')
        end
        
        % plot picks in other spaces
        switch disp_type
            
            case {'twtt' '~depth'}
                
                if strcmp(disp_type, 'twtt')
                    
                    if (any(p_pkdepth) && any(ishandle(p_pkdepth)))
                        delete(p_pkdepth(logical(p_pkdepth) & ishandle(p_pkdepth)))
                    end
                    p_pkdepth ...
                            = zeros(1, pk.num_layer);
                    if depth_avail
                        for ii = 1:pk.num_layer
                            tmp1= ind_decim(~isnan(pk.layer(ii).ind_y(ind_decim)) & ~isnan(ind_surf(ind_decim)));
                            tmp2= pk.layer(ii).ind_y(tmp1) - ind_surf(tmp1) + 1;
                            tmp2((tmp2 < 1) | (tmp2 > num_sample_trim)) ...
                                = NaN;
                            p_pkdepth(ii) ...
                                = plot(block.dist_lin(tmp1(~isnan(tmp2))), (1e6 .* block.twtt(round(tmp2(~isnan(tmp2))))), 'r.', 'markersize', 12, 'visible', 'off');
                        end
                    end
                    
                else
                    if (any(p_pk) && any(ishandle(p_pk)))
                        delete(p_pk(logical(p_pk) & ishandle(p_pk)))
                    end
                    p_pk    = zeros(1, pk.num_layer);
                    for ii = 1:pk.num_layer
                        tmp1= pk.layer(ii).ind_y(ind_decim);
                        p_pk(ii) ...
                            = plot(block.dist_lin(ind_decim(~isnan(tmp1))), (1e6 .* block.twtt(round(tmp1(~isnan(tmp1))))), 'r.', 'markersize', 12, 'visible', 'off');
                    end
                    
                end
                
                if (any(p_pkflat) && any(ishandle(p_pkflat)))
                    delete(p_pkflat(logical(p_pkflat) & ishandle(p_pkflat)))
                end
                p_pkflat    = zeros(1, pk.num_layer);
                if flat_done
                    warning('off', 'MATLAB:interp1:NaNinY')
                    ind_y_curr ...
                            = NaN(pk.num_layer, block.num_trace);
                    for ii = 1:pk.num_layer
                        ind_y_curr(ii, :) ...
                            = pk.layer(ii).ind_y;
                    end
                    for ii = find(~sum(isnan(ind_y_flat(:, ind_decim_flat))))
                        [~, tmp4] ...
                            = unique(ind_y_flat(:, ind_decim_flat(ii)));
                        ind_y_flat_mean(~isnan(ind_y_curr(:, ind_decim_flat(ii))), ii) ...
                            = interp1(ind_y_flat(tmp4, ind_decim_flat(ii)), tmp4, ind_y_curr(~isnan(ind_y_curr(:, ind_decim_flat(ii))), ind_decim_flat(ii)), 'nearest', 'extrap');
                    end
                    warning('on', 'MATLAB:interp1:NaNinY')
                    for ii = 1:pk.num_layer
                        p_pkflat(ii) ...
                            = plot(block.dist_lin(ind_decim_flat(~isnan(ind_y_flat_mean(ii, :)))), (1e6 .* block.twtt(ind_y_flat_mean(ii, ~isnan(ind_y_flat_mean(ii, :))))), 'r.', 'markersize', 12, 'visible', 'off');
                    end
                end
                
            case 'flat'
                
                warning('off', 'MATLAB:interp1:NaNinY')
                for ii = 1:pk.num_layer
                    pk.layer(ii).ind_y ...
                            = NaN(1, block.num_trace);
                    tmp1    = round(interp1(ind_decim_flat, ind_y_flat_mean(ii, :), 1:block.num_trace, 'linear', 'extrap'));
                    tmp1((tmp1 < 1) | (tmp1 > num_sample_trim)) ...
                            = NaN;
                    tmp2    = 1:block.num_trace;
                    tmp2    = tmp2(~isnan(tmp1));
                    try %#ok<TRYNC>
                        pk.layer(ii).ind_y(~isnan(tmp1)) ...
                            = ind_y_flat(sub2ind([num_sample_trim block.num_trace], tmp1(~isnan(tmp1)), tmp2));
                    end
                    pk.layer(ii).ind_y((pk.layer(ii).ind_y < 1) | (pk.layer(ii).ind_y > num_sample_trim)) ...
                            = NaN;
                    for jj = 1:(num_decim_flat - 1)
                        if all(isnan(ind_y_flat_mean(ii, jj:(jj + 1))))
                            pk.layer(ii).ind_y(ind_decim_flat(jj):ind_decim_flat(jj + 1)) ...
                            = NaN; % deal with breaks in layers
                        end
                    end
                end
                warning('on', 'MATLAB:interp1:NaNinY')
                
                % plot picked layers in twtt space
                if (any(p_pk) && any(ishandle(p_pk)))
                    delete(p_pk(logical(p_pk) & ishandle(p_pk)))
                end
                if (any(p_pkdepth) && any(ishandle(p_pkdepth)))
                    delete(p_pkdepth(logical(p_pkdepth) & ishandle(p_pkdepth)))
                end
                
                [p_pk, p_pkdepth] ...
                            = deal(zeros(1, pk.num_layer));
                for ii = 1:pk.num_layer
                    tmp1    = pk.layer(ii).ind_y(ind_decim);
                    tmp2    = block.dist_lin(ind_decim);
                    p_pk(ii)= plot(tmp2(~isnan(tmp1)), (1e6 .* block.twtt(round(tmp1(~isnan(tmp1))))), 'r.', 'markersize', 12, 'visible', 'off');
                end
                if depth_avail
                    for ii = 1:pk.num_layer
                        tmp1= ind_decim(~isnan(pk.layer(ii).ind_y(ind_decim)) & ~isnan(ind_surf(ind_decim)));
                        tmp2= pk.layer(ii).ind_y(tmp1) - ind_surf(tmp1) + 1;
                        tmp2((tmp2 < 1) | (tmp2 > num_sample_trim)) ...
                            = NaN;
                        p_pkdepth(ii) ...
                            = plot(block.dist_lin(tmp1(~isnan(tmp2))), (1e6 .* block.twtt(round(tmp2(~isnan(tmp2))))), 'r.', 'markersize', 12, 'visible', 'off');
                    end
                end
        end
        
        set(pk_check, 'value', 1)
        
        pk_sort % sort layers
        
        if (tmp3 <= pk.num_layer)
            curr_layer      = find((pk.num_layer == tmp1), 1); % highlight last picked layer, here tmp1 is the sorted indices from pk_sort (a bit obtuse...sorry)
        elseif curr_layer
            curr_layer      = find((curr_layer == tmp1), 1); % highlight same current layer, roughly
        else
            curr_layer      = pk.num_layer; % just starting out
        end
        if (isempty(curr_layer) || any(curr_layer < 1) || any(curr_layer > pk.num_layer))
            curr_layer      = 1;
        end
        
        set(layer_list, 'string', [num2cell(1:pk.num_layer) 'surface' 'bed'], 'value', curr_layer)
        
        pk_done             = true;
        match_done          = false;
        pk_smooth
        pk_select
    end

%% Propagate layer from pick

    function pk_prop(source, eventdata)
        
        pk.layer(curr_layer).ind_y ...
                            = NaN(1, block.num_trace);
        switch disp_type
            case 'twtt'
                [~, pk.layer(curr_layer).ind_y(ind_x_pk)] ...
                            = max(block.amp((ind_y_pk - pk.num_win):(ind_y_pk + pk.num_win), ind_x_pk)); % y index of nearest min/max
                pk.layer(curr_layer).ind_y(ind_x_pk) ...
                            = ind_y_pk - ((pk.num_win + 1) - pk.layer(curr_layer).ind_y(ind_x_pk)); % correct y index of min/max because search was done in a narrow window
            case '~depth'
                tmp1        = NaN(1, num_decim);
                tmp2        = [interp1(ind_decim, 1:num_decim, ind_x_pk, 'nearest', 'extrap') (ind_y_pk - ind_surf(ind_x_pk) + 1)];
                [~, tmp1(tmp2(1))] ...
                            = max(amp_depth((tmp2(2) - pk.num_win):(tmp2(2) + pk.num_win), tmp2(1)));
                tmp1(tmp2(1)) ...
                            = tmp2(2) - ((pk.num_win + 1) - tmp1(tmp2(1)));
            case 'flat'
                try %#ok<TRYNC>
                    [~, ind_y_flat_mean(curr_layer, ind_x_pk)] ...
                            = max(amp_flat_mean((ind_y_pk - pk.num_win):(ind_y_pk + pk.num_win), ind_x_pk));
                    ind_y_flat_mean(curr_layer, ind_x_pk) ...
                            = ind_y_pk - ((pk.num_win + 1) - ind_y_flat_mean(curr_layer, ind_x_pk));
                end
        end
        
        % loop for left of ind_x_pk
        switch disp_type
            case {'twtt' 'flat'}
                for ii = (ind_x_pk - 1):-1:1
                    try
                        switch disp_type
                            case 'twtt'
                                [~, pk.layer(curr_layer).ind_y(ii)] ...
                                    = max(block.amp((pk.layer(curr_layer).ind_y(ii + 1) - pk.num_win):(pk.layer(curr_layer).ind_y(ii + 1) + pk.num_win), ii));
                                pk.layer(curr_layer).ind_y(ii) ...
                                    = pk.layer(curr_layer).ind_y(ii + 1) - ((pk.num_win + 1) - pk.layer(curr_layer).ind_y(ii));
                            case 'flat'
                                [~, ind_y_flat_mean(curr_layer, ii)] ...
                                    = max(amp_flat_mean((ind_y_flat_mean(curr_layer, (ii + 1)) - pk.num_win):(ind_y_flat_mean(curr_layer, (ii + 1)) + pk.num_win), ii));
                                ind_y_flat_mean(curr_layer, ii) ...
                                    = ind_y_flat_mean(curr_layer, (ii + 1)) - ((pk.num_win + 1) - ind_y_flat_mean(curr_layer, ii));
                        end
                    catch %#ok<*CTCH>
                        break
                    end
                end
            case '~depth'
                for ii = (tmp2(1) - 1):-1:1
                    try
                        [~, tmp1(ii)] ...
                            = max(amp_depth((tmp1(ii + 1) - pk.num_win):(tmp1(ii + 1) + pk.num_win), ii));
                        tmp1(ii) ...
                            = tmp1(ii + 1) - ((pk.num_win + 1) - tmp1(ii));
                    catch
                        break
                    end
                end
        end
        
        switch disp_type
            case {'twtt' '~depth'}
                tmp4        = block.num_trace;
            case 'flat'
                tmp4        = num_decim_flat;
        end
        
        % loop for right of ind_x_pk
        switch disp_type
            case {'twtt' 'flat'}
                for ii = (ind_x_pk + 1):tmp4
                    try
                        switch disp_type
                            case 'twtt'
                                [~, pk.layer(curr_layer).ind_y(ii)] ...
                                    = max(block.amp((pk.layer(curr_layer).ind_y(ii - 1) - pk.num_win):(pk.layer(curr_layer).ind_y(ii - 1) + pk.num_win), ii));
                                pk.layer(curr_layer).ind_y(ii) ...
                                    = pk.layer(curr_layer).ind_y(ii - 1) - ((pk.num_win + 1) - pk.layer(curr_layer).ind_y(ii));
                            case 'flat'
                                [~, ind_y_flat_mean(curr_layer, ii)] ...
                                    = max(amp_flat_mean((ind_y_flat_mean(curr_layer, (ii - 1)) - pk.num_win):(ind_y_flat_mean(curr_layer, (ii - 1)) + pk.num_win), ii));
                                ind_y_flat_mean(curr_layer, ii) ...
                                    = ind_y_flat_mean(curr_layer, (ii - 1)) - ((pk.num_win + 1) - ind_y_flat_mean(curr_layer, ii));
                        end
                    catch %#ok<CTCH>
                        break
                    end
                end
            case '~depth'
                for ii = (tmp2(1) + 1):num_decim
                    try
                        [~, tmp1(ii)] ...
                            = max(amp_depth((tmp1(ii - 1) - pk.num_win):(tmp1(ii - 1) + pk.num_win), ii));
                        tmp1(ii) ...
                            = tmp1(ii - 1) - ((pk.num_win + 1) - tmp1(ii));
                    catch
                        break
                    end
                end
        end
        
        if strcmp(disp_type, '~depth')
            pk.layer(curr_layer).ind_y ...
                            = interp1(ind_decim, (tmp1 + ind_surf(ind_decim) - 1), 1:block.num_trace, 'linear', 'extrap');
        end
        
        % plot new layer
        switch disp_type
            case 'twtt'
                p_pk(curr_layer) ...
                            = plot(block.dist_lin(ind_decim(~isnan(pk.layer(curr_layer).ind_y(ind_decim)))), (1e6 .* block.twtt(pk.layer(curr_layer).ind_y(ind_decim(~isnan(pk.layer(curr_layer).ind_y(ind_decim)))))), '.', 'color', [1 0.7 0.7], 'markersize', 12);
            case '~depth'
                tmp1        = ind_decim(~isnan(pk.layer(curr_layer).ind_y(ind_decim)) & ~isnan(ind_surf(ind_decim)));
                tmp2        = pk.layer(curr_layer).ind_y(tmp1) - ind_surf(tmp1) + 1;
                tmp2((tmp2 < 1) | (tmp2 > num_sample_trim)) ...
                            = NaN;
                p_pkdepth(curr_layer) ...
                            = plot(block.dist_lin(tmp1(~isnan(tmp2))), (1e6 .* block.twtt(tmp2(~isnan(tmp2)))), '.', 'color', [1 0.7 0.7], 'markersize', 12);
            case 'flat'
                p_pkflat(curr_layer) ...
                            = plot(block.dist_lin(ind_decim_flat(~isnan(ind_y_flat_mean(curr_layer, :)))), (1e6 .* block.twtt(ind_y_flat_mean(curr_layer, ~isnan(ind_y_flat_mean(curr_layer, :))))), '.', 'color', [1 0.7 0.7], 'markersize', 12);
        end
    end

%% Manually pick a layer

    function pk_man(source, eventdata)
        
        if ~any(strcmp(disp_type, {'twtt' '~depth'}))
            set(status_box, 'string', 'Layers can only be traced manually in twtt or ~depth.')
            return
        end
        
        axes(ax_radar)
        ii                  = 0;
        [ind_x_pk, ind_y_pk, button] ...
                            = deal([]);
        
        set(pkgui, 'keypressfcn', [], 'windowbuttondownfcn', [])
        
        while true
            
            set(status_box, 'string', 'Left-click: pick; U: undo; return: done...')
            [tmp1, tmp2, button] ...
                            = ginput(1); % pick manually, press enter when done
            if ~ii % dummy initialization
                ii          = 1;
                continue
            end
            
            if (button == 1)
                
                [ind_x_pk(ii), ind_y_pk(ii)] ...
                            = deal(tmp1, tmp2);
                ind_x_pk(ii)= interp1(block.dist_lin(ind_decim), ind_decim, ind_x_pk(ii), 'nearest', 'extrap'); % raw picks must be indices, not dimensionalized vectors (horizontal)
                switch disp_type
                    case 'twtt'
                        ind_y_pk(ii) ...
                            = interp1(block.twtt, 1:num_sample_trim, (1e-6 .* ind_y_pk(ii)), 'nearest', 'extrap'); % interpolate traveltime pick onto traveltime vector
                    case '~depth'
                        ind_y_pk(ii) ...
                            = interp1(block.twtt, 1:num_sample_trim, (1e-6 .* (ind_y_pk(ii) + (1e6 * (block.twtt_surf(ind_x_pk(ii)) - block.twtt(1))))), 'nearest', 'extrap');
                end
                if (ii > 1)
                    delete(tmp3) % get rid of old plot handle
                end
                switch disp_type
                    case 'twtt'
                        tmp3= plot(block.dist_lin(ind_x_pk), (1e6 .* block.twtt(ind_y_pk)), 'rx', 'markersize', 12); % original picks
                    case '~depth'
                        tmp3= plot(block.dist_lin(ind_x_pk), (1e6 .* block.twtt(ind_y_pk - ind_surf(ind_x_pk) + 1)), 'rx', 'markersize', 12); % original picks
                end
                ii          = ii + 1;
                
            elseif (strcmpi(char(button), 'U') && (ii > 1))
                
                ind_x_pk    = ind_x_pk(1:(end - 1));
                ind_y_pk    = ind_y_pk(1:(end - 1));
                if (logical(tmp3) && ishandle(tmp3))
                    delete(tmp3)
                end
                if (ii > 2)
                    switch disp_type
                        case 'twtt'
                            tmp3 ...
                                = plot(block.dist_lin(ind_x_pk), (1e6 .* block.twtt(ind_y_pk)), 'rx', 'markersize', 12); % original picks
                        case '~depth'
                            tmp3 ...
                                = plot(block.dist_lin(ind_x_pk), (1e6 .* block.twtt(ind_y_pk - ind_surf(ind_x_pk) + 1)), 'rx', 'markersize', 12); % original picks
                    end
                end
                ii          = ii - 1;
                
            elseif (double(get(pkgui, 'currentcharacter')) == 13)
                
                if (ii < 4)
                    if (logical(tmp3) && ishandle(tmp3))
                        delete(tmp3)
                    end
                    set(status_box, 'string', 'Not enough picked points to make a manual layer. Start over.')
                    set(pkgui, 'keypressfcn', @keypress, 'windowbuttondownfcn', @mouse_click)
                    break
                end
                
                if ~issorted(ind_x_pk) % resort picks
                    [ind_x_pk, tmp1] ...
                            = sort(ind_x_pk);
                    ind_y_pk= ind_y_pk(tmp1);
                end
                if (length(unique(ind_x_pk)) < length(ind_x_pk)) % don't keep picks that are accidentally at the same horizontal index, otherwise spline will fail
                    [ind_x_pk, tmp1] ...
                            = unique(ind_x_pk);
                    ind_y_pk= ind_y_pk(tmp1);
                end
                
                pk.num_layer= pk.num_layer + 1;
                
                pk.layer(pk.num_layer).ind_y ...
                            = NaN(1, block.num_trace);
                pk.layer(pk.num_layer).ind_y(ind_x_pk(1):decim:ind_x_pk(end)) ...
                            = round(spline(ind_x_pk, ind_y_pk, ind_x_pk(1):decim:ind_x_pk(end))); % interpolate spline through picks
                for ii = interp1(ind_decim, 1:num_decim, ind_x_pk, 'nearest', 'extrap') % check for local maxima
                    [~, tmp2] ...
                            = max(amp_mean((pk.layer(pk.num_layer).ind_y(ind_decim(ii)) - pk.num_win):(pk.layer(pk.num_layer).ind_y(ind_decim(ii)) + pk.num_win), ii));
                    pk.layer(pk.num_layer).ind_y(ind_decim(ii)) ...
                            = pk.layer(pk.num_layer).ind_y(ind_decim(ii)) - pk.num_win - 1 + tmp2; % adjust local maximum index
                end
                pk.layer(pk.num_layer).ind_y(ind_x_pk(1):ind_x_pk(end)) ...
                            = round(spline(ind_x_pk, pk.layer(pk.num_layer).ind_y(ind_x_pk), ind_x_pk(1):ind_x_pk(end))); % interpolate spline through improved picks
                
                tmp1        = find(~isnan(pk.layer(pk.num_layer).ind_y(ind_decim)));
                delete(tmp3)
                p_pk(pk.num_layer) ...
                            = plot(block.dist_lin(ind_decim(tmp1)), (1e6 .* block.twtt(round(pk.layer(pk.num_layer).ind_y(ind_decim(tmp1))))), 'r.', 'markersize', 12, 'visible', 'off');
                
                if depth_avail
                    tmp1    = ind_decim(~isnan(pk.layer(pk.num_layer).ind_y(ind_decim)) & ~isnan(ind_surf(ind_decim)));
                    tmp2    = pk.layer(pk.num_layer).ind_y(tmp1) - ind_surf(tmp1) + 1;
                    tmp2((tmp2 < 1) | (tmp2 > num_sample_trim)) ...
                            = NaN;
                    p_pkdepth(pk.num_layer) ...
                            = plot(block.dist_lin(tmp1(~isnan(tmp2))), (1e6 .* block.twtt(round(tmp2(~isnan(tmp2))))), 'r.', 'markersize', 12, 'visible', 'off');
                else
                    p_pkdepth(pk.num_layer) ...
                            = 0;
                end
                
                [ind_y_flat_mean, ind_y_flat_smooth] ...
                            = deal([ind_y_flat_mean; NaN(1, num_decim_flat)], [ind_y_flat_smooth; NaN(1, num_decim_flat)]);
                if flat_done
                    warning('off', 'MATLAB:interp1:NaNinY')
                    for ii = interp1(ind_decim_flat, 1:num_decim_flat, ind_x_pk(1), 'nearest', 'extrap'):interp1(ind_decim_flat, 1:num_decim_flat, ind_x_pk(end), 'nearest', 'extrap')
                        if find(~sum(isnan(ind_y_flat(:, ind_decim_flat(ii)))))
                            [~, tmp1] ...
                                = unique(ind_y_flat(:, ind_decim_flat(ii)));
                            ind_y_flat_mean(pk.num_layer, ii) ...
                                = interp1(ind_y_flat(tmp1, ind_decim_flat(ii)), tmp1, pk.layer(pk.num_layer).ind_y(ind_decim_flat(ii)), 'nearest', 'extrap');
                        end
                    end
                    ind_y_flat_smooth ...
                            = [ind_y_flat_smooth; NaN(1, num_decim_flat)]; %#ok<AGROW>
                    warning('on', 'MATLAB:interp1:NaNinY')
                    
                    p_pkflat(pk.num_layer) ...
                            = plot(block.dist_lin(ind_decim_flat(~isnan(ind_y_flat_mean(pk.num_layer, :)))), (1e6 .* block.twtt(ind_y_flat_mean(pk.num_layer, ~isnan(ind_y_flat_mean(pk.num_layer, :))))), 'r.', 'markersize', 12, 'visible', 'off');
                else
                    p_pkflat(pk.num_layer) ...
                            = 0;
                end
                
                [smooth_done, p_pksmooth, p_pksmoothdepth, p_pksmoothflat] ...
                            = deal([smooth_done false], [p_pksmooth 0], [p_pksmoothdepth 0], [p_pksmoothflat 0]);
                
                pk_done     = true;
                pk_sort
                curr_layer  = find(tmp1 == pk.num_layer); % tmp1 from pk_sort
                pk_smooth
                set([pk_check smooth_check], 'value', 1)
                set(layer_list, 'string', [num2cell(1:pk.num_layer) 'surface' 'bed'], 'value', curr_layer)
                pk_select
                show_pk
                set(status_box, 'string', 'Manual layer successfully picked.')
                set(pkgui, 'keypressfcn', @keypress, 'windowbuttondownfcn', @mouse_click)
                return
            end
        end
        set(pkgui, 'keypressfcn', @keypress, 'windowbuttondownfcn', @mouse_click)
    end

%% Sort layers from top to bottom based on their mean vertical index

    function pk_sort(source, eventdata)
        tmp1                = zeros(pk.num_layer, 1);
        for ii = 1:pk.num_layer
            tmp1(ii)        = nanmean(pk.layer(ii).ind_y);
        end
        [~, tmp1]           = sort(tmp1);
        [pk.layer, smooth_done, p_pk, p_pkdepth, p_pkflat, p_pksmooth, p_pksmoothdepth, p_pksmoothflat, ind_y_flat_mean, ind_y_flat_smooth] ...
                            = deal(pk.layer(tmp1), smooth_done(tmp1), p_pk(tmp1), p_pkdepth(tmp1), p_pkflat(tmp1), p_pksmooth(tmp1), p_pksmoothdepth(tmp1), p_pksmoothflat(tmp1), ind_y_flat_mean(tmp1, :), ind_y_flat_smooth(tmp1, :));
    end

%% Choose/highlight the current layer

    function pk_select(source, eventdata)
        
        if ~(pk.num_layer || surf_avail || bed_avail)
            set(status_box, 'string', 'No picked layers to select.')
            return
        end
        
        curr_layer          = get(layer_list, 'value');
        
        if (logical(p_bed) && ishandle(p_bed))
            set(p_bed, 'markersize', 12)
        end
        if (logical(p_beddepth) && ishandle(p_beddepth))
            set(p_beddepth, 'markersize', 12)
        end
        if (logical(p_bedflat) && ishandle(p_bedflat))
            set(p_bedflat, 'markersize', 12)
        end
        if (logical(p_surf) && ishandle(p_surf))
            set(p_surf, 'markersize', 12)
        end
        if (logical(p_surfflat) && ishandle(p_surfflat))
            set(p_surfflat, 'markersize', 12)
        end
        if (any(p_pk) && any(ishandle(p_pk)))
            set(p_pk(logical(p_pk) & ishandle(p_pk)), 'markersize', 12)
        end
        if (any(p_pkdepth) && any(ishandle(p_pkdepth)))
            set(p_pkdepth(logical(p_pkdepth) & ishandle(p_pkdepth)), 'markersize', 12)
        end
        if (any(p_pkflat) && any(ishandle(p_pkflat)))
            set(p_pkflat(logical(p_pkflat) & ishandle(p_pkflat)), 'markersize', 12)
        end
        if (any(p_pksmooth) && any(ishandle(p_pksmooth)))
            set(p_pksmooth(logical(p_pksmooth) & ishandle(p_pksmooth)), 'markersize', 12)
        end
        if (any(p_pksmoothdepth) && any(ishandle(p_pksmoothdepth)))
            set(p_pksmoothdepth(logical(p_pksmoothdepth) & ishandle(p_pksmoothdepth)), 'markersize', 12)
        end
        if (any(p_pksmoothflat) && any(ishandle(p_pksmoothflat)))
            set(p_pksmoothflat(logical(p_pksmoothflat) & ishandle(p_pksmoothflat)), 'markersize', 12)
        end
        
        if (curr_layer == (pk.num_layer + 1))
        
            if (logical(p_surf) && ishandle(p_surf))
                set(p_surf, 'markersize', 24)
            end
            if (logical(p_surfflat) && ishandle(p_surfflat))
                set(p_surfflat, 'markersize', 24)
            end
            set(status_box, 'string', 'Surface selected.')
            
        elseif (curr_layer == (pk.num_layer + 2))
            
            if (logical(p_bed) && ishandle(p_bed))
                set(p_bed, 'markersize', 24)
            end
            if (logical(p_beddepth) && ishandle(p_beddepth))
                set(p_beddepth, 'markersize', 24)
            end
            if (logical(p_bedflat) && ishandle(p_bedflat))
                set(p_bedflat, 'markersize', 24)
            end
            set(status_box, 'string', 'Bed selected.')
            
        else
            
            if (logical(p_pk(curr_layer)) && ishandle(p_pk(curr_layer)))
                set(p_pk(curr_layer), 'markersize', 24)
            end
            if (logical(p_pkdepth(curr_layer)) && ishandle(p_pkdepth(curr_layer)))
                set(p_pkdepth(curr_layer), 'markersize', 24)
            end
            if (logical(p_pkflat(curr_layer)) && ishandle(p_pkflat(curr_layer)))
                set(p_pkflat(curr_layer), 'markersize', 24)
            end
            if smooth_done(curr_layer)
                if (logical(p_pksmooth(curr_layer)) && ishandle(p_pksmooth(curr_layer)))
                    set(p_pksmooth(curr_layer), 'markersize', 24)
                end
                if (logical(p_pksmoothdepth(curr_layer)) && ishandle(p_pksmoothdepth(curr_layer)))
                    set(p_pksmoothdepth(curr_layer), 'markersize', 24)
                end
                if (logical(p_pksmoothflat(curr_layer)) && ishandle(p_pksmoothflat(curr_layer)))
                    set(p_pksmoothflat(curr_layer), 'markersize', 24)
                end
            end
            set(status_box, 'string', ['Layer #' num2str(curr_layer) ' selected.'])
            
        end
    end

%% Choose/select a layer interactively

    function pk_select_gui(source, eventdata)
        switch disp_type
            case {'twtt' '~depth'}
                ind_x_pk    = interp1(block.dist_lin(ind_decim), ind_decim, ind_x_pk, 'nearest', 'extrap');
            case 'flat'
                ind_x_pk    = ind_decim_flat(interp1(block.dist_lin(ind_decim_flat), 1:num_decim_flat, ind_x_pk, 'nearest', 'extrap'));
            otherwise
                return
        end
        switch disp_type
            case {'twtt' 'flat'}
                ind_y_pk   = interp1(block.twtt, 1:num_sample_trim, (1e-6 * ind_y_pk), 'nearest', 'extrap');
            case '~depth'
                ind_y_pk   = interp1(block.twtt, 1:num_sample_trim, (1e-6 * (ind_y_pk + (1e6 * (block.twtt_surf(ind_x_pk) - block.twtt(1))))), 'nearest', 'extrap');
        end
        tmp1                = NaN(pk.num_layer, 1);
        for ii = 1:pk.num_layer
            switch disp_type
                case {'twtt' '~depth'}
                    tmp1(ii)= pk.layer(ii).ind_y(ind_x_pk);
                case 'flat'
                    tmp1(ii)= ind_y_flat_mean(ii, interp1(ind_decim_flat, 1:num_decim_flat, ind_x_pk, 'nearest', 'extrap'));
            end
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
            case '~depth'
                if bed_avail
                    tmp1    = [tmp1; NaN; (ind_bed(ind_x_pk) - ind_surf(ind_x_pk) + 1)];
                end
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
        [tmp1, tmp3]        = unique(tmp1);
        if (length(tmp1) > 1)
            curr_layer      = interp1(tmp1, tmp2(tmp3), ind_y_pk, 'nearest', 'extrap');
        else
            curr_layer      = tmp2;
        end
        if (isempty(curr_layer) || any(curr_layer < 1) || (length(curr_layer) > 1))
            curr_layer      = 1;
        end
        set(layer_list, 'value', curr_layer)
        pk_select
    end

%% Focus on a layer vertically

    function pk_focus(source, eventdata)
        if ~pk_done
            set(status_box, 'string', 'No picked layers to focus on.')
            return
        end
        axes(ax_radar)
        switch disp_type
            case 'twtt'
                if (curr_layer == (pk.num_layer + 1))
                    ylim(1e6 .* [min(block.twtt_surf(ind_decim)) max(block.twtt_surf(ind_decim))])
                elseif (curr_layer == (pk.num_layer + 2))
                    ylim(1e6 .* [min(block.twtt_bed(ind_decim)) max(block.twtt_bed(ind_decim))])
                else
                    ylim(1e6 .* block.twtt(round([min(pk.layer(curr_layer).ind_y(ind_decim)) max(pk.layer(curr_layer).ind_y(ind_decim))])))
                end
            case '~depth'
                if (curr_layer == (pk.num_layer + 2))
                    ylim(1e6 .* [min(block.twtt_bed(ind_decim) - block.twtt(ind_surf(ind_decim) + 1)) max(block.twtt_bed(ind_decim) - block.twtt(ind_surf(ind_decim) + 1))])
                elseif (curr_layer <= pk.num_layer)
                    ylim(1e6 .* block.twtt(round([min(pk.layer(curr_layer).ind_y(ind_decim) - ind_surf(ind_decim) + 1) max(pk.layer(curr_layer).ind_y(ind_decim) - ind_surf(ind_decim) + 1)])))
                end
            case 'flat'
                if (curr_layer == (pk.num_layer + 1))
                    ylim(1e6 .* [min(block.twtt(ind_surf_flat(ind_decim(~isnan(ind_surf_flat(ind_decim)))))) max(block.twtt(ind_surf_flat(ind_decim(~isnan(ind_surf_flat(ind_decim))))))])
                elseif (curr_layer == (pk.num_layer + 2))
                    ylim(1e6 .* [min(block.twtt(ind_bed_flat(ind_decim(~isnan(ind_bed_flat(ind_decim)))))) max(block.twtt(ind_bed_flat(ind_decim(~isnan(ind_bed_flat(ind_decim))))))])
                else
                    ylim(1e6 .* block.twtt(round([min(ind_y_flat_mean(curr_layer, :)) max(ind_y_flat_mean(curr_layer, :))])))
                end
            otherwise
                return
        end
        tmp1                = get(ax_radar, 'ylim');
        [tmp1(1), tmp1(2)]  = deal((tmp1(1) - diff(tmp1)), (tmp1(2) + diff(tmp1)));
        if (tmp1(1) < (1e6 * twtt_min_ref))
            tmp1(1)         = 1e6 * twtt_min_ref;
        end
        if (tmp1(2) > (1e6 * twtt_max_ref))
            tmp1(2)         = 1e6 * twtt_max_ref;
        end
        ylim(tmp1)
        [twtt_min, twtt_max]= deal((1e-6 * tmp1(1)), (1e-6 * tmp1(2)));
        if ((1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref))) < get(twtt_min_slide, 'min'))
            set(twtt_min_slide, 'value', get(twtt_min_slide, 'min'))
        elseif ((1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref))) > get(twtt_min_slide, 'max'))
            set(twtt_min_slide, 'value', get(twtt_min_slide, 'max'))
        else
            set(twtt_min_slide, 'value', (1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref))))
        end
        if ((1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref))) < get(twtt_max_slide, 'min'))
            set(twtt_max_slide, 'value', get(twtt_max_slide, 'min'))
        elseif ((1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref))) > get(twtt_max_slide, 'max'))
            set(twtt_max_slide, 'value', get(twtt_max_slide, 'max'))
        else
            set(twtt_max_slide, 'value', (1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref))))
        end
        set(twtt_min_edit, 'string', sprintf('%3.1f', (1e6 * twtt_min)))
        set(twtt_max_edit, 'string', sprintf('%3.1f', (1e6 * twtt_max)))

        if (curr_layer == (pk.num_layer + 1))
            set(status_box, 'string', 'Focused on surface.')
        elseif (curr_layer == (pk.num_layer + 2))
            set(status_box, 'string', 'Focused on bed.')
        else
            set(status_box, 'string', ['Focused on layer #' num2str(curr_layer) '.'])
        end
        narrow_cb
    end

%% Switch to previous layer in list

    function pk_last(source, eventdata)
        if (curr_layer > 1)
            curr_layer      = curr_layer - 1;
            set(layer_list, 'value', curr_layer)
            pk_select
        end
    end

%% Switch to next layer in the list

    function pk_next(source, eventdata)
        if (curr_layer < pk.num_layer)
            curr_layer      = curr_layer + 1;
            set(layer_list, 'value', curr_layer)
            pk_select
        end
    end

%% Delete layer

    function pk_del(source, eventdata)
        if ~(pk_done && pk.num_layer && curr_layer)
            set(status_box, 'string', 'No picked layers to delete yet.')
            return
        end
        set(status_box, 'string', 'Delete current layer? Y: yes; otherwise: no.')
        waitforbuttonpress
        if ~strcmpi(get(pkgui, 'currentcharacter'), 'Y')
            set(status_box, 'string', 'Layer deletion cancelled.')
            return
        end
        pk_del_breakout

        if (curr_layer == (pk.num_layer + 1))
            set(status_box, 'string', 'Surface deleted.')
        elseif (curr_layer == (pk.num_layer + 2))
            set(status_box, 'string', 'Bed deleted.')
        else
            set(status_box, 'string', ['Layer #' num2str(curr_layer) ' deleted.'])
        end
        pause(0.1)
        pk_select
    end

    function pk_del_breakout(source, eventdata)
        if (curr_layer == (pk.num_layer + 1))
            if (logical(p_surf) && ishandle(p_surf))
                delete(p_surf)
            end
            if (logical(p_surfflat) && ishandle(p_surfflat))
                delete(p_surfflat)
            end
            [ind_surf, ind_surf_flat, block.twtt_surf] ...
                                = deal(NaN(1, block.num_trace));
            do_surfbed          = true;
            surf_avail          = false;
        elseif (curr_layer == (pk.num_layer + 2))
            if (logical(p_bed) && ishandle(p_bed))
                delete(p_bed)
            end
            if (logical(p_beddepth) && ishandle(p_beddepth))
                delete(p_beddepth)
            end
            if (logical(p_bedflat) && ishandle(p_bedflat))
                delete(p_bedflat)
            end
            [ind_bed, ind_bed_flat, block.twtt_bed] ...
                                = deal(NaN(1, block.num_trace));
            do_surfbed          = true;
            bed_avail           = false;            
        else
            tmp1                = setdiff(1:pk.num_layer, curr_layer);
            if (logical(p_pk(curr_layer)) && ishandle(p_pk(curr_layer)))
                delete(p_pk(curr_layer))
            end
            if (logical(p_pkdepth(curr_layer)) && ishandle(p_pkdepth(curr_layer)))
                delete(p_pkdepth(curr_layer))
            end
            if (logical(p_pkflat(curr_layer)) && ishandle(p_pkflat(curr_layer)))
                delete(p_pkflat(curr_layer))
            end
            if (logical(p_pksmooth(curr_layer)) && ishandle(p_pksmooth(curr_layer)))
                delete(p_pksmooth(curr_layer))
            end
            if (logical(p_pksmoothdepth(curr_layer)) && ishandle(p_pksmoothdepth(curr_layer)))
                delete(p_pksmoothdepth(curr_layer))
            end
            if (logical(p_pksmoothflat(curr_layer)) && ishandle(p_pksmoothflat(curr_layer)))
                delete(p_pksmoothflat(curr_layer))
            end
            [pk.num_layer, pk.layer, smooth_done, p_pk, p_pkdepth, p_pkflat, p_pksmooth, p_pksmoothdepth, p_pksmoothflat, ind_y_flat_mean, ind_y_flat_smooth] ...
                            = deal((pk.num_layer - 1), pk.layer(tmp1), smooth_done(tmp1), p_pk(tmp1), p_pkdepth(tmp1), p_pkflat(tmp1), p_pksmooth(tmp1), p_pksmoothdepth(tmp1), p_pksmoothflat(tmp1), ind_y_flat_mean(tmp1, :), ind_y_flat_smooth(tmp1, :));
            curr_layer      = curr_layer - 1;
            if (isempty(curr_layer) || any(curr_layer < 1))
                curr_layer  = 1;
            end
            if pk.num_layer
                set(layer_list, 'string', [num2cell(1:pk.num_layer) 'surface' 'bed'], 'value', curr_layer)
            else
                set(layer_list, 'string', {'surface' 'bed'}, 'value', 1)
            end
            match_done      = false;
        end
    end

%% Adjust/edit current layer

    function pk_adj(source, eventdata)
        
        if ~(pk_done || surf_avail || bed_avail)
            set(status_box, 'string', 'No picked layers to adjust yet.')
            return
        end
        
        axes(ax_radar)
        tmp5                = 0;
        
        set(pkgui, 'keypressfcn', [], 'windowbuttondownfcn', [])
        while true
            
            set(status_box, 'string', 'L: delete left; R: delete right; C: cut (left start); U: undo; Q: quit.')
            
            % get pick and convert to indices
            [ind_x_pk, ~, button] ...
                            = ginput(1);
            switch disp_type
                case {'twtt' '~depth'}
                    ind_x_pk= interp1(block.dist_lin(ind_decim), ind_decim, ind_x_pk, 'nearest', 'extrap');
                case 'flat'
                    ind_x_pk= interp1(block.dist_lin(ind_decim_flat), ind_decim_flat, ind_x_pk, 'nearest', 'extrap');
            end
            
            if (strcmpi(char(button), 'L') || strcmpi(char(button), 'R') || strcmpi(char(button), 'C'))
                
                switch char(button)
                    case {'l' 'L'}
                        tmp1= cell(1, 2);
                        [tmp1{1}, tmp1{2}] ...
                            = deal(1:ind_x_pk, 1:interp1(ind_decim_flat, 1:num_decim_flat, ind_x_pk, 'nearest', 'extrap'));
                    case {'r' 'R'}
                        tmp1= cell(1, 2);
                        [tmp1{1}, tmp1{2}] ...
                            = deal(ind_x_pk:block.num_trace, interp1(ind_decim_flat, 1:num_decim_flat, ind_x_pk, 'nearest', 'extrap'):num_decim_flat);
                    case {'c' 'C'}
                        set(status_box, 'string', 'Now choose right end of cut...')
                        [tmp1, ~] = ginput(1);
                        switch disp_type
                            case {'twtt' '~depth'}
                                tmp1 = interp1(block.dist_lin(ind_decim), ind_decim, tmp1, 'nearest', 'extrap');
                            case 'flat'
                                tmp1 = ind_decim_flat(interp1(block.dist_lin(ind_decim_flat), 1:num_decim_flat, tmp1, 'nearest', 'extrap'));
                        end
                        tmp4= tmp1; % in case of undo
                        tmp1= cell(1, 2);
                        [tmp1{1}, tmp1{2}] ...
                            = deal(ind_x_pk:tmp4, interp1(ind_decim_flat, 1:num_decim_flat, ind_x_pk, 'nearest', 'extrap'):interp1(ind_decim_flat, 1:num_decim_flat, tmp4, 'nearest', 'extrap'));
                end
                
                [tmp2, tmp3]= deal(button, cell(1, 4));
                
                if (curr_layer == (pk.num_layer + 1))
                    [tmp3{1}, tmp3{3}] ...
                            = deal(ind_surf(tmp1{1}), block.twtt_surf(tmp1{1}));
                    [block.twtt_surf(tmp1{1}), ind_surf(tmp1{1})] ...
                            = deal(NaN);
                    if flat_done
                        tmp3{2} ...
                            = ind_surf_flat(tmp1{1});
                        ind_surf_flat(tmp1{1}) ...
                            = NaN;
                    end
                elseif (curr_layer == (pk.num_layer + 2))
                    [tmp3{1}, tmp3{3}] ...
                            = deal(ind_bed(tmp1{1}), block.twtt_bed(tmp1{1}));
                    [block.twtt_bed(tmp1{1}), ind_bed(tmp1{1})] ...
                            = deal(NaN);
                    if flat_done
                        tmp3{2} ...
                            = ind_bed_flat(tmp1{1});
                        ind_bed_flat(tmp1{1}) ...
                            = NaN;
                    end
                else
                    tmp3{1} = pk.layer(curr_layer).ind_y(tmp1{1});
                    pk.layer(curr_layer).ind_y(tmp1{1}) ...
                            = NaN;
                    tmp3{2} = ind_y_flat_mean(curr_layer, tmp1{2});
                    ind_y_flat_mean(curr_layer, tmp1{2}) ...
                            = NaN;
                    tmp3{3} = ind_y_flat_smooth(curr_layer, tmp1{2});
                    ind_y_flat_smooth(curr_layer, tmp1{2}) ...
                            = NaN;
                    if smooth_done(curr_layer)
                        tmp3{4} ...
                            = deal(pk.layer(curr_layer).ind_y_smooth(tmp1{1}));
                        pk.layer(curr_layer).ind_y_smooth(tmp1{1}) ...
                            = NaN;
                    end
                end
                
                pk_fix
                
                tmp5        = ind_x_pk;
                switch char(button)
                    case {'l' 'L'}
                        if (curr_layer == (pk.num_layer + 1))
                            set(status_box, 'string', ['Surface cut left at ' num2str(block.dist_lin(ind_x_pk), '%3.1f') ' km.'])
                        elseif (curr_layer == (pk.num_layer + 2))
                            set(status_box, 'string', ['Bed cut left at ' num2str(block.dist_lin(ind_x_pk), '%3.1f') ' km.'])
                        else
                            set(status_box, 'string', ['Layer #' num2str(curr_layer) ' cut left at ' num2str(block.dist_lin(ind_x_pk), '%3.1f') ' km.'])
                        end
                    case {'r' 'R'}
                        if (curr_layer == (pk.num_layer + 1))
                            set(status_box, 'string', ['Surface cut right at ' num2str(block.dist_lin(ind_x_pk), '%3.1f') ' km.'])
                        elseif (curr_layer == (pk.num_layer + 2))
                            set(status_box, 'string', ['Bed cut right at ' num2str(block.dist_lin(ind_x_pk), '%3.1f') ' km.'])
                        else
                            set(status_box, 'string', ['Layer #' num2str(curr_layer) ' cut right at ' num2str(block.dist_lin(ind_x_pk), '%3.1f') ' km.'])
                        end
                    case {'c' 'C'}
                        if (curr_layer == (pk.num_layer + 1))
                            set(status_box, 'string', ['Surface cut between ' num2str(block.dist_lin(ind_x_pk), '%3.1f') ' and ' num2str(block.dist_lin(tmp1(find(~isnan(tmp1), 1, 'last'))), '%3.1f') ' km.'])
                        elseif (curr_layer == (pk.num_layer + 2))
                            set(status_box, 'string', ['Bed cut between ' num2str(block.dist_lin(ind_x_pk), '%3.1f') ' and ' num2str(block.dist_lin(tmp1(find(~isnan(tmp1), 1, 'last'))), '%3.1f') ' km.'])
                        else
                            set(status_box, 'string', ['Layer #' num2str(curr_layer) ' cut between ' num2str(block.dist_lin(ind_x_pk), '%3.1f') ' and ' num2str(block.dist_lin(tmp1(find(~isnan(tmp1), 1, 'last'))), '%3.1f') ' km.'])
                        end
                end
                
            elseif strcmpi(char(button), 'U') % undo adjustment done in current set
                
                if ~tmp5
                    continue
                end
                
                if (strcmpi(char(tmp2), 'L') || strcmpi(char(tmp2), 'R') || strcmpi(char(tmp2), 'C'))
                    
                    tmp1        = cell(1, 2);
                    switch char(tmp2)
                        case {'l' 'L'}
                            [tmp1{1}, tmp1{2}] ...
                                = deal(1:tmp5, 1:interp1(ind_decim_flat, 1:num_decim_flat, tmp5, 'nearest', 'extrap'));
                        case {'r' 'R'}
                            [tmp1{1}, tmp1{2}] ...
                                = deal(tmp5:block.num_trace, interp1(ind_decim_flat, 1:num_decim_flat, tmp5, 'nearest', 'extrap'):num_decim_flat);
                        case {'c' 'C'}
                            [tmp1{1}, tmp1{2}] ...
                                = deal(tmp5:tmp4, interp1(ind_decim_flat, 1:num_decim_flat, tmp5, 'nearest', 'extrap'):interp1(ind_decim_flat, 1:num_decim_flat, tmp4, 'nearest', 'extrap'));
                    end
                    if (curr_layer == (pk.num_layer + 1))
                        [ind_surf(tmp1{1}), block.twtt_surf(tmp1{1})] ...
                            = deal(tmp3{1}, tmp3{3});
                        if flat_done
                            ind_surf_flat(tmp1{1}) ...
                                = tmp3{2};
                        end
                    elseif (curr_layer == (pk.num_layer + 2))
                        [ind_bed(tmp1{1}), block.twtt_bed(tmp1{1})] ...
                            = deal(tmp3{1}, tmp3{3});
                        if flat_done
                            ind_bed_flat(tmp1{1}) ...
                                = tmp3{2};
                        end
                    else
                        pk.layer(curr_layer).ind_y(tmp1{1}) ...
                            = tmp3{1};
                        ind_y_flat_mean(curr_layer, tmp1{2}) ...
                            = tmp3{2};
                        ind_y_flat_smooth(curr_layer, tmp1{2}) ...
                            = tmp3{3};
                        if smooth_done(curr_layer)
                            pk.layer(curr_layer).ind_y_smooth(tmp1{1}) ...
                                = tmp3{4};
                        end
                    end
                end
                
                tmp3        = 0;
                pk_fix
                set(status_box, 'string', 'Undid previous adjustment.')
                
            elseif strcmpi(char(button), 'Q')

                if (curr_layer == (pk.num_layer + 1))
                    set(status_box, 'string', 'Done adjusting surface.')
                elseif (curr_layer == (pk.num_layer + 2))
                    set(status_box, 'string', 'Done adjusting bed.')
                else
                    pk_sort
                    tmp1    = find((tmp1 == curr_layer), 1); % tmp1 now from pk_sort
                    if ((tmp1 ~= curr_layer) && ~isempty(tmp1))
                        set(status_box, 'string', ['Done adjusting layer #' num2str(curr_layer) ' (now layer #' num2str(tmp1) ').'])
                        curr_layer = tmp1;
                    else
                        set(status_box, 'string', ['Done adjusting layer #' num2str(curr_layer) '.'])
                    end
                    if (isempty(curr_layer) || any(curr_layer < 1))
                        curr_layer = 1;
                    end
                    set(layer_list, 'value', curr_layer)
                    pk_select
                    match_done ...
                            = false;
                end
                set(pkgui, 'keypressfcn', @keypress, 'windowbuttondownfcn', @mouse_click)
                return
            end
        end
    end

%% Fix layer based on interactive adjustments

    function pk_fix(source, eventdata)
        
        if (curr_layer == (pk.num_layer + 1))
            tmp1            = ind_decim(~isnan(ind_surf(ind_decim)));
        elseif (curr_layer == (pk.num_layer + 2))
            tmp1            = ind_decim(~isnan(ind_bed(ind_decim)));
        else
            tmp1            = ind_decim(~isnan(pk.layer(curr_layer).ind_y(ind_decim)));
        end
        if isempty(tmp1)
            pk_del
            set(status_box, 'string', 'Layer now empty so it was deleted.')
            return
        end
        
        if (curr_layer == (pk.num_layer + 1))
            
            if (any(p_pkdepth) && any(ishandle(p_pkdepth)))
                delete(p_pkdepth(logical(p_pkdepth) & ishandle(p_pkdepth)))
            end
            if (any(p_pksmoothdepth) && any(ishandle(p_pksmoothdepth)))
                delete(p_pksmoothdepth(logical(p_pksmoothdepth) & ishandle(p_pksmoothdepth)))
            end
            if (logical(p_surf) && ishandle(p_surf))
                delete(p_surf)
                p_surf      = plot(block.dist_lin(tmp1), (1e6 .* block.twtt_surf(tmp1)), 'm.', 'markersize', 24, 'visible', 'off');
            end
            if flat_done
                if (logical(p_surfflat) && ishandle(p_surfflat))
                    delete(p_surfflat)
                    p_surfflat ...
                            = plot(block.dist_lin(ind_decim_flat(~isnan(ind_surf_flat(ind_decim_flat)))), (1e6 .* block.twtt(ind_surf_flat(ind_decim_flat(~isnan(ind_surf_flat(ind_decim_flat)))))), 'm.', 'markersize', 24, 'visible', 'off');
                end
            end
            if ~isempty(find(~isnan(ind_surf(ind_decim)), 1))
                amp_depth   = NaN(num_sample_trim, num_decim, 'single');
                for ii = find(~isnan(ind_surf(ind_decim)))
                    amp_depth(1:(num_sample_trim - ind_surf(ind_decim(ii)) + 1), ii) ...
                            = amp_mean(ind_surf(ind_decim(ii)):num_sample_trim, ii); % shift data up to surface
                end
                set(disp_check(2), 'visible', 'on')
                depth_avail = true;
                if bed_avail
                    if (logical(p_beddepth) && ishandle(p_beddepth))
                        delete(p_beddepth)
                    end
                    tmp1    = ind_decim(~isnan(ind_bed(ind_decim)) & ~isnan(ind_surf(ind_decim)));
                    p_beddepth ...
                            = plot(block.dist_lin(tmp1), (1e6 .* (block.twtt_bed(tmp1) - block.twtt_surf(tmp1) + block.twtt(1))), 'm.', 'markersize', 12, 'visible', 'off');
                end
                amp_depth   = NaN(num_sample_trim, num_decim, 'single');
                for ii = find(~isnan(ind_surf(ind_decim)))
                    amp_depth(1:(num_sample_trim - ind_surf(ind_decim(ii)) + 1), ii) ...
                            = amp_mean(ind_surf(ind_decim(ii)):num_sample_trim, ii); % shift data up to surface
                end
                set(disp_check(2), 'visible', 'on')
                depth_avail = true;
                for ii = 1:pk.num_layer
                    if ~isempty(find((~isnan(pk.layer(ii).ind_y(ind_decim)) & ~isnan(ind_surf(ind_decim))), 1))
                        tmp1= ind_decim(~isnan(pk.layer(ii).ind_y(ind_decim)) & ~isnan(ind_surf(ind_decim)));
                        tmp2= pk.layer(ii).ind_y(tmp1) - ind_surf(tmp1) + 1;
                        tmp2((tmp2 < 1) | (tmp2 > num_sample_trim)) ...
                            = NaN;
                        p_pkdepth(ii) ...
                            = plot(block.dist_lin(tmp1(~isnan(tmp2))), (1e6 .* block.twtt(round(tmp2(~isnan(tmp2))))), 'r.', 'markersize', 12, 'visible', 'off');
                    end
                    if ~isempty(find((~isnan(pk.layer(ii).ind_y_smooth(ind_decim)) & ~isnan(ind_surf(ind_decim))), 1))
                        tmp1= ind_decim(~isnan(pk.layer(ii).ind_y_smooth(ind_decim)) & ~isnan(ind_surf(ind_decim)));
                        tmp2= pk.layer(ii).ind_y_smooth(tmp1) - ind_surf(tmp1) + 1;
                        tmp2((tmp2 < 1) | (tmp2 > num_sample_trim)) ...
                            = NaN;
                        p_pksmoothdepth(ii) ...
                            = plot(block.dist_lin(tmp1(~isnan(tmp2))), (1e6 .* block.twtt(round(tmp2(~isnan(tmp2))))), 'g.', 'markersize', 12, 'visible', 'off');
                    end
                end
                if strcmp(disp_type, '~depth')
                    plot_depth
                end
            end
            do_surfbed      = true;
            show_surfbed
            
        elseif (curr_layer == (pk.num_layer + 2))
            
            if (logical(p_bed) && ishandle(p_bed))
                delete(p_bed)
                p_bed       = plot(block.dist_lin(tmp1), (1e6 .* block.twtt_bed(tmp1)), 'm.', 'markersize', 24, 'visible', 'off');
            end
            if (logical(p_beddepth) && ishandle(p_beddepth))
                delete(p_beddepth)
                if (bed_avail && depth_avail)
                    tmp1    = ind_decim(~isnan(ind_bed(ind_decim)) & ~isnan(ind_surf(ind_decim)));
                    p_beddepth ...
                            = plot(block.dist_lin(tmp1), (1e6 .* (block.twtt_bed(tmp1) - block.twtt_surf(tmp1) + block.twtt(1))), 'm.', 'markersize', 12, 'visible', 'off');
                end
            end
            if flat_done
                if (logical(p_bedflat) && ishandle(p_bedflat))
                    delete(p_bedflat)
                    p_bedflat ...
                            = plot(block.dist_lin(ind_decim_flat(~isnan(ind_bed_flat(ind_decim_flat)))), (1e6 .* block.twtt(ind_bed_flat(ind_decim_flat(~isnan(ind_bed_flat(ind_decim_flat)))))), 'm.', 'markersize', 24, 'visible', 'off');
                end
            end
            do_surfbed      = true;
            show_surfbed
            
        else
            
            if (logical(p_pk(curr_layer)) && ishandle(p_pk(curr_layer)))
                delete(p_pk(curr_layer))
                p_pk(curr_layer) ...
                            = plot(block.dist_lin(tmp1), (1e6 .* block.twtt(round(pk.layer(curr_layer).ind_y(tmp1)))), 'r.', 'markersize', 24, 'visible', 'off');
            end
            if (logical(p_pkdepth(curr_layer)) && ishandle(p_pkdepth(curr_layer)))
                delete(p_pkdepth(curr_layer))
                if depth_avail
                    tmp1    = ind_decim(~isnan(pk.layer(curr_layer).ind_y(ind_decim)) & ~isnan(ind_surf(ind_decim)));
                    tmp2    = pk.layer(curr_layer).ind_y(tmp1) - ind_surf(tmp1) + 1;
                    tmp2((tmp2 < 1) | (tmp2 > num_sample_trim)) ...
                            = NaN;
                    p_pkdepth(curr_layer) ...
                            = plot(block.dist_lin(tmp1(~isnan(tmp2))), (1e6 .* block.twtt(round(tmp2(~isnan(tmp2))))), 'r.', 'markersize', 24, 'visible', 'off');
                end
            end
            if flat_done
                if (logical(p_pkflat(curr_layer)) && ishandle(p_pkflat(curr_layer)))
                    delete(p_pkflat(curr_layer))
                    p_pkflat(curr_layer) ...
                            = plot(block.dist_lin(ind_decim_flat(~isnan(ind_y_flat_mean(curr_layer, :)))), (1e6 .* block.twtt(ind_y_flat_mean(curr_layer, ~isnan(ind_y_flat_mean(curr_layer, :))))), 'r.', 'markersize', 24, 'visible', 'off');
                end
            end
            show_pk
            if smooth_done(curr_layer)
                if (logical(p_pksmooth(curr_layer)) && ishandle(p_pksmooth(curr_layer)))
                    delete(p_pksmooth(curr_layer))
                end
                if (logical(p_pksmoothdepth(curr_layer)) && ishandle(p_pksmoothdepth(curr_layer)))
                    delete(p_pksmoothdepth(curr_layer))
                end
                if (logical(p_pksmoothflat(curr_layer)) && ishandle(p_pksmoothflat(curr_layer)))
                    delete(p_pksmoothflat(curr_layer))
                end
                tmp1        = ind_decim(~isnan(pk.layer(curr_layer).ind_y_smooth(ind_decim)));
                p_pksmooth(curr_layer) ...
                            = plot(block.dist_lin(tmp1), (1e6 .* block.twtt(round(pk.layer(curr_layer).ind_y_smooth(tmp1)))), 'g.', 'markersize', 24, 'visible', 'off');
                if depth_avail
                    tmp1    = ind_decim(~isnan(pk.layer(curr_layer).ind_y_smooth(ind_decim)) & ~isnan(ind_surf(ind_decim)));
                    tmp2    = pk.layer(curr_layer).ind_y_smooth(tmp1) - ind_surf(tmp1) + 1;
                    tmp2((tmp2 < 1) | (tmp2 > num_sample_trim)) ...
                            = NaN;
                    p_pksmoothdepth(curr_layer) ...
                            = plot(block.dist_lin(tmp1(~isnan(tmp2))), (1e6 .* block.twtt(round(tmp2(~isnan(tmp2))))), 'g.', 'markersize', 24, 'visible', 'off');
                end
                if flat_done
                    tmp1    = find(~isnan(ind_y_flat_smooth(curr_layer, :)));
                    p_pksmoothflat(curr_layer) ...
                            = plot(block.dist_lin(ind_decim_flat(tmp1)), (1e6 .* block.twtt(round(ind_y_flat_smooth(curr_layer, tmp1)))), 'g.', 'markersize', 24, 'visible', 'off');
                end
                show_smooth
            end
        end
        
        pk_select
        tmp2                = button;
    end

%% Merge two layers

    function pk_merge(source, eventdata)
        
        if ~any(strcmp(disp_type, {'twtt' '~depth' 'flat'}))
            set(status_box, 'string', 'Must be in twtt, ~depth or flat to merge layers.')
            return
        end
        if ~(pk_done || surf_avail || bed_avail)
            set(status_box, 'string', 'No picked layers to merge yet.')
            return
        end
        if ((pk.num_layer < 2) && ((pk.num_layer <= 1) && ~(surf_avail || bed_avail)))
            set(status_box, 'string', 'Not enough layers to merge.')
            return
        end
        if ~all(smooth_done)
            set(status_box, 'string', 'All layers must be smoothed prior to merging...')
            pk_smooth
            pk_merge
        end
        
        if (curr_layer == (pk.num_layer + 1))
            set(status_box, 'string', 'Pick layer to merge with surface (Q: cancel)...')
        elseif (curr_layer == (pk.num_layer + 2))
            set(status_box, 'string', 'Pick layer to merge with bed (Q: cancel)...')
        else
            set(status_box, 'string', ['Pick layer to merge with layer #' num2str(curr_layer) ' (Q: cancel)...'])
        end
        
        axes(ax_radar)
        
        % get pick and convert to indices
        [ind_x_pk, ind_y_pk, button] ...
                            = ginput(1);
        if strcmpi(char(button), 'Q')
            set(status_box, 'string', 'Layer merging cancelled.')
            return
        end
        
        switch disp_type
            case {'twtt' '~depth'}
                ind_x_pk    = interp1(block.dist_lin(ind_decim), ind_decim, ind_x_pk, 'nearest', 'extrap');
            case 'flat'
                ind_x_pk    = ind_decim_flat(interp1(block.dist_lin(ind_decim_flat), 1:num_decim_flat, ind_x_pk, 'nearest', 'extrap'));
        end
        switch disp_type
            case {'twtt' 'flat'}
                ind_y_pk    = interp1(block.twtt, 1:num_sample_trim, (1e-6 * ind_y_pk), 'nearest', 'extrap');
            case '~depth'
                ind_y_pk    = interp1(block.twtt, 1:num_sample_trim, (1e-6 * (ind_y_pk + (1e6 * (block.twtt_surf(ind_x_pk) - block.twtt(1))))), 'nearest', 'extrap');
        end
        
        % get current layer positions at ind_x_pk, depending on what we're working with
        tmp1                = NaN(pk.num_layer, 1);
        for ii = 1:pk.num_layer
            switch disp_type
                case {'twtt' '~depth'}
                    tmp1(ii)= pk.layer(ii).ind_y(ind_x_pk);
                case 'flat'
                    tmp1(ii)= ind_y_flat_mean(ii, interp1(ind_decim_flat, 1:num_decim_flat, ind_x_pk, 'nearest', 'extrap'));
            end
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
            case '~depth'
                if bed_avail
                    tmp1    = [tmp1; (ind_bed(ind_x_pk) - ind_surf(ind_x_pk) + 1)];
                end
            case 'flat'
                if surf_avail
                    tmp1    = [tmp1; ind_surf_flat(interp1(ind_decim_flat, 1:num_decim_flat, ind_x_pk, 'nearest', 'extrap'))];
                else
                    tmp1    = [tmp1; NaN];
                end
                if bed_avail
                    tmp1    = [tmp1; ind_bed_flat(interp1(ind_decim_flat, 1:num_decim_flat, ind_x_pk, 'nearest', 'extrap'))];
                end
        end
        tmp2                = find(~isnan(tmp1));
        tmp1                = tmp1(tmp2);
        [tmp1, tmp3]        = unique(tmp1); % test for repeating tmp1 values
        if (length(tmp1) > 1)
            tmp1            = interp1(tmp1, tmp2(tmp3), ind_y_pk, 'nearest', 'extrap');
        else
            tmp1            = tmp2;
        end
        if (tmp1 == curr_layer)
            set(status_box, 'string', 'Aborted merging because picked layer is the same as current layer.')
            return
        end
        if isempty(tmp1)
            set(status_box, 'string', 'No layer picked to be merged. Pick more precisely.')
            return
        end
        if (((tmp1 == (pk.num_layer + 1)) && (curr_layer == (pk.num_layer + 2))) || ((tmp1 == (pk.num_layer + 2)) && (curr_layer == (pk.num_layer + 1))))
            set(status_box, 'string', 'Cannot merge surface and bed.')
            return
        end
        
        % force curr_layer to be surface or bed
        if (any(tmp1 == (pk.num_layer + [1 2])) && (curr_layer <= pk.num_layer))
            [curr_layer, tmp1] ...
                            = deal(tmp1, curr_layer);
        end
        
        if (curr_layer == (pk.num_layer + 1))
            
            % replace NaN values in surface with those in the second layer
            ind_surf(isnan(ind_surf)) ...
                            = pk.layer(tmp1).ind_y_smooth(isnan(ind_surf));
            block.twtt_surf(isnan(block.twtt_surf) & ~isnan(pk.layer(tmp1).ind_y_smooth)) ...
                            = block.twtt(round(pk.layer(tmp1).ind_y_smooth(isnan(block.twtt_surf) & ~isnan(pk.layer(tmp1).ind_y_smooth))));
            if flat_done
                ind_surf_flat= NaN(1, block.num_trace);
                for ii = find(~sum(isnan(ind_y_flat)))
                    [~, tmp2] ...
                            = unique(ind_y_flat(:, ii));
                    ind_surf_flat(ii) ...
                            = interp1(ind_y_flat(tmp2, ii), tmp2, ind_surf(ii), 'nearest', 'extrap');
                end
            end
            
            % fix plots
            tmp2            = ind_decim(~isnan(ind_surf(ind_decim)));
            if (any([p_pk(tmp1) p_surf]) && any(ishandle([p_pk(tmp1) p_surf])))
                delete([p_pk(tmp1) p_surf])
                p_surf      = plot(block.dist_lin(tmp2), (1e6 .* block.twtt_surf(tmp2)), 'm.', 'markersize', 24, 'visible', 'off');
            end
            if flat_done
                if (any([p_pkflat(tmp1) p_surfflat]) && any(ishandle([p_pkflat(tmp1) p_surfflat])))
                    delete([p_pkflat(tmp1) p_surfflat])
                    p_surfflat ...
                            = plot(block.dist_lin(ind_decim_flat(~isnan(ind_surf_flat(ind_decim_flat)))), (1e6 .* block.twtt(ind_surf_flat(ind_decim_flat(~isnan(ind_surf_flat(ind_decim_flat)))))), 'm.', 'markersize', 24, 'visible', 'off');
                end
            end
            
            if (logical(p_pksmooth(tmp1)) && ishandle(p_pksmooth(tmp1)))
                delete(p_pksmooth(tmp1))
            end
            if (logical(p_pksmoothflat(tmp1)) && ishandle(p_pksmoothflat(tmp1)))
                delete(p_pksmoothflat(tmp1))
            end
            if (any(p_pkdepth) && any(ishandle(p_pkdepth)))
                delete(p_pkdepth(logical(p_pkdepth) & ishandle(p_pkdepth)))
            end
            if (any(p_pksmoothdepth) && any(ishandle(p_pksmoothdepth)))
                delete(p_pksmoothdepth(logical(p_pksmoothdepth) & ishandle(p_pksmoothdepth)))
            end
            
            amp_depth       = NaN(num_sample_trim, num_decim, 'single');
            for ii = find(~isnan(ind_surf(ind_decim)))
                amp_depth(1:(num_sample_trim - ind_surf(ind_decim(ii)) + 1), ii) ...
                            = amp_mean(ind_surf(ind_decim(ii)):num_sample_trim, ii); % shift data up to surface
            end
            set(disp_check(2), 'visible', 'on')
            depth_avail     = true;
            [p_pkdepth, p_pksmoothdepth] ...
                            = deal(zeros(1, pk.num_layer));
            for ii = 1:pk.num_layer
                if ~isempty(find((~isnan(pk.layer(ii).ind_y(ind_decim)) & ~isnan(ind_surf(ind_decim))), 1))
                    tmp2    = ind_decim(~isnan(pk.layer(ii).ind_y(ind_decim)) & ~isnan(ind_surf(ind_decim)));
                    tmp3    = pk.layer(ii).ind_y(tmp2) - ind_surf(tmp2) + 1;
                    tmp3((tmp3 < 1) | (tmp3 > num_sample_trim)) ...
                            = NaN;
                    p_pkdepth(ii) ...
                            = plot(block.dist_lin(tmp2(~isnan(tmp3))), (1e6 .* block.twtt(round(tmp3(~isnan(tmp3))))), 'r.', 'markersize', 12, 'visible', 'off');
                end
                if ~isempty(find((~isnan(pk.layer(ii).ind_y_smooth(ind_decim)) & ~isnan(ind_surf(ind_decim))), 1))
                    tmp2    = ind_decim(~isnan(pk.layer(ii).ind_y_smooth(ind_decim)) & ~isnan(ind_surf(ind_decim)));
                    tmp3    = pk.layer(ii).ind_y_smooth(tmp2) - ind_surf(tmp2) + 1;
                    tmp3((tmp3 < 1) | (tmp3 > num_sample_trim)) ...
                            = NaN;
                    p_pksmoothdepth(ii) ...
                            = plot(block.dist_lin(tmp2(~isnan(tmp3))), (1e6 .* block.twtt(round(tmp3(~isnan(tmp3))))), 'g.', 'markersize', 12, 'visible', 'off');
                end
            end
            
            tmp3            = setdiff(1:pk.num_layer, tmp1);
            [pk.layer, smooth_done, p_pk, p_pkdepth, p_pkflat, p_pksmooth, p_pksmoothdepth, p_pksmoothflat, ind_y_flat_mean, ind_y_flat_smooth, pk.num_layer] ...
                            = deal(pk.layer(tmp3), smooth_done(tmp3), p_pk(tmp3), p_pkdepth(tmp3), p_pkflat(tmp3), p_pksmooth(tmp3), p_pksmoothdepth(tmp3), p_pksmoothflat(tmp3), ind_y_flat_mean(tmp3, :), ind_y_flat_smooth(tmp3, :), (pk.num_layer - 1));
            
            if ~isempty(find(~isnan(ind_surf(ind_decim)), 1))
                amp_depth   = NaN(num_sample_trim, num_decim, 'single');
                for ii = find(~isnan(ind_surf(ind_decim)))
                    amp_depth(1:(num_sample_trim - ind_surf(ind_decim(ii)) + 1), ii) ...
                            = amp_mean(ind_surf(ind_decim(ii)):num_sample_trim, ii); % shift data up to surface
                end
                set(disp_check(2), 'visible', 'on')
                depth_avail = true;
                if bed_avail
                    if (logical(p_beddepth) && ishandle(p_beddepth))
                        delete(p_beddepth)
                    end
                    tmp2    = ind_decim(~isnan(ind_bed(ind_decim)) & ~isnan(ind_surf(ind_decim)));
                    p_beddepth ...
                            = plot(block.dist_lin(tmp2), (1e6 .* (block.twtt_bed(tmp2) - block.twtt_surf(tmp2) + block.twtt(1))), 'm.', 'markersize', 12, 'visible', 'off');
                end
            end
            
            do_surfbed      = true;
            set(status_box, 'string', ['Surface and layer #' num2str(tmp1) ' merged.'])
            
            if (curr_layer > tmp1)
                curr_layer  = curr_layer - 1;
            end
            pk_sort
            curr_layer      = find((tmp1 == curr_layer), 1); % using tmp1 from pk_sort
            if (isempty(curr_layer) || any(curr_layer < 1))
                curr_layer  = 1;
            end
            
            if strcmp(disp_type, '~depth')
                plot_depth
            end
            
        elseif (curr_layer == (pk.num_layer + 2))
            
            % replace NaN values in bed with those in the second layer
            ind_bed(isnan(ind_bed)) ...
                            = pk.layer(tmp1).ind_y_smooth(isnan(ind_bed));
            block.twtt_bed(isnan(block.twtt_bed) & ~isnan(pk.layer(tmp1).ind_y_smooth)) ...
                            = block.twtt(round(pk.layer(tmp1).ind_y_smooth(isnan(block.twtt_bed) & ~isnan(pk.layer(tmp1).ind_y_smooth))));
            if flat_done
                ind_bed_flat= NaN(1, block.num_trace);
                for ii = find(~sum(isnan(ind_y_flat)))
                    [~, tmp2] ...
                            = unique(ind_y_flat(:, ii));
                    ind_bed_flat(ii) ...
                            = interp1(ind_y_flat(tmp2, ii), tmp2, ind_bed(ii), 'nearest', 'extrap');
                end
            end
            
            % fix plots
            tmp2            = ind_decim(~isnan(ind_bed(ind_decim)));
            if (any([p_pk(tmp1) p_bed]) && any(ishandle([p_pk(tmp1) p_bed])))
                delete([p_pk(tmp1) p_bed])
                p_bed       = plot(block.dist_lin(tmp2), (1e6 .* block.twtt_bed(tmp2)), 'm.', 'markersize', 24, 'visible', 'off');
            end
            if (bed_avail && depth_avail)
                if (any([p_pkdepth(tmp1) p_beddepth]) && any(ishandle([p_pkdepth(tmp1) p_beddepth])))
                    delete([p_pkdepth(tmp1) p_beddepth])
                    tmp2    = ind_decim(~isnan(ind_bed(ind_decim)) & ~isnan(ind_surf(ind_decim)));
                    p_beddepth ...
                            = plot(block.dist_lin(tmp2), (1e6 .* (block.twtt_bed(tmp2) - block.twtt_surf(tmp2) + block.twtt(1))), 'm.', 'markersize', 12, 'visible', 'off');
                end
            end
            if flat_done
                if (any([p_pkflat(tmp1) p_bedflat]) && any(ishandle([p_pkflat(tmp1) p_bedflat])))
                    delete([p_pkflat(tmp1) p_bedflat])
                    p_bedflat ...
                            = plot(block.dist_lin(ind_decim_flat(~isnan(ind_bed_flat(ind_decim_flat)))), (1e6 .* block.twtt(ind_bed_flat(ind_decim_flat(~isnan(ind_bed_flat(ind_decim_flat)))))), 'm.', 'markersize', 24, 'visible', 'off');
                end
            end
            
            if (logical(p_pksmooth(tmp1)) && ishandle(p_pksmooth(tmp1)))
                delete(p_pksmooth(tmp1))
            end
            if (logical(p_pksmoothdepth(tmp1)) && ishandle(p_pksmoothdepth(tmp1)))
                delete(p_pksmoothdepth(tmp1))
            end
            if (logical(p_pksmoothflat(tmp1)) && ishandle(p_pksmoothflat(tmp1)))
                delete(p_pksmoothflat(tmp1))
            end
            
            tmp3            = setdiff(1:pk.num_layer, tmp1);
            [pk.layer, smooth_done, p_pk, p_pkdepth, p_pkflat, p_pksmooth, p_pksmoothdepth, p_pksmoothflat, ind_y_flat_mean, ind_y_flat_smooth, pk.num_layer] ...
                            = deal(pk.layer(tmp3), smooth_done(tmp3), p_pk(tmp3), p_pkdepth(tmp3), p_pkflat(tmp3), p_pksmooth(tmp3), p_pksmoothdepth(tmp3), p_pksmoothflat(tmp3), ind_y_flat_mean(tmp3, :), ind_y_flat_smooth(tmp3, :), (pk.num_layer - 1));
            
            do_surfbed      = true;                        
            set(status_box, 'string', ['Bed and layer #' num2str(tmp1) ' merged.'])
            
            if (curr_layer > tmp1)
                curr_layer  = curr_layer - 1;
            end
            pk_sort
            curr_layer      = find((tmp1 == curr_layer), 1); % using tmp1 from pk_sort
            if (isempty(curr_layer) || any(curr_layer < 1))
                curr_layer  = 1;
            end
            
        else
            
            % replace NaN values in first layer with those in the second layer
            pk.layer(curr_layer).ind_y(isnan(pk.layer(curr_layer).ind_y)) ...
                            = pk.layer(tmp1).ind_y(isnan(pk.layer(curr_layer).ind_y));
            pk.layer(curr_layer).ind_y_smooth(isnan(pk.layer(curr_layer).ind_y_smooth)) ...
                            = pk.layer(tmp1).ind_y_smooth(isnan(pk.layer(curr_layer).ind_y_smooth));
            ind_y_flat_mean(curr_layer, isnan(ind_y_flat_mean(curr_layer, :))) ...
                            = ind_y_flat_mean(tmp1, isnan(ind_y_flat_mean(curr_layer, :)));
            ind_y_flat_smooth(curr_layer, isnan(ind_y_flat_smooth(curr_layer, :))) ...
                            = ind_y_flat_smooth(tmp1, isnan(ind_y_flat_smooth(curr_layer, :)));
            
            % fix plots
            tmp2            = ind_decim(~isnan(pk.layer(curr_layer).ind_y(ind_decim)));
            if (any(p_pk([curr_layer tmp1])) && any(ishandle(p_pk([curr_layer tmp1]))))
                delete(p_pk([curr_layer tmp1]))
                p_pk(curr_layer) ...
                            = plot(block.dist_lin(tmp2), (1e6 .* block.twtt(round(pk.layer(curr_layer).ind_y(tmp2)))), 'r.', 'markersize', 24, 'visible', 'off');
            end
            if depth_avail
                if (any(p_pkdepth([curr_layer tmp1])) && any(ishandle(p_pkdepth([curr_layer tmp1]))))
                    delete(p_pkdepth([curr_layer tmp1]))
                    tmp2    = ind_decim(~isnan(pk.layer(curr_layer).ind_y(ind_decim)) & ~isnan(ind_surf(ind_decim)));
                    tmp3    = pk.layer(curr_layer).ind_y(tmp2) - ind_surf(tmp2) + 1;
                    tmp3((tmp3 < 1) | (tmp3 > num_sample_trim)) ...
                            = NaN;
                    p_pkdepth(curr_layer) ...
                            = plot(block.dist_lin(tmp2(~isnan(tmp3))), (1e6 .* block.twtt(round(tmp3(~isnan(tmp3))))), 'r.', 'markersize', 24, 'visible', 'off');
                end
            end
            if flat_done
                if (any(p_pkflat([curr_layer tmp1])) && any(ishandle(p_pkflat([curr_layer tmp1]))))
                    delete(p_pkflat([curr_layer tmp1]))
                    p_pkflat(curr_layer) ...
                            = plot(block.dist_lin(ind_decim_flat(~isnan(ind_y_flat_mean(curr_layer, :)))), (1e6 .* block.twtt(ind_y_flat_mean(curr_layer, ~isnan(ind_y_flat_mean(curr_layer, :))))), 'r.', 'markersize', 24, 'visible', 'off');
                end
            end
            
            if (logical(p_pksmooth(curr_layer)) && ishandle(p_pksmooth(curr_layer)))
                delete(p_pksmooth(curr_layer))
            end
            if (logical(p_pksmooth(tmp1)) && ishandle(p_pksmooth(tmp1)))
                delete(p_pksmooth(tmp1))
            end
            if (logical(p_pksmoothdepth(curr_layer)) && ishandle(p_pksmoothdepth(curr_layer)))
                delete(p_pksmoothdepth(curr_layer))
            end
            if (logical(p_pksmoothdepth(tmp1)) && ishandle(p_pksmoothdepth(tmp1)))
                delete(p_pksmoothdepth(tmp1))
            end
            if (logical(p_pksmoothflat(curr_layer)) && ishandle(p_pksmoothflat(curr_layer)))
                delete(p_pksmoothflat(curr_layer))
            end
            if (logical(p_pksmoothflat(tmp1)) && ishandle(p_pksmoothflat(tmp1)))
                delete(p_pksmoothflat(tmp1))
            end
            
            smooth_done(curr_layer) ...
                            = false;
            tmp3            = setdiff(1:pk.num_layer, tmp1);
            [pk.layer, smooth_done, p_pk, p_pkdepth, p_pkflat, p_pksmooth, p_pksmoothdepth, p_pksmoothflat, ind_y_flat_mean, ind_y_flat_smooth, pk.num_layer] ...
                            = deal(pk.layer(tmp3), smooth_done(tmp3), p_pk(tmp3), p_pkdepth(tmp3), p_pkflat(tmp3), p_pksmooth(tmp3), p_pksmoothdepth(tmp3), p_pksmoothflat(tmp3), ind_y_flat_mean(tmp3, :), ind_y_flat_smooth(tmp3, :), (pk.num_layer - 1));
            
            set(status_box, 'string', ['Layers #' num2str(curr_layer) ' and #' num2str(tmp1) ' merged.'])
            
            if (curr_layer > tmp1)
                curr_layer  = curr_layer - 1;
            end
            pk_sort
            curr_layer      = find((tmp1 == curr_layer), 1); % using tmp1 from pk_sort
            if (isempty(curr_layer) || any(curr_layer < 1))
                curr_layer  = 1;
            end
            
            pk_smooth
        end
        
        set(layer_list, 'string', [num2cell(1:pk.num_layer) 'surface' 'bed'], 'value', curr_layer)
        match_done          = false;
        pk_select
        show_surfbed
        show_pk
        show_smooth
    end

%% Smooth picked layers

    function pk_smooth(source, eventdata)
        
        if ~pk_done
            set(status_box, 'string', 'No picked layers to smooth yet.')
            return
        end
        if all(smooth_done)
            set(status_box, 'string', 'No layers need smoothing.')
            return
        end
        
        set(status_box, 'string', 'Smoothing picked layers...')
        
        tmp1                = [];
        warning('off', 'MATLAB:interp1:NaNinY')
        for ii = find(~smooth_done)
            pk.layer(ii).ind_y_smooth ...
                            = round(smooth_lowess(pk.layer(ii).ind_y, round(pk.length_smooth / nanmean(diff(block.dist_lin))))');
            pk.layer(ii).ind_y_smooth((pk.layer(ii).ind_y_smooth < 1) | (pk.layer(ii).ind_y_smooth > num_sample_trim)) ...
                            = NaN;
            if flat_done
                ind_y_flat_smooth(ii, :) ...
                            = round(smooth_lowess(ind_y_flat_mean(ii, :), round(pk.length_smooth / (nanmean(diff(block.dist_lin)) * decim_flat)))');
                ind_y_flat_smooth(ii, ((ind_y_flat_smooth(ii, :) < 1) | (ind_y_flat_smooth(ii, :) > num_sample_trim))) ...
                            = NaN;
            end
            if all(isnan(pk.layer(ii).ind_y_smooth))
                tmp1        = [tmp1 ii]; %#ok<AGROW>
            end
            if flat_done
                if all(isnan(ind_y_flat_smooth(ii, :)))
                    tmp1    = [tmp1 ii]; %#ok<AGROW>
                end
            end
        end
        warning('on', 'MATLAB:interp1:NaNinY')
        
        % remove layers that are empty when smoothed
        axes(ax_radar)
        if ~isempty(tmp1)
            for ii = tmp1
                if (logical(p_pk(ii)) && ishandle(p_pk(ii)))
                    delete(p_pk(ii))
                end
                if (logical(p_pkdepth(ii)) && ishandle(p_pkdepth(ii)))
                    delete(p_pkdepth(ii))
                end
                if (logical(p_pkflat(ii)) && ishandle(p_pkflat(ii)))
                    delete(p_pkflat(tmp1))
                end
                if (logical(p_pksmooth(ii)) && ishandle(p_pksmooth(ii)))
                    delete(p_pksmooth(ii))
                end
                if (logical(p_pksmoothdepth(ii)) && ishandle(p_pksmoothdepth(ii)))
                    delete(p_pksmoothdepth(ii))
                end
                if (logical(p_pksmoothflat(ii)) && ishandle(p_pksmoothflat(ii)))
                    delete(p_pksmoothflat(ii))
                end
            end
            tmp2            = setdiff(1:pk.num_layer, tmp1);
            if (length(tmp2) > 1)
                curr_layer  = interp1(1:length(tmp2), tmp2, curr_layer, 'nearest', 'extrap');
            else
                curr_layer  = tmp2;
            end
            [pk.num_layer, pk.layer, p_pk, p_pkdepth, smooth_done, p_pksmooth, p_pksmoothdepth, p_pkflat, p_pksmoothflat] ...
                            = deal(length(tmp2), pk.layer(tmp2), p_pk(tmp2), p_pkdepth(tmp2), smooth_done(tmp2), p_pksmooth(tmp2), p_pksmoothdepth(tmp2), p_pkflat(tmp2), p_pksmoothflat(tmp2));
            set(layer_list, 'string', [num2cell(1:pk.num_layer) 'surface' 'bed'], 'value', curr_layer)
            set(status_box, 'string', ['Layer(s) #' num2str(tmp1) ' deleted because they were empty after smoothing.'])
        end
        
        for ii = find(~smooth_done)
            if (logical(p_pksmooth(ii)) && ishandle(p_pksmooth(ii)))
                delete(p_pksmooth(ii))
            end
            if (logical(p_pksmoothdepth(ii)) && ishandle(p_pksmoothdepth(ii)))
                delete(p_pksmoothdepth(ii))
            end
            if (logical(p_pksmoothflat(ii)) && ishandle(p_pksmoothflat(ii)))
                delete(p_pksmoothflat(ii))
            end
            tmp1            = pk.layer(ii).ind_y_smooth(ind_decim);
            tmp2            = block.dist_lin(ind_decim);
            p_pksmooth(ii)  = plot(tmp2(~isnan(tmp1)), (1e6 .* block.twtt(round(tmp1(~isnan(tmp1))))), 'g.', 'markersize', 12, 'visible', 'off');
            if depth_avail
                tmp2        = ind_decim(~isnan(pk.layer(ii).ind_y_smooth(ind_decim)) & ~isnan(ind_surf(ind_decim)));
                tmp3        = pk.layer(ii).ind_y_smooth(tmp2) - ind_surf(tmp2) + 1;
                tmp3((tmp3 < 1) | (tmp3 > num_sample_trim)) ...
                            = NaN;
                p_pksmoothdepth(ii) ...
                            = plot(block.dist_lin(tmp2(~isnan(tmp3))), (1e6 .* block.twtt(round(tmp3(~isnan(tmp3))))), 'g.', 'markersize', 12, 'visible', 'off');
            end
            if (flat_done && any(~isnan(ind_y_flat_smooth(ii, :))))
                p_pksmoothflat(ii) ...
                            = plot(block.dist_lin(ind_decim_flat(~isnan(ind_y_flat_smooth(ii, :)))), (1e6 .* block.twtt(round(ind_y_flat_smooth(ii, ~isnan(ind_y_flat_smooth(ii, :)))))), 'g.', 'markersize', 12, 'visible', 'off');
            end
        end
        set(smooth_check, 'value', 1)
        set(status_box, 'string', ['Smoothed ' num2str(length(find(~smooth_done))) ' layer(s).'])
        smooth_done(~smooth_done) ...
                            = true;
        match_done          = false;
        show_smooth
    end

%% Match picked layers in current block to reference layers

    function pk_match(source, eventdata)
        
        if ~pk.num_layer
            set(status_box, 'string', 'No picked layers to match yet.')
            return
        end
        if ~all(smooth_done)
            set(status_box, 'string', 'Not all layers have been smoothed yet.')
            return
        end
        if ~strcmp(disp_type, 'twtt')
            set(status_box, 'string', 'Layers must be matched in twtt space. Switching...')
            disp_type       = 'twtt';
            set(disp_group, 'selectedobject', disp_check(1))
            plot_twtt
        end
        if ~ref_done
            set(status_box, 'string', 'First picked block of the sub-transect? Y: yes; otherwise: no.')
            waitforbuttonpress
            if strcmpi(get(pkgui, 'currentcharacter'), 'Y')
                pk.ind_match= (1:pk.num_layer)';
                match_done  = true;
                set(status_box, 'string', 'Layer numbers set.')
                pk.ind_match_max ...
                            = pk.num_layer;
                return
            else
                set(status_box, 'string', 'No reference layers loaded to match to.')
                return
            end
        end
        if strcmp(ref_start_or_end, 'end')
            set(status_box, 'string', 'Cannot match to right-hand-side reference layers.')
            return
        end
        set(status_box, 'string', 'Matching layers...')
        
        % everything starts NaN
        pk.ind_match        = NaN(pk.num_layer, 1);
        
        % get mean indices of reference layers and their standard deviation
        tmp1                = NaN(pk_ref.num_layer, block.ind_overlap(1));
        for ii = 1:pk_ref.num_layer
            tmp1(ii, :)     = pk_ref.layer(ii).twtt_smooth(pk_ref.ind_overlap(2):end);
        end
        
        % fade layer colors
        set(p_pk(logical(p_pk) & ishandle(p_pk)), 'color', [1 0.7 0.7])
        set(p_pksmooth(logical(p_pksmooth) & ishandle(p_pksmooth)), 'color', [0.7 1 0.7])
        set(p_ref(logical(p_ref) & ishandle(p_ref)), 'color', [1 1 0.7])
        if depth_avail
            set(p_pkdepth(logical(p_pkdepth) & ishandle(p_pkdepth)), 'color', [1 0.7 0.7])
            set(p_pksmoothdepth(logical(p_pksmoothdepth) & ishandle(p_pksmoothdepth)), 'color', [0.7 1 0.7])
            set(p_refdepth(logical(p_refdepth) & ishandle(p_refdepth)), 'color', [1 1 0.7])
        end
        if flat_done
            set(p_pkflat(logical(p_pkflat) & ishandle(p_pkflat)), 'color', [1 0.7 0.7])
            set(p_pksmoothflat(logical(p_pksmoothflat) & ishandle(p_pksmoothflat)), 'color', [0.7 1 0.7])
        end
        
        % search for matching layers near reference layers
        for ii = 1:pk.num_layer
            
            tmp2            = NaN(pk_ref.num_layer, block.ind_overlap(1));
            tmp2(:, ~isnan(pk.layer(ii).ind_y_smooth(1:block.ind_overlap(1)))) ...
                            = tmp1(:, ~isnan(pk.layer(ii).ind_y_smooth(1:block.ind_overlap(1)))) - repmat(block.twtt(round(pk.layer(ii).ind_y_smooth(~isnan(pk.layer(ii).ind_y_smooth(1:block.ind_overlap(1))))))', pk_ref.num_layer, 1);
            tmp3            = NaN(1, pk_ref.num_layer);
            for jj = 1:pk_ref.num_layer
                tmp3(jj)    = abs(nanmean(tmp2(jj, :)));
            end
            
            if (length(find(tmp3 < pk.twtt_match)) == 1) % found a single matching layer within the threshold
                pk.ind_match(ii) ...
                    = pk_ref.ind_match(logical(tmp3 < pk.twtt_match));
                % brighten layer colors once successful
                if (logical(p_pk(ii)) && ishandle(p_pk(ii)))
                    set(p_pk(ii), 'color', 'r')
                end
                if (logical(p_pksmooth(ii)) && ishandle(p_pksmooth(ii)))
                    set(p_pksmooth(ii), 'color', 'g')
                end
                if (logical(p_ref(logical(tmp3 < pk.twtt_match))) && ishandle(p_ref(logical(tmp3 < pk.twtt_match))))
                    set(p_ref(logical(tmp3 < pk.twtt_match)), 'color', 'y')
                end
                if depth_avail
                    if (logical(p_pkdepth(ii)) && ishandle(p_pkdepth(ii)))
                        set(p_pkdepth(ii), 'color', 'r')
                    end
                    if (logical(p_pksmoothdepth(ii)) && ishandle(p_pksmoothdepth(ii)))
                        set(p_pksmoothdepth(ii), 'color', 'g')
                    end
                    if (logical(p_refdepth(logical(tmp3 < pk.twtt_match))) && ishandle(p_refdepth(logical(tmp3 < pk.twtt_match))))
                        set(p_refdepth(logical(tmp3 < pk.twtt_match)), 'color', 'y')
                    end
                end
                if flat_done
                    if (logical(p_pkflat(ii)) && ishandle(p_pkflat(ii)))
                        set(p_pkflat(ii), 'color', 'r')
                    end
                    if (logical(p_pksmoothflat(ii)) && ishandle(p_pksmoothflat(ii)))
                        set(p_pksmoothflat(ii), 'color', 'g')
                    end
                end
            end
        end
        
        set(status_box, 'string', ['Matched ' num2str(length(find(~isnan(pk.ind_match)))) '/' num2str(pk.num_layer) ' new layers from ' num2str(length(find(sum(~isnan(tmp1), 2)))) ' overlapping reference layers.'])
        
        % add new layer numbers for those that didn't match any reference layers
        if any(isnan(pk.ind_match))
            tmp1            = nanmax(pk.ind_match);
            if isfield(pk_ref, 'ind_match_max') % if available, take care not to repeat layer numbers
                if ((pk_ref.ind_match_max > tmp1) || isnan(tmp1))
                    tmp1    = pk_ref.ind_match_max;
                end
            end
            pk.ind_match(isnan(pk.ind_match)) ...
                            = (tmp1 + 1):(tmp1 + length(find(isnan(pk.ind_match))));
        end
        
        pk.ind_match_max    = max(pk.ind_match);
        
        match_done          = true;
    end

%% Assign current layer to surface or bed

    function pk_surfbed(source, eventdata)
        
        if ~(pk_done && pk.num_layer && curr_layer)
            set(status_box, 'string', 'No picked layers to assign to surface or bed.')
            return
        end
        set(status_box, 'string', 'Assign current layer? S: surface; B: bed; otherwise: no.')
        [~, ~, button]      = ginput(1);
        
        if strcmpi(char(button), 'S')
            
            ind_surf        = round(pk.layer(curr_layer).ind_y_smooth);
            block.twtt_surf = NaN(1, block.num_trace);
            block.twtt_surf(~isnan(pk.layer(curr_layer).ind_y_smooth)) ...
                            = block.twtt(round(pk.layer(curr_layer).ind_y_smooth(~isnan(pk.layer(curr_layer).ind_y_smooth))));
            pk_del_breakout
            surf_avail      = true;
            
            if (logical(p_surf) && ishandle(p_surf))
                delete(p_surf)
            end
            if (logical(p_surfflat) && ishandle(p_surfflat))
                delete(p_surfflat)
            end
            p_surf          = plot(block.dist_lin(ind_decim), (1e6 .* block.twtt_surf(ind_decim)), 'm.', 'markersize', 12, 'visible', 'off');
            if flat_done
                if (surf_avail  && ~all(isnan(ind_surf_flat(ind_decim))))
                    tmp1    = ind_surf_flat(ind_decim_flat);
                    p_surfflat ...
                            = plot(block.dist_lin(ind_decim_flat(~isnan(tmp1))), (1e6 .* block.twtt(tmp1(~isnan(tmp1)))), 'm.', 'markersize', 12, 'visible', 'off');
                end
            end
            
            if ~isempty(find(~isnan(ind_surf(ind_decim)), 1))
                if (logical(p_beddepth) && ishandle(p_beddepth))
                    delete(p_beddepth)
                end
                if (any(p_pkdepth) && any(ishandle(p_pkdepth)))
                    delete(p_pkdepth(logical(p_pkdepth) & ishandle(p_pkdepth)))
                end
                if (any(p_pksmoothdepth) && any(ishandle(p_pksmoothdepth)))
                    delete(p_pksmoothdepth(logical(p_pksmoothdepth) & ishandle(p_pksmoothdepth)))
                end
                amp_depth   = NaN(num_sample_trim, num_decim, 'single');
                for ii = find(~isnan(ind_surf(ind_decim)))
                    amp_depth(1:(num_sample_trim - ind_surf(ind_decim(ii)) + 1), ii) ...
                            = amp_mean(ind_surf(ind_decim(ii)):num_sample_trim, ii); % shift data up to surface
                end
                set(disp_check(2), 'visible', 'on')
                depth_avail = true;
                [p_pkdepth, p_pksmoothdepth] ...
                            = deal(zeros(1, pk.num_layer));
                for ii = 1:pk.num_layer
                    if ~isempty(find((~isnan(pk.layer(ii).ind_y(ind_decim)) & ~isnan(ind_surf(ind_decim))), 1))
                        tmp1= ind_decim(~isnan(pk.layer(ii).ind_y(ind_decim)) & ~isnan(ind_surf(ind_decim)));
                        tmp2= pk.layer(ii).ind_y(tmp1) - ind_surf(tmp1) + 1;
                        tmp2((tmp2 < 1) | (tmp2 > num_sample_trim)) ...
                            = NaN;
                        p_pkdepth(ii) ...
                            = plot(block.dist_lin(tmp1(~isnan(tmp2))), (1e6 .* block.twtt(round(tmp2(~isnan(tmp2))))), 'r.', 'markersize', 12, 'visible', 'off');
                    end
                    if ~isempty(find((~isnan(pk.layer(ii).ind_y_smooth(ind_decim)) & ~isnan(ind_surf(ind_decim))), 1))
                        tmp1= ind_decim(~isnan(pk.layer(ii).ind_y_smooth(ind_decim)) & ~isnan(ind_surf(ind_decim)));
                        tmp2= pk.layer(ii).ind_y_smooth(tmp1) - ind_surf(tmp1) + 1;
                        tmp2((tmp2 < 1) | (tmp2 > num_sample_trim)) ...
                            = NaN;
                        p_pksmoothdepth(ii) ...
                            = plot(block.dist_lin(tmp1(~isnan(tmp2))), (1e6 .* block.twtt(round(tmp2(~isnan(tmp2))))), 'g.', 'markersize', 12, 'visible', 'off');
                    end
                end
                if bed_avail
                    tmp1    = ind_decim(~isnan(ind_bed(ind_decim)) & ~isnan(ind_surf(ind_decim)));
                    p_beddepth ...
                            = plot(block.dist_lin(tmp1), (1e6 .* (block.twtt_bed(tmp1) - block.twtt_surf(tmp1) + block.twtt(1))), 'm.', 'markersize', 12, 'visible', 'off');
                end
                if strcmp(disp_type, '~depth')
                    plot_depth
                end
            end
            
            do_surfbed      = true;
            curr_layer      = pk.num_layer + 1;
            set(layer_list, 'string', [num2cell(1:pk.num_layer) 'surface' 'bed'], 'value', (pk.num_layer + 1))
            show_surfbed
            set(status_box, 'string', 'Layer assigned to surface.')
            
        elseif strcmpi(char(button), 'B')
            
            ind_bed         = round(pk.layer(curr_layer).ind_y_smooth);
            block.twtt_bed  = NaN(1, block.num_trace);
            block.twtt_bed(~isnan(pk.layer(curr_layer).ind_y_smooth)) ...
                            = block.twtt(round(pk.layer(curr_layer).ind_y_smooth(~isnan(pk.layer(curr_layer).ind_y_smooth))));
            pk_del_breakout
            bed_avail       = true;
            if surf_avail
                depth_avail = true;
            end
            
            if (logical(p_bed) && ishandle(p_bed))
                delete(p_bed)
            end
            if (logical(p_beddepth) && ishandle(p_beddepth))
                delete(p_beddepth)
            end            
            if (logical(p_bedflat) && ishandle(p_bedflat))
                delete(p_bedflat)
            end
            p_bed           = plot(block.dist_lin(ind_decim), (1e6 .* block.twtt_bed(ind_decim)), 'm.', 'markersize', 12, 'visible', 'off');
            if (bed_avail && depth_avail)
                tmp1        = ind_decim(~isnan(ind_bed(ind_decim)) & ~isnan(ind_surf(ind_decim)));
                p_beddepth  = plot(block.dist_lin(tmp1), (1e6 .* (block.twtt_bed(tmp1) - block.twtt_surf(tmp1) + block.twtt(1))), 'm.', 'markersize', 12, 'visible', 'off');
            end
            if flat_done
                if (surf_avail  && ~all(isnan(ind_surf_flat(ind_decim))))
                    tmp1    = ind_bed_flat(ind_decim_flat);
                    p_bedflat ...
                            = plot(block.dist_lin(ind_decim_flat(~isnan(tmp1))), (1e6 .* block.twtt(tmp1(~isnan(tmp1)))), 'm.', 'markersize', 12, 'visible', 'off');
                end
            end
            
            do_surfbed      = true;
            set(layer_list, 'string', [num2cell(1:pk.num_layer) 'surface' 'bed'], 'value', (pk.num_layer + 2))
            show_surfbed
            set(status_box, 'string', 'Layer assigned to bed.')
            
        else
            set(status_box, 'string', 'Layer not assigned.')
        end
    end

%% Save layer picks

    function pk_save(source, eventdata)
        
        % want everything done before saving
        if ~match_done
            set(status_box, 'string', 'Layers not matched yet, or need re-matching.')
            return
        end
        
        set(status_box, 'string', 'Saving picks...')
        pause(0.1)
        
        reset_xy
        pause(0.1)
        
        if isempty(path_pk)
            if ispc
                if (exist([path_data '..\pk'], 'dir') || exist([path_data '..\..\pk'], 'dir'))
                    path_pk = [path_data(1:strfind(path_data, '\block')) 'pk\'];
                end
            else
                if (exist([path_data '../pk'], 'dir') || exist([path_data '../../pk'], 'dir'))
                    path_pk = [path_data(1:strfind(path_data, '/block')) 'pk/'];
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
            [file_pk, path_pk] = deal('', tmp1);
            set(status_box, 'string', 'Saving cancelled.')
            return
        end
            
        set(status_box, 'string', 'Saving picks...')
        
        tmp1                = pk;
        
        % save many variables in pk structure for easy reference independent of data later on
        [pk.lat, pk.lon, pk.x, pk.y, pk.num_sample, pk.num_trace, pk.file_in, pk.file_block, pk.twtt_min_ref, pk.twtt_max_ref, pk.dist, pk.block.dist_lin, pk.ind_overlap, pk.elev_air, pk.time] ...
                            = deal(block.lat, block.lon, block.x, block.y, block.num_sample, block.num_trace, block.file_in, file_data(1:(end - 4)), twtt_min_ref, twtt_max_ref, block.dist, block.dist_lin, block.ind_overlap, block.elev_air, block.time);
        if isfield(block, 'elev_air_gimp')
            pk.elev_air_gimp= block.elev_air_gimp;
        else
            pk.elev_air_gimp= NaN(1, block.num_trace);
        end
        pk.twtt_surf        = block.twtt_surf;
        pk.twtt_bed         = block.twtt_bed;
        
        % get traveltimes and echo intensities from indices, and adjust indices as appropriate assuming trimming has occurred
        pk.elev_surf        = block.elev_air - (block.twtt_surf .* (speed_vacuum / 2)); % ice-sheet surface elevation (used to calculate layer elevations)
        pk.elev_surf_gimp   = pk.elev_air_gimp - (block.twtt_surf .* (speed_vacuum / 2));
        pk.elev_bed         = pk.elev_surf - ((block.twtt_bed - block.twtt_surf) .* (speed_ice / 2)); % bed elevation
        pk.elev_bed_gimp    = pk.elev_surf_gimp - ((block.twtt_bed - block.twtt_surf) .* (speed_ice / 2));
        
        [pk.int_surf, pk.int_bed] ...
                            = deal(NaN(1, block.num_trace));
        if surf_avail
            pk.int_surf(~isnan(ind_surf)) ...
                            = block.amp(sub2ind([num_sample_trim block.num_trace], ind_surf(~isnan(ind_surf)), find(~isnan(ind_surf))));
        end
        if bed_avail
            pk.int_bed(~isnan(ind_bed)) ...
                            = block.amp(sub2ind([num_sample_trim block.num_trace], ind_bed(~isnan(ind_bed)), find(~isnan(ind_bed))));
        end
        
        for ii = 1:pk.num_layer
            [pk.layer(ii).twtt, pk.layer(ii).twtt_smooth, pk.layer(ii).int, pk.layer(ii).int_smooth, pk.layer(ii).depth, pk.layer(ii).depth_smooth, pk.layer(ii).elev, pk.layer(ii).elev_smooth] ...
                            = deal(NaN(1, block.num_trace));
            [pk.layer(ii).twtt(~isnan(pk.layer(ii).ind_y)), pk.layer(ii).twtt_smooth(~isnan(pk.layer(ii).ind_y_smooth))] ...
                            = deal(block.twtt(round(pk.layer(ii).ind_y(~isnan(pk.layer(ii).ind_y))))', block.twtt(round(pk.layer(ii).ind_y_smooth(~isnan(pk.layer(ii).ind_y_smooth))))');
            [pk.layer(ii).twtt_ice, pk.layer(ii).twtt_ice_smooth] ...
                            = deal((pk.layer(ii).twtt - block.twtt_surf), (pk.layer(ii).twtt_smooth - block.twtt_surf));
            [pk.layer(ii).depth, pk.layer(ii).depth_smooth] ...
                            = deal((pk.layer(ii).twtt_ice .* (speed_ice / 2)), (pk.layer(ii).twtt_ice_smooth .* (speed_ice / 2)));
            [pk.layer(ii).elev, pk.layer(ii).elev_smooth] ...
                            = deal((pk.elev_surf - pk.layer(ii).depth), (pk.elev_surf - pk.layer(ii).depth_smooth));
            [pk.layer(ii).elev_gimp, pk.layer(ii).elev_smooth_gimp] ...
                            = deal((pk.elev_surf_gimp - pk.layer(ii).depth), (pk.elev_surf_gimp - pk.layer(ii).depth_smooth));
            [pk.layer(ii).int(~isnan(pk.layer(ii).ind_y)), pk.layer(ii).int_smooth(~isnan(pk.layer(ii).ind_y_smooth))] ...
                            = deal(block.amp(sub2ind([num_sample_trim block.num_trace], round(pk.layer(ii).ind_y(~isnan(pk.layer(ii).ind_y))), find(~isnan(pk.layer(ii).ind_y)))), ...
                                   block.amp(sub2ind([num_sample_trim block.num_trace], round(pk.layer(ii).ind_y_smooth(~isnan(pk.layer(ii).ind_y_smooth))), find(~isnan(pk.layer(ii).ind_y_smooth)))));
            [pk.layer(ii).ind_y, pk.layer(ii).ind_y_smooth] ...
                            = deal((pk.layer(ii).ind_y + pk.ind_trim_start), (pk.layer(ii).ind_y_smooth + pk.ind_trim_start)); % adjust due to trimming
        end
        
        pk                  = orderfields(pk);
        pk.layer            = orderfields(pk.layer);
        
        save([path_pk file_pk], '-v7.3', 'pk')
        
        % make a simple figure that also gets saved
        set(0, 'DefaultFigureWindowStyle', 'default')
        pkfig               = figure('position', [10 10 1600 1000]);
        imagesc(block.dist_lin((1 + decim):decim:(block.num_trace - decim)), (1e6 .* block.twtt), amp_mean, [db_min db_max])
        hold on
        colormap(bone)
        for ii = 1:pk.num_layer
            plot(block.dist_lin(ind_decim), (1e6 .* pk.layer(ii).twtt_smooth(ind_decim)), 'g.', 'markersize', 12)
        end
        set(gca, 'fontsize', 16)
        xlabel('Distance (km)')
        ylabel('Traveltime ({\mu}s)')
        title(file_pk(1:(end - 4)), 'fontweight', 'bold', 'interpreter', 'none')
        grid on
        box on
        
        pause(0.5)
        print(pkfig, '-dpng', [path_pk file_pk(1:(end - 4)) '.png'])
        
        pk                  = tmp1;
        set(0, 'DefaultFigureWindowStyle', 'docked')
        
        set(status_box, 'string', ['Picks saved as ' file_pk(1:(end - 4)) '.'])
        
        if do_surfbed
            set(status_box, 'string', 'Surface/bed edited so re-saving block...')
            pause(0.1)
            if trim_done
                tmp1        = load([path_data file_data], 'block');
                tmp2        = block;
                [block.amp, block.twtt] ...
                            = deal(tmp1.block.amp, tmp1.block.twtt);
                if phase_avail
                    block.phase_diff_filt ...
                            = tmp1.block.phase_diff_filt;
                end
                if aresp_avail
                    block.slope_aresp ...
                            = tmp1.block.slope_aresp;
                end
                tmp1        = 0;
            end
            save([path_data file_data], '-v7.3', 'block')
            pause(0.1)
            if trim_done
                block       = tmp2;
                tmp2        = 0;
            end
            set(status_box, 'string', ['Block saved as ' file_data(1:(end - 4)) ' in ' path_data '.'])
        end
        
    end

%% Update minimum twtt

    function slide_twtt_min(source, eventdata)
        if (((1e6 * twtt_max_ref) - (get(twtt_min_slide, 'value') - (1e6 * twtt_min_ref))) < (1e6 * twtt_max))
            if get(twttfix_check, 'value')
                tmp1        = twtt_max - twtt_min;
            end
            twtt_min        = ((1e6 * twtt_max_ref) - (get(twtt_min_slide, 'value') - (1e6 * twtt_min_ref))) * 1e-6;
            if get(twttfix_check, 'value')
                twtt_max    = twtt_min + tmp1;
                if (twtt_max > twtt_max_ref)
                    twtt_max= twtt_max_ref;
                    twtt_min= twtt_max - tmp1;
                    if (twtt_min < twtt_min_ref)
                        twtt_min ...
                            = twtt_min_ref;
                    end
                    if ((1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref))) < get(twtt_min_slide, 'min'))
                        set(twtt_min_slide, 'value', get(twtt_min_slide, 'min'))
                    elseif ((1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref))) > get(twtt_min_slide, 'max'))
                        set(twtt_min_slide, 'value', get(twtt_min_slide, 'max'))
                    else
                        set(twtt_min_slide, 'value', (1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref))))
                    end
                end
                set(twtt_max_edit, 'string', sprintf('%3.1f', (1e6 * twtt_max)))
                if ((1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref))) < get(twtt_max_slide, 'min'))
                    set(twtt_max_slide, 'value', get(twtt_max_slide, 'min'))
                elseif ((1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref))) > get(twtt_max_slide, 'max'))
                    set(twtt_max_slide, 'value', get(twtt_max_slide, 'max'))
                else
                    set(twtt_max_slide, 'value', (1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref))))
                end
            end
            set(twtt_min_edit, 'string', sprintf('%3.1f', (1e6 * twtt_min)))
            update_twtt_range
        else
            if ((1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref))) < get(twtt_min_slide, 'min'))
                set(twtt_min_slide, 'value', get(twtt_min_slide, 'min'))
            elseif ((1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref))) > get(twtt_min_slide, 'max'))
                set(twtt_min_slide, 'value', get(twtt_min_slide, 'max'))
            else
                set(twtt_min_slide, 'value', (1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref))))
            end
        end
        set(twtt_min_slide, 'enable', 'off')
        drawnow
        set(twtt_min_slide, 'enable', 'on')
    end

%% Update maximum twtt

    function slide_twtt_max(source, eventdata)
        if (((1e6 * twtt_max_ref) - (get(twtt_max_slide, 'value') - (1e6 * twtt_min_ref))) > (1e6 * twtt_min))
            if get(twttfix_check, 'value')
                tmp1        = twtt_max - twtt_min;
            end
            twtt_max        = ((1e6 * twtt_max_ref) - (get(twtt_max_slide, 'value') - (1e6 * twtt_min_ref))) * 1e-6;
            if get(twttfix_check, 'value')
                twtt_min    = twtt_max - tmp1;
                if (twtt_min < twtt_min_ref)
                    twtt_min= twtt_min_ref;
                    twtt_max= twtt_min + tmp1;
                    if (twtt_max > twtt_max_ref)
                        twtt_max ...
                            = twtt_max_ref;
                    end
                    if ((1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref))) < get(twtt_max_slide, 'min'))
                        set(twtt_max_slide, 'value', get(twtt_max_slide, 'min'))
                    elseif ((1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref))) > get(twtt_max_slide, 'max'))
                        set(twtt_max_slide, 'value', get(twtt_max_slide, 'max'))
                    else
                        set(twtt_max_slide, 'value', (1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref))))
                    end
                end
                set(twtt_min_edit, 'string', sprintf('%3.1f', (1e6 * twtt_min)))
                if ((1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref))) < get(twtt_min_slide, 'min'))
                    set(twtt_min_slide, 'value', get(twtt_min_slide, 'min'))
                elseif ((1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref))) > get(twtt_min_slide, 'max'))
                    set(twtt_min_slide, 'value', get(twtt_min_slide, 'max'))
                else
                    set(twtt_min_slide, 'value', (1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref))))
                end
            end
            set(twtt_max_edit, 'string', sprintf('%3.1f', (1e6 * twtt_max)))
            update_twtt_range
        else
            if ((1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref))) < get(twtt_max_slide, 'min'))
                set(twtt_max_slide, 'value', get(twtt_max_slide, 'min'))
            elseif ((1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref))) > get(twtt_max_slide, 'max'))
                set(twtt_max_slide, 'value', get(twtt_max_slide, 'max'))
            else
                set(twtt_max_slide, 'value', (1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref))))
            end
        end
        set(twtt_max_slide, 'enable', 'off')
        drawnow
        set(twtt_max_slide, 'enable', 'on')
    end

%% Reset minimum twtt

    function reset_twtt_min(source, eventdata)
        if ((1e6 * twtt_max_ref) < get(twtt_min_slide, 'min'))
            set(twtt_min_slide, 'value', get(twtt_min_slide, 'min'))
        elseif ((1e6 * twtt_max_ref) > get(twtt_min_slide, 'max'))
            set(twtt_min_slide, 'value', get(twtt_min_slide, 'max'))
        else
            set(twtt_min_slide, 'value', (1e6 * twtt_max_ref))
        end
        set(twtt_min_edit, 'string', sprintf('%3.1f', (1e6 * twtt_min_ref)))
        twtt_min              = twtt_min_ref;
        update_twtt_range
    end

%% Reset maximum twtt

    function reset_twtt_max(source, eventdata)
        if ((1e6 * twtt_min_ref) < get(twtt_max_slide, 'min'))
            set(twtt_max_slide, 'value', get(twtt_max_slide, 'min'))
        elseif ((1e6 * twtt_min_ref) > get(twtt_max_slide, 'max'))
            set(twtt_max_slide, 'value', get(twtt_max_slide, 'max'))
        else
            set(twtt_max_slide, 'value', (1e6 * twtt_min_ref))
        end
        set(twtt_max_edit, 'string', sprintf('%3.1f', (1e6 * twtt_max_ref)))
        twtt_max              = twtt_max_ref;
        update_twtt_range
    end

%% Update twtt range

    function update_twtt_range(source, eventdata)
        axes(ax_radar)
        ylim(1e6 .* [twtt_min twtt_max])
        narrow_cb
    end

%% Update minimum dB/phase/ARESP

    function slide_db_min(source, eventdata)
        switch disp_type
            case {'twtt' '~depth' 'flat'}
                tmp1        = [db_min db_max];
                tmp2        = [db_min_ref db_max_ref];
            case 'phase'
                tmp1        = [phase_diff_min phase_diff_max];
                tmp2        = [phase_diff_min_ref phase_diff_max_ref];
            case 'ARESP'
                tmp1        = [aresp_min aresp_max];
                tmp2        = [aresp_min_ref aresp_max_ref];
        end
        if (get(cb_min_slide, 'value') < tmp1(2))
            if get(cbfix_check1, 'value')
                tmp3        = diff(tmp1);
            end
            tmp1(1)         = get(cb_min_slide, 'value');
            if get(cbfix_check1, 'value')
                tmp1(2)     = tmp1(1) + tmp3;
                if (tmp1(2) > tmp2(2))
                    tmp1(2) = tmp2(2);
                    tmp1(1) = tmp1(2) - tmp3;
                    if (tmp1(1) < tmp2(1))
                        tmp1(1) ...
                            = tmp2(1);
                    end
                    if (tmp1(1) < get(cb_min_slide, 'min'))
                        set(cb_min_slide, 'value', get(cb_min_slide, 'min'))
                    else
                        set(cb_min_slide, 'value', tmp1(1))
                    end
                end
                set(cb_max_edit, 'string', sprintf('%3.0f', tmp1(2)))
                if (tmp1(2) > get(cb_max_slide, 'max'))
                    set(cb_max_slide, 'value', get(cb_max_slide, 'max'))
                else
                    set(cb_max_slide, 'value', tmp1(2))
                end
            end
            set(cb_min_edit, 'string', sprintf('%3.0f', tmp1(1)))
        else
            if (tmp1(1) < get(cb_min_slide, 'min'))
                set(cb_min_slide, 'value', get(cb_min_slide, 'min'))
            else
                set(cb_min_slide, 'value', tmp1(1))
            end
        end
        switch disp_type
            case {'twtt' '~depth' 'flat'}
                [db_min, db_max] ...
                    = deal(tmp1(1), tmp1(2));
            case 'phase'
                [phase_diff_min, phase_diff_max] ...
                    = deal(tmp1(1), tmp1(2));
            case 'aresp'
                [aresp_min, aresp_max] ...
                    = deal(tmp1(1), tmp1(2));
        end
        update_db_range
        set(cb_min_slide, 'enable', 'off')
        drawnow
        set(cb_min_slide, 'enable', 'on')
    end

%% Update maximum dB/phase/ARESP

    function slide_db_max(source, eventdata)
        switch disp_type
            case {'twtt' '~depth' 'flat'}
                tmp1        = [db_min db_max];
                tmp2        = [db_min_ref db_max_ref];
            case 'phase'
                tmp1        = [phase_diff_min phase_diff_max];
                tmp2        = [phase_diff_min_ref phase_diff_max_ref];
            case 'ARESP'
                tmp1        = [aresp_min aresp_max];
                tmp2        = [aresp_min_ref aresp_max_ref];
        end
        if (get(cb_max_slide, 'value') > tmp1(1))
            if get(cbfix_check1, 'value')
                tmp3        = diff(tmp1);
            end
            tmp1(2)         = get(cb_max_slide, 'value');
            if get(cbfix_check1, 'value')
                tmp1(1)     = tmp1(2) - tmp1;
                if (tmp1(1) < tmp2(1))
                    tmp1(1) = tmp2(1);
                    tmp1(2) = tmp1(1) + tmp1;
                    if (tmp1(2) > tmp2(2))
                        tmp1(2) ...
                            = tmp2(2);
                    end
                    if (tmp1(2) > get(cb_max_slide, 'max'))
                        set(cb_max_slide, 'value', get(cb_max_slide, 'max'))
                    else
                        set(cb_max_slide, 'value', tmp1(2))
                    end
                end
                set(cb_min_edit, 'string', sprintf('%3.0f', tmp1(1)))
                if (tmp1(1) < get(cb_min_slide, 'min'))
                    set(cb_min_slide, 'value', get(cb_min_slide, 'min'))
                else
                    set(cb_min_slide, 'value', tmp1(1))
                end
            end
            set(cb_max_edit, 'string', sprintf('%3.0f', tmp1(2)))
        else
            if (tmp1(2) > get(cb_max_slide, 'max'))
                set(cb_max_slide, 'value', get(cb_max_slide, 'max'))
            else
                set(cb_max_slide, 'value', tmp1(2))
            end
        end
        switch disp_type
            case {'twtt' '~depth' 'flat'}
                [db_min, db_max] ...
                    = deal(tmp1(1), tmp1(2));
            case 'phase'
                [phase_diff_min, phase_diff_max] ...
                    = deal(tmp1(1), tmp1(2));
            case 'aresp'
                [aresp_min, aresp_max] ...
                    = deal(tmp1(1), tmp1(2));
        end
        update_db_range
        set(cb_max_slide, 'enable', 'off')
        drawnow
        set(cb_max_slide, 'enable', 'on')
    end

%% Reset minimum dB/phase/ARESP

    function reset_db_min(source, eventdata)
        switch disp_type
            case {'twtt' '~depth' 'flat'}
                [tmp1, db_min] ...
                            = deal(db_min_ref);
            case 'phase'
                [tmp1, phase_diff_min] ...
                            = deal(phase_diff_min_ref);
            case 'ARESP'
                [tmp1, aresp_min_ref] ...
                            = deal(aresp_min_ref);
        end
        if (tmp1 < get(cb_min_slide, 'min'))
            set(cb_min_slide, 'value', get(cb_min_slide, 'min'))
        else
            set(cb_min_slide, 'value', tmp1)
        end
        set(cb_min_edit, 'string', num2str(tmp1))
        update_db_range
    end

%% Reset maximum dB/phase/ARESP

    function reset_db_max(source, eventdata)
        switch disp_type
            case {'twtt' '~depth' 'flat'}
                [tmp1, db_max] ...
                            = deal(db_max_ref);
            case 'phase'
                [tmp1, phase_diff_max] ...
                            = deal(phase_diff_max_ref);
            case 'ARESP'
                [tmp1, aresp_max_ref] ...
                            = deal(aresp_max_ref);
        end
        if (tmp1 > get(cb_max_slide, 'max'))
            set(cb_max_slide, 'value', get(cb_max_slide, 'max'))
        else
            set(cb_max_slide, 'value', tmp1)
        end
        set(cb_max_edit, 'string', num2str(tmp1))
        update_db_range
    end

%% Update dB/phase/ARESP

    function update_db_range(source, eventdata)
        axes(ax_radar)
        switch disp_type
            case {'twtt' '~depth' 'flat'}
                caxis([db_min db_max])
            case 'phase'
                caxis([phase_diff_min phase_diff_max])
            case 'ARESP'
                caxis([aresp_min aresp_max])
        end
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

%% Reset both distance (x) and two-way traveltime (y)

    function reset_xy(source, eventdata)
        reset_twtt_min
        reset_twtt_max
        reset_dist_min
        reset_dist_max
    end

%% Adjust slider limits after panning or zooming

    function panzoom(source, eventdata)
        if ~load_done
            set(status_box, 'string', 'No data loaded yet.')
            return
        end
        tmp1                = get(ax_radar, 'xlim');
        if (tmp1(1) < dist_min_ref)
            reset_dist_min
        else
            if (tmp1(1) < get(dist_min_slide, 'min'))
                tmp1(1)     = dist_min_ref;
            end
            set(dist_min_slide, 'value', tmp1(1))
            set(dist_min_edit, 'string', sprintf('%3.1f', tmp1(1)))
            dist_min        = tmp1(1);
        end
        if (tmp1(2) > dist_max_ref)
            reset_dist_max
        else
            if (tmp1(2) > get(dist_max_slide, 'max'))
                tmp1(2)     = dist_max_ref;
            end
            set(dist_max_slide, 'value', tmp1(2))
            set(dist_max_edit, 'string', sprintf('%3.1f', tmp1(2)))
            dist_max        = tmp1(2);
        end
        tmp1                = get(ax_radar, 'ylim');
        if (tmp1(1) < (1e6 * twtt_min_ref))
            reset_twtt_min
        else
            if (((1e6 * twtt_max_ref) - (tmp1(1) - (1e6 * twtt_min_ref))) < get(twtt_min_slide, 'min'))
                set(twtt_min_slide, 'value', get(twtt_min_slide, 'min'))
            elseif (((1e6 * twtt_max_ref) - (tmp1(1) - (1e6 * twtt_min_ref))) > get(twtt_min_slide, 'max'))
                set(twtt_min_slide, 'value', get(twtt_min_slide, 'max'))
            else
                set(twtt_min_slide, 'value', ((1e6 * twtt_max_ref) - (tmp1(1) - (1e6 * twtt_min_ref))))
            end
            set(twtt_min_edit, 'string', sprintf('%3.1f', tmp1(1)))
            twtt_min        = 1e-6 * tmp1(1);
        end
        if (tmp1(2) > (1e6 * twtt_max_ref))
            reset_twtt_max
        else
            if (((1e6 * twtt_max_ref) - (tmp1(2) - (1e6 * twtt_min_ref))) < get(twtt_max_slide, 'min'))
                set(twtt_max_slide, 'value', get(twtt_max_slide, 'min'))
            elseif (((1e6 * twtt_max_ref) - (tmp1(2) - (1e6 * twtt_min_ref))) > get(twtt_max_slide, 'max'))
                set(twtt_max_slide, 'value', get(twtt_max_slide, 'max'))
            else
                set(twtt_max_slide, 'value', ((1e6 * twtt_max_ref) - (tmp1(2) - (1e6 * twtt_min_ref))))
            end
            set(twtt_max_edit, 'string', sprintf('%3.1f', tmp1(2)))
            twtt_max        = 1e-6 * tmp1(2);
        end
        narrow_cb
    end

%%  Arrow pan/zoom functions

    function pan_left(source, eventdata)
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
    end

    function pan_right(source, eventdata)
        tmp1                = dist_max - dist_min;
        tmp2                = dist_max + (0.25 * tmp1);
        if (tmp2 > dist_max_ref)
            dist_max        = dist_max_ref;
        else
            dist_max        = tmp2;
        end
        dist_min    = dist_max - tmp1;
        if (dist_min < dist_min_ref)
            dist_min        = dist_min_ref;
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
    end

    function pan_up(source, eventdata)
        tmp1                = twtt_max - twtt_min;
        tmp2                = twtt_min - (0.25 * tmp1);
        if (tmp2 < twtt_min_ref)
            twtt_min        = twtt_min_ref;
        else
            twtt_min        = tmp2;
        end
        twtt_max            = twtt_min + tmp1;
        set(twtt_min_edit, 'string', sprintf('%3.1f', (1e6 * twtt_min)))
        set(twtt_max_edit, 'string', sprintf('%3.1f', (1e6 * twtt_max)))
        if ((1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref))) < get(twtt_min_slide, 'min'))
            set(twtt_min_slide, 'value', get(twtt_min_slide, 'min'))
        elseif ((1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref))) > get(twtt_min_slide, 'max'))
            set(twtt_min_slide, 'value', get(twtt_min_slide, 'max'))
        else
            set(twtt_min_slide, 'value', (1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref))))
        end
        if ((1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref))) < get(twtt_max_slide, 'min'))
            set(twtt_max_slide, 'value', get(twtt_max_slide, 'min'))
        elseif ((1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref))) > get(twtt_max_slide, 'max'))
            set(twtt_max_slide, 'value', get(twtt_max_slide, 'max'))
        else
            set(twtt_max_slide, 'value', (1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref))))
        end
        update_twtt_range
    end

    function pan_down(source, eventdata)
        tmp1                = twtt_max - twtt_min;
        tmp2                = twtt_max + (0.25 * tmp1);
        if (tmp2 > twtt_max_ref)
            twtt_max        = twtt_max_ref;
        else
            twtt_max        = tmp2;
        end
        twtt_min            = twtt_max - tmp1;
        set(twtt_min_edit, 'string', sprintf('%3.1f', (1e6 * twtt_min)))
        set(twtt_max_edit, 'string', sprintf('%3.1f', (1e6 * twtt_max)))
        if ((1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref))) < get(twtt_min_slide, 'min'))
            set(twtt_min_slide, 'value', get(twtt_min_slide, 'min'))
        elseif ((1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref))) > get(twtt_min_slide, 'max'))
            set(twtt_min_slide, 'value', get(twtt_min_slide, 'max'))
        else
            set(twtt_min_slide, 'value', (1e6 * (twtt_max_ref - (twtt_min - twtt_min_ref))))
        end
        if ((1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref))) < get(twtt_max_slide, 'min'))
            set(twtt_max_slide, 'value', get(twtt_max_slide, 'min'))
        elseif ((1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref))) > get(twtt_max_slide, 'max'))
            set(twtt_max_slide, 'value', get(twtt_max_slide, 'max'))
        else
            set(twtt_max_slide, 'value', (1e6 * (twtt_max_ref - (twtt_max - twtt_min_ref))))
        end
        update_twtt_range
    end

    function zoom_in(source, eventdata)
        tmp1                = dist_max - dist_min;
        tmp2                = [(dist_min + (0.25 * tmp1)) (dist_max - (0.25 * tmp1))];
        tmp4                = twtt_max - twtt_min;
        tmp4                = [(twtt_min + (0.25 * tmp4)) (twtt_max - (0.25 * tmp4))];
        set(ax_radar, 'xlim', tmp2, 'ylim', (1e6 .* tmp4))
        panzoom
    end

    function zoom_out(source, eventdata)
        tmp1                = dist_max - dist_min;
        tmp2                = [(dist_min - (0.25 * tmp1)) (dist_max + (0.25 * tmp1))];
        tmp4                = twtt_max - twtt_min;
        tmp4                = [(twtt_min - (0.25 * tmp4)) (twtt_max + (0.25 * tmp4))];
        set(ax_radar, 'xlim', tmp2, 'ylim', (1e6 .* tmp4))
        panzoom
    end

%% Switch display type 

    function disp_radio(~, eventdata)
        disp_type           = get(eventdata.NewValue, 'string');
        switch disp_type
            case 'twtt'
                plot_twtt
            case '~depth'
                plot_depth
            case 'phase'
                plot_phase_diff
            case 'ARESP'
                plot_aresp
            case 'flat'
                plot_flat
        end
        set(disp_check, 'enable', 'off')
        drawnow
        set(disp_check, 'enable', 'on')
    end

%% Plot traveltime

    function plot_twtt(source, eventdata)
        if ~load_done
            set(status_box, 'string', 'Data not yet loaded.')
            return
        end
        if (logical(p_data) && ishandle(p_data))
            delete(p_data)
        end
        axes(ax_radar)
        if (get(cmap_list, 'value') ~= 1)
            set(cmap_list, 'value', 1)
            change_cmap
        end
        p_data              = imagesc(block.dist_lin(ind_decim), (1e6 .* block.twtt), amp_mean, [db_min db_max]);
        disp_type           = 'twtt';
        narrow_cb
        show_surfbed
        show_phase
        show_aresp
        show_man
        show_ref
        show_pk
        show_smooth
        set(cbl, 'string', '(dB)')
        set(cb_min_slide, 'min', db_min_ref, 'max', db_max_ref, 'value', db_min)
        set(cb_max_slide, 'min', db_min_ref, 'max', db_max_ref, 'value', db_max)
        set(cb_min_edit, 'string', sprintf('%3.0f', db_min))
        set(cb_max_edit, 'string', sprintf('%3.0f', db_max))
    end

%% Plot adjusted traveltime (~depth)

    function plot_depth(source, eventdata)
        if ~load_done
            set(status_box, 'string', 'Data not yet loaded.')
            return
        end
        if (logical(p_data) && ishandle(p_data))
            delete(p_data)
        end
        axes(ax_radar)
        if (get(cmap_list, 'value') ~= 1)
            set(cmap_list, 'value', 1)
            change_cmap
        end
        p_data              = imagesc(block.dist_lin(ind_decim), (1e6 .* block.twtt), amp_depth, [db_min db_max]);
        disp_type           = '~depth';
        narrow_cb
        show_surfbed
        show_phase
        show_aresp
        show_man
        show_ref
        show_pk
        show_smooth
        set(cbl, 'string', '(dB)')
        set(cb_min_slide, 'min', db_min_ref, 'max', db_max_ref, 'value', db_min)
        set(cb_max_slide, 'min', db_min_ref, 'max', db_max_ref, 'value', db_max)
        set(cb_min_edit, 'string', sprintf('%3.0f', db_min))
        set(cb_max_edit, 'string', sprintf('%3.0f', db_max))
    end

%% Plot phase difference

    function plot_phase_diff(source, eventdata)
        if ~phase_avail
            set(disp_group, 'selectedobject', disp_check(1))
            disp_type       = 'twtt';
            plot_twtt
            return
        end
        if (logical(p_data) && ishandle(p_data)) % get rid of old plotted data
            delete(p_data)
        end
        axes(ax_radar) %#ok<*MAXES>
        if (get(cmap_list, 'value') ~= 2)
            set(cmap_list, 'value', 2)
            change_cmap
        end
        p_data              = imagesc(block.dist_lin(1:decim:size(block.phase_diff_filt, 2)), (1e6 .* block.twtt), block.phase_diff_filt(:, 1:decim:end), [phase_diff_min phase_diff_max]);
        disp_type           = 'phase';
        if get(cbfix_check2, 'value')
            set(cbfix_check2, 'value', 0)
        end
        show_surfbed
        show_phase
        show_aresp
        show_man
        show_ref
        show_pk
        show_smooth
        set(cbl, 'string', '(rad)')
        set(cb_min_slide, 'min', phase_diff_min_ref, 'max', phase_diff_max_ref, 'value', phase_diff_min)
        set(cb_max_slide, 'min', phase_diff_min_ref, 'max', phase_diff_max_ref, 'value', phase_diff_max)
        set(cb_min_edit, 'string', sprintf('%2.1f', phase_diff_min))
        set(cb_max_edit, 'string', sprintf('%2.1f', phase_diff_max))
    end

%% Plot ARESP-derived layer slope

    function plot_aresp(source, eventdata)
        if ~aresp_avail
            set(disp_group, 'selectedobject', disp_check(1))
            disp_type       = 'twtt';
            plot_twtt
            return
        end
        if (logical(p_data) && ishandle(p_data)) % get rid of old plotted data
            delete(p_data)
        end
        axes(ax_radar) %#ok<*MAXES>
        if (get(cmap_list, 'value') ~= 2)
            set(cmap_list, 'value', 2)
            change_cmap
        end
        p_data              = imagesc(block.dist_lin(1:decim:size(block.slope_aresp, 2)), (1e6 .* block.twtt), atand(block.slope_aresp(:, 1:decim:end)), [aresp_min aresp_max]);
        disp_type           = 'ARESP';
        if get(cbfix_check2, 'value')
            set(cbfix_check2, 'value', 0)
        end
        show_surfbed
        show_phase
        show_aresp
        show_man
        show_ref
        show_pk
        show_smooth
        set(cbl, 'string', '(\theta)')
        set(cb_min_slide, 'min', aresp_min_ref, 'max', aresp_max_ref, 'value', aresp_min)
        set(cb_max_slide, 'min', aresp_min_ref, 'max', aresp_max_ref, 'value', aresp_max)
        set(cb_min_edit, 'string', sprintf('%2.1f', aresp_min))
        set(cb_max_edit, 'string', sprintf('%2.1f', aresp_max))
    end

%% Plot layer-flattened radargram

    function plot_flat(source, eventdata)
        if ~flat_done
            set(disp_group, 'selectedobject', disp_check(1))
            disp_type       = 'twtt';
            plot_twtt
            return
        end
        if (logical(p_data) && ishandle(p_data)) % get rid of old plotted data
            delete(p_data)
        end
        axes(ax_radar)
        if (get(cmap_list, 'value') ~= 1)
            set(cmap_list, 'value', 1)
            change_cmap
        end
        p_data              = imagesc(block.dist_lin(ind_decim_flat), (1e6 .* block.twtt), amp_flat_mean, [db_min db_max]);
        disp_type           = 'flat';
        narrow_cb
        show_surfbed
        show_phase
        show_aresp
        show_man
        show_ref
        show_pk
        show_smooth
        set(cbl, 'string', '(dB)')
        set(cb_min_slide, 'min', db_min_ref, 'max', db_max_ref, 'value', db_min)
        set(cb_max_slide, 'min', db_min_ref, 'max', db_max_ref, 'value', db_max)
        set(cb_min_edit, 'string', sprintf('%3.0f', db_min))
        set(cb_max_edit, 'string', sprintf('%3.0f', db_max))
    end

%% Show surface/bed

    function show_surfbed(source, eventdata)
        if (surf_avail || bed_avail)
            if get(surfbed_check, 'value')
                switch disp_type
                    case {'twtt' 'phase' 'ARESP'}
                        if (logical(p_bed) && ishandle(p_bed))
                            set(p_bed, 'visible', 'on')
                            uistack(p_bed, 'top')
                        end
                        if (logical(p_surf) && ishandle(p_surf))
                            set(p_surf, 'visible', 'on')
                            uistack(p_surf, 'top')
                        end
                        if (logical(p_beddepth) && ishandle(p_beddepth))
                            set(p_beddepth, 'visible', 'off')
                        end
                        if (logical(p_bedflat) && ishandle(p_bedflat))
                            set(p_bedflat, 'visible', 'off')
                        end
                        if (logical(p_surfflat) && ishandle(p_surfflat))
                            set(p_surfflat, 'visible', 'off')
                        end
                    case '~depth'
                        if (logical(p_beddepth) && ishandle(p_beddepth))
                            set(p_beddepth, 'visible', 'on')
                            uistack(p_beddepth, 'top')
                        end
                        if (logical(p_bed) && ishandle(p_bed))
                            set(p_bed, 'visible', 'off')
                        end
                        if (logical(p_surf) && ishandle(p_surf))
                            set(p_surf, 'visible', 'off')
                        end
                        if (logical(p_bedflat) && ishandle(p_bedflat))
                            set(p_bedflat, 'visible', 'off')
                        end
                        if (logical(p_surfflat) && ishandle(p_surfflat))
                            set(p_surfflat, 'visible', 'off')
                        end
                    case 'flat'
                        if (logical(p_bedflat) && ishandle(p_bedflat))
                            set(p_bedflat, 'visible', 'on')
                            uistack(p_bedflat, 'top')
                        end
                        if (logical(p_surfflat) && ishandle(p_surfflat))
                            set(p_surfflat, 'visible', 'on')
                            uistack(p_surfflat, 'top')
                        end
                        if (logical(p_bed) && ishandle(p_bed))
                            set(p_bed, 'visible', 'off')
                        end
                        if (logical(p_surf) && ishandle(p_surf))
                            set(p_surf, 'visible', 'off')
                        end
                        if (logical(p_beddepth) && ishandle(p_beddepth))
                            set(p_beddepth, 'visible', 'off')
                        end
                end
            else
                if (logical(p_bed) && ishandle(p_bed))
                    set(p_bed, 'visible', 'off')
                end
                if (logical(p_surf) && ishandle(p_surf))
                    set(p_surf, 'visible', 'off')
                end
                if (logical(p_beddepth) && ishandle(p_beddepth))
                    set(p_beddepth, 'visible', 'off')
                end
                if (logical(p_bedflat) && ishandle(p_bedflat))
                    set(p_bedflat, 'visible', 'off')
                end
                if (logical(p_surfflat) && ishandle(p_surfflat))
                    set(p_surfflat, 'visible', 'off')
                end
            end
        elseif get(surfbed_check, 'value')
            set(surfbed_check, 'value', 0)
        end
    end

%% Show phase-tracked layers

    function show_phase(source, eventdata)
        if phase_done
            if (get(phase_check, 'value') && any(strcmp(disp_type, {'twtt' '~depth'})))
                switch disp_type
                    case 'twtt'
                        if (any(p_phase) && any(ishandle(p_phase)))
                            set(p_phase(logical(p_phase) & ishandle(p_phase)), 'visible', 'on')
                            uistack(p_phase(logical(p_phase) & ishandle(p_phase)), 'top')
                        end
                        if (any(p_startphase) && any(ishandle(p_startphase)))
                            set(p_startphase(logical(p_startphase) & ishandle(p_startphase)), 'visible', 'on')
                            uistack(p_startphase(logical(p_startphase) & ishandle(p_startphase)), 'top')
                        end
                        if (any(p_startphasedepth) && any(ishandle(p_startphasedepth)))
                            set(p_startphasedepth(logical(p_startphasedepth) & ishandle(p_startphasedepth)), 'visible', 'off')
                        end
                        if (any(p_phasedepth) && any(ishandle(p_phasedepth)))
                            set(p_phasedepth(logical(p_phasedepth) & ishandle(p_phasedepth)), 'visible', 'off')
                        end
                    case '~depth'
                        if (any(p_phasedepth) && any(ishandle(p_phasedepth)))
                            set(p_phasedepth(logical(p_phasedepth) & ishandle(p_phasedepth)), 'visible', 'on')
                            uistack(p_phasedepth(logical(p_phasedepth) & ishandle(p_phasedepth)), 'top')
                        end
                        if (any(p_startphasedepth) && any(ishandle(p_startphasedepth)))
                            set(p_startphasedepth(logical(p_startphasedepth) & ishandle(p_startphasedepth)), 'visible', 'on')
                            uistack(p_startphasedepth(logical(p_startphasedepth) & ishandle(p_startphasedepth)), 'top')
                        end
                        if (any(p_startphase) && any(ishandle(p_startphase)))
                            set(p_startphase(logical(p_startphase) & ishandle(p_startphase)), 'visible', 'off')
                        end
                        if (any(p_phase) && any(ishandle(p_phase)))
                            set(p_phase(logical(p_phase) & ishandle(p_phase)), 'visible', 'off')
                        end
                end
            else
                if (any(p_startphase) && any(ishandle(p_startphase)))
                    set(p_startphase(logical(p_startphase) & ishandle(p_startphase)), 'visible', 'off')
                end
                if (any(p_phase) && any(ishandle(p_phase)))
                    set(p_phase(logical(p_phase) & ishandle(p_phase)), 'visible', 'off')
                end
                if (any(p_startphasedepth) && any(ishandle(p_startphasedepth)))
                    set(p_startphasedepth(logical(p_startphasedepth) & ishandle(p_startphasedepth)), 'visible', 'off')
                end
                if (any(p_phasedepth) && any(ishandle(p_phasedepth)))
                    set(p_phasedepth(logical(p_phasedepth) & ishandle(p_phasedepth)), 'visible', 'off')
                end
            end
        elseif get(phase_check, 'value')
            set(phase_check, 'value', 0)
        end
    end

%% Show ARESP-tracked layers

    function show_aresp(source, eventdata)
        if aresp_done
            if (get(aresp_check, 'value') && any(strcmp(disp_type, {'twtt' '~depth'})))
                switch disp_type
                    case 'twtt'
                        if (any(p_aresp) && any(ishandle(p_aresp)))
                            set(p_aresp(logical(p_aresp) & ishandle(p_aresp)), 'visible', 'on')
                            uistack(p_aresp(logical(p_aresp) & ishandle(p_aresp)), 'top')
                        end
                        if (any(p_startaresp) && any(ishandle(p_startaresp)))
                            set(p_startaresp(logical(p_startaresp) & ishandle(p_startaresp)), 'visible', 'on')
                            uistack(p_startaresp(logical(p_startaresp) & ishandle(p_startaresp)), 'top')
                        end
                        if (any(p_startarespdepth) && any(ishandle(p_startarespdepth)))
                            set(p_startarespdepth(logical(p_startarespdepth) & ishandle(p_startarespdepth)), 'visible', 'off')
                        end
                        if (any(p_arespdepth) && any(ishandle(p_arespdepth)))
                            set(p_arespdepth(logical(p_arespdepth) & ishandle(p_arespdepth)), 'visible', 'off')
                        end
                    case '~depth'
                        if (any(p_arespdepth) && any(ishandle(p_arespdepth)))
                            set(p_arespdepth(logical(p_arespdepth) & ishandle(p_arespdepth)), 'visible', 'on')
                            uistack(p_arespdepth(logical(p_arespdepth) & ishandle(p_arespdepth)), 'top')
                        end
                        if (any(p_startarespdepth) && any(ishandle(p_startarespdepth)))
                            set(p_startarespdepth(logical(p_startarespdepth) & ishandle(p_startarespdepth)), 'visible', 'on')
                            uistack(p_startarespdepth(logical(p_startarespdepth) & ishandle(p_startarespdepth)), 'top')
                        end
                        if (any(p_startaresp) && any(ishandle(p_startaresp)))
                            set(p_startaresp(logical(p_startaresp) & ishandle(p_startaresp)), 'visible', 'off')
                        end
                        if (any(p_aresp) && any(ishandle(p_aresp)))
                            set(p_aresp(logical(p_aresp) & ishandle(p_aresp)), 'visible', 'off')
                        end
                end
            else
                if (any(p_startaresp) && any(ishandle(p_startaresp)))
                    set(p_startaresp(logical(p_startaresp) & ishandle(p_startaresp)), 'visible', 'off')
                end
                if (any(p_aresp) && any(ishandle(p_aresp)))
                    set(p_aresp(logical(p_aresp) & ishandle(p_aresp)), 'visible', 'off')
                end
                if (any(p_startarespdepth) && any(ishandle(p_startarespdepth)))
                    set(p_startarespdepth(logical(p_startarespdepth) & ishandle(p_startarespdepth)), 'visible', 'off')
                end
                if (any(p_arespdepth) && any(ishandle(p_arespdepth)))
                    set(p_arespdepth(logical(p_arespdepth) & ishandle(p_arespdepth)), 'visible', 'off')
                end
            end
        elseif get(aresp_check, 'value')
            set(aresp_check, 'value', 0)
        end
    end

%% Show picked layers

    function show_pk(source, eventdata)
        if pk_done
            if (get(pk_check, 'value') && any(strcmp(disp_type, {'twtt' '~depth' 'flat'})))
                switch disp_type
                    case 'twtt'
                        if (any(p_pk) && any(ishandle(p_pk)))
                            set(p_pk(logical(p_pk) & ishandle(p_pk)), 'visible', 'on')
                            uistack(p_pk(logical(p_pk) & ishandle(p_pk)), 'top')
                        end
                        if (any(p_pkdepth) && any(ishandle(p_pkdepth)))
                            set(p_pkdepth(logical(p_pkdepth) & ishandle(p_pkdepth)), 'visible', 'off')
                        end
                        if (any(p_pkflat) && any(ishandle(p_pkflat)))
                            set(p_pkflat(logical(p_pkflat) & ishandle(p_pkflat)), 'visible', 'off')
                        end
                    case '~depth'
                        if (any(p_pkdepth) && any(ishandle(p_pkdepth)))
                            set(p_pkdepth(logical(p_pkdepth) & ishandle(p_pkdepth)), 'visible', 'on')
                            uistack(p_pkdepth(logical(p_pkdepth) & ishandle(p_pkdepth)), 'top')
                        end
                        if (any(p_pk) && any(ishandle(p_pk)))
                            set(p_pk(logical(p_pk) & ishandle(p_pk)), 'visible', 'off')
                        end
                        if (any(p_pkflat) && any(ishandle(p_pkflat)))
                            set(p_pkflat(logical(p_pkflat) & ishandle(p_pkflat)), 'visible', 'off')
                        end
                    case 'flat'
                        if (any(p_pkflat) && any(ishandle(p_pkflat)))
                            set(p_pkflat(logical(p_pkflat) & ishandle(p_pkflat)), 'visible', 'on')
                            uistack(p_pkflat(logical(p_pkflat) & ishandle(p_pkflat)), 'top')
                        end
                        if (any(p_pk) && any(ishandle(p_pk)))
                            set(p_pk(logical(p_pk) & ishandle(p_pk)), 'visible', 'off')
                        end
                        if (any(p_pkdepth) && any(ishandle(p_pkdepth)))
                            set(p_pkdepth(logical(p_pkdepth) & ishandle(p_pkdepth)), 'visible', 'off')
                        end
                end
            else
                if (any(p_pk) && any(ishandle(p_pk)))
                    set(p_pk(logical(p_pk) & ishandle(p_pk)), 'visible', 'off')
                end
                if (any(p_pkdepth) && any(ishandle(p_pkdepth)))
                    set(p_pkdepth(logical(p_pkdepth) & ishandle(p_pkdepth)), 'visible', 'off')
                end
                if (any(p_pkflat) && any(ishandle(p_pkflat)))
                    set(p_pkflat(logical(p_pkflat) & ishandle(p_pkflat)), 'visible', 'off')
                end
            end
        elseif get(pk_check, 'value')
            set(pk_check, 'value', 0)
        end
    end

%% Show smoothed layers

    function show_smooth(source, eventdata)
        if any(smooth_done)
            if (get(smooth_check, 'value') && any(strcmp(disp_type, {'twtt' '~depth' 'flat'})))
                switch disp_type
                    case 'twtt'
                        if (any(p_pksmooth) && any(ishandle(p_pksmooth)))
                            set(p_pksmooth(logical(p_pksmooth) & ishandle(p_pksmooth)), 'visible', 'on')
                            uistack(p_pksmooth(logical(p_pksmooth) & ishandle(p_pksmooth)), 'top')
                        end
                        if (any(p_pksmoothdepth) && any(ishandle(p_pksmoothdepth)))
                            set(p_pksmoothdepth(logical(p_pksmoothdepth) & ishandle(p_pksmoothdepth)), 'visible', 'off')
                        end
                        if (any(p_pksmoothflat) && any(ishandle(p_pksmoothflat)))
                            set(p_pksmoothflat(logical(p_pksmoothflat) & ishandle(p_pksmoothflat)), 'visible', 'off')
                        end
                    case '~depth'
                        if (any(p_pksmoothdepth) && any(ishandle(p_pksmoothdepth)))
                            set(p_pksmoothdepth(logical(p_pksmoothdepth) & ishandle(p_pksmoothdepth)), 'visible', 'on')
                            uistack(p_pksmoothdepth(logical(p_pksmoothdepth) & ishandle(p_pksmoothdepth)), 'top')
                        end
                        if (any(p_pksmooth) && any(ishandle(p_pksmooth)))
                            set(p_pksmooth(logical(p_pksmooth) & ishandle(p_pksmooth)), 'visible', 'off')
                        end
                        if (any(p_pksmoothflat) && any(ishandle(p_pksmoothflat)))
                            set(p_pksmoothflat(logical(p_pksmoothflat) & ishandle(p_pksmoothflat)), 'visible', 'off')
                        end
                    case 'flat'
                        if (any(p_pksmoothflat) && any(ishandle(p_pksmoothflat)))
                            set(p_pksmoothflat(logical(p_pksmoothflat) & ishandle(p_pksmoothflat)), 'visible', 'on')
                            uistack(p_pksmoothflat(logical(p_pksmoothflat) & ishandle(p_pksmoothflat)), 'top')
                        end
                        if (any(p_pksmooth) && any(ishandle(p_pksmooth)))
                            set(p_pksmooth(logical(p_pksmooth) & ishandle(p_pksmooth)), 'visible', 'off')
                        end
                        if (any(p_pksmoothdepth) && any(ishandle(p_pksmoothdepth)))
                            set(p_pksmoothdepth(logical(p_pksmoothdepth) & ishandle(p_pksmoothdepth)), 'visible', 'off')
                        end
                end
            else
                if (any(p_pksmooth) && any(ishandle(p_pksmooth)))
                    set(p_pksmooth(logical(p_pksmooth) & ishandle(p_pksmooth)), 'visible', 'off')
                end
                if (any(p_pksmoothdepth) && any(ishandle(p_pksmoothdepth)))
                    set(p_pksmoothdepth(logical(p_pksmoothdepth) & ishandle(p_pksmoothdepth)), 'visible', 'off')
                end
                if (any(p_pksmoothflat) && any(ishandle(p_pksmoothflat)))
                    set(p_pksmoothflat(logical(p_pksmoothflat) & ishandle(p_pksmoothflat)), 'visible', 'off')
                end
            end
        elseif get(smooth_check, 'value')
            set(smooth_check, 'value', 0)
        end
    end

%% Show reference layers

    function show_ref(source, eventdata)
        if ref_done
            if (get(ref_check, 'value') && any(strcmp(disp_type, {'twtt' '~depth'})))
                switch disp_type
                    case 'twtt'
                        if (any(p_ref) && any(ishandle(p_ref)))
                            set(p_ref(logical(p_ref) & ishandle(p_ref)), 'visible', 'on')
                            uistack(p_ref(logical(p_ref) & ishandle(p_ref)), 'top')
                        end
                        if (any(p_refdepth) && any(ishandle(p_refdepth)))
                            set(p_refdepth(logical(p_refdepth) & ishandle(p_refdepth)), 'visible', 'off')
                        end
                    case '~depth'
                        if (any(p_refdepth) && any(ishandle(p_refdepth)))
                            set(p_refdepth(logical(p_refdepth) & ishandle(p_refdepth)), 'visible', 'on')
                            uistack(p_refdepth(logical(p_refdepth) & ishandle(p_refdepth)), 'top')
                        end
                        if (any(p_ref) && any(ishandle(p_ref)))
                            set(p_ref(logical(p_ref) & ishandle(p_ref)), 'visible', 'off')
                        end
                end
            else
                if (any(p_ref) && any(ishandle(p_ref)))
                    set(p_ref(logical(p_ref) & ishandle(p_ref)), 'visible', 'off')
                end
                if (any(p_refdepth) && any(ishandle(p_refdepth)))
                    set(p_refdepth(logical(p_refdepth) & ishandle(p_refdepth)), 'visible', 'off')
                end
            end
        elseif get(ref_check, 'value')
            set(ref_check, 'value', 0)
        end
    end

%% Show manual layers

    function show_man(source, eventdata)
        if pk.num_man
            if (get(man_check, 'value') && any(strcmp(disp_type, {'twtt' '~depth'})))
                switch disp_type
                    case 'twtt'
                        if (any(p_man(:)) && any(ishandle(p_man(:))))
                            set(p_man(logical(p_man(:)) & ishandle(p_man(:))), 'visible', 'on')
                            uistack(p_man(logical(p_man(:)) & ishandle(p_man(:))), 'top')
                        end
                        if (any(p_mandepth(:)) && any(ishandle(p_mandepth(:))))
                            set(p_mandepth(logical(p_mandepth(:)) & ishandle(p_mandepth(:))), 'visible', 'off')
                        end
                    case '~depth'
                        if (any(p_mandepth(:)) && any(ishandle(p_mandepth(:))))
                            set(p_mandepth(logical(p_mandepth(:)) & ishandle(p_mandepth(:))), 'visible', 'on')
                            uistack(p_mandepth(logical(p_mandepth(:)) & ishandle(p_mandepth(:))), 'top')
                        end
                        if (any(p_man(:)) && any(ishandle(p_man(:))))
                            set(p_man(logical(p_man(:)) & ishandle(p_man(:))), 'visible', 'off')
                        end
                end
            else
                if (any(p_man(:)) && any(ishandle(p_man(:))))
                    set(p_man(logical(p_man(:)) & ishandle(p_man(:))), 'visible', 'off')
                end
                if (any(p_mandepth(:)) && any(ishandle(p_mandepth(:))))
                    set(p_mandepth(logical(p_mandepth(:)) & ishandle(p_mandepth(:))), 'visible', 'off')
                end
            end
        elseif get(man_check, 'value')
            set(man_check, 'value', 0)
        end
    end

%% Switch to a chunk

    function plot_chunk(source, eventdata)
        curr_chunk          = get(chunk_list, 'value');
        axes(ax_radar) %#ok<*MAXES>
        
        if (curr_chunk <= num_chunk)
            xlim(dist_chunk([curr_chunk (curr_chunk + 1)]))
        else
            xlim(block.dist_lin([1 end]))
        end
        tmp1                = get(ax_radar, 'xlim');
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
    end

%% Adjust number of indices to display

    function adj_decim(source, eventdata)
        decim               = abs(round(str2double(get(decim_edit, 'string'))));
        if load_done
            if (decim > 1)
                ind_decim   = (1 + ceil(decim / 2)):decim:(block.num_trace - ceil(decim / 2));
            else
                ind_decim   = 1:block.num_trace;
            end
            num_decim       = length(ind_decim);
            if (decim > 1)
                amp_mean    = NaN(num_sample_trim, num_decim, 'single');
                tmp1        = floor(decim / 2);
                for ii = 1:num_decim
                    amp_mean(:, ii) ...
                            = nanmean(block.amp(:, (ind_decim(ii) - tmp1):(ind_decim(ii) + tmp1)), 2);
                end
            else
                amp_mean    = single(block.amp);
            end
            if surf_avail
                if ~isempty(find(~isnan(ind_surf(ind_decim)), 1))
                    amp_depth   = NaN(num_sample_trim, num_decim, 'single');
                    for ii = find(~isnan(ind_surf(ind_decim)))
                        amp_depth(1:(num_sample_trim - ind_surf(ind_decim(ii)) + 1), ii) ...
                            = amp_mean(ind_surf(ind_decim(ii)):num_sample_trim, ii); % shift data up to surface
                    end
                    depth_avail ...
                            = true;
                    set(disp_check(2), 'visible', 'on')
                else
                    depth_avail ...
                            = false;
                    set(disp_check(2), 'visible', 'off')
                end
            end
            if (logical(p_bed) && ishandle(p_bed))
                delete(p_bed)
            end
            if (logical(p_beddepth) && ishandle(p_beddepth))
                delete(p_beddepth)
            end
            if (logical(p_surf) && ishandle(p_surf))
                delete(p_surf)
            end
            if (any(p_pkdepth) && any(ishandle(p_pkdepth)))
                delete(p_pkdepth(logical(p_pkdepth) & ishandle(p_pkdepth)))
            end            
            if (any(p_pk) && any(ishandle(p_pk)))
                delete(p_pk(logical(p_pk) & ishandle(p_pk)))
            end
            if (any(p_pkdepth) && any(ishandle(p_pkdepth)))
                delete(p_pkdepth(logical(p_pkdepth) & ishandle(p_pkdepth)))
            end
            if (any(p_pksmooth) && any(ishandle(p_pksmooth)))
                delete(p_pksmooth(logical(p_pksmooth) & ishandle(p_pksmooth)))
            end
            if (any(p_pksmoothdepth) && any(ishandle(p_pksmoothdepth)))
                delete(p_pksmoothdepth(logical(p_pksmoothdepth) & ishandle(p_pksmoothdepth)))
            end
            if surf_avail
                p_surf      = plot(block.dist_lin(ind_decim), (1e6 .* block.twtt_surf(ind_decim)), 'm.', 'markersize', 12, 'visible', 'off');
            end
            if bed_avail
                p_bed       = plot(block.dist_lin(ind_decim), (1e6 .* block.twtt_bed(ind_decim)), 'm.', 'markersize', 12, 'visible', 'off');
            end
            if (bed_avail && depth_avail)
                tmp1        = ind_decim(~isnan(ind_bed(ind_decim)) & ~isnan(ind_surf(ind_decim)));
                p_beddepth  = plot(block.dist_lin(tmp1), (1e6 .* (block.twtt_bed(tmp1) - block.twtt_surf(tmp1) + block.twtt(1))), 'm.', 'markersize', 12, 'visible', 'off');
            end
            if pk.num_layer
                [p_pk, p_pkdepth, p_pksmooth, p_pksmoothdepth] ...
                            = deal(zeros(1, pk.num_layer));
                tmp1        = [];
                for ii = 1:pk.num_layer
                    if ~isempty(find(~isnan(pk.layer(ii).ind_y(ind_decim)), 1))
                        p_pk(ii) ...
                            = plot(block.dist_lin(ind_decim(~isnan(pk.layer(ii).ind_y(ind_decim)))), (1e6 .* block.twtt(round(pk.layer(ii).ind_y(ind_decim(~isnan(pk.layer(ii).ind_y(ind_decim))))))), 'r.', 'markersize', 12, 'visible', 'off');
                        tmp2= ind_decim(~isnan(pk.layer(ii).ind_y(ind_decim)) & ~isnan(ind_surf(ind_decim)));
                        tmp3= pk.layer(ii).ind_y(tmp2) - ind_surf(tmp2) + 1;
                        tmp3((tmp3 < 1) | (tmp3 > num_sample_trim)) ...
                            = NaN;
                        p_pkdepth(ii) ...
                            = plot(block.dist_lin(tmp2(~isnan(tmp3))), (1e6 .* block.twtt(round(tmp3(~isnan(tmp3))))), 'r.', 'markersize', 12, 'visible', 'off');
                    else
                        tmp1= [tmp1 ii]; %#ok<AGROW>
                        continue
                    end
                    if (smooth_done(ii) && ~isempty(find(~isnan(pk.layer(ii).ind_y_smooth(ind_decim)), 1)))
                        p_pksmooth(ii) = plot(block.dist_lin(ind_decim(~isnan(pk.layer(ii).ind_y_smooth(ind_decim)))), (1e6 .* block.twtt(round(pk.layer(ii).ind_y_smooth(ind_decim(~isnan(pk.layer(ii).ind_y_smooth(ind_decim))))))), 'g.', 'markersize', 12, 'visible', 'off');
                        tmp2= ind_decim(~isnan(pk.layer(ii).ind_y_smooth(ind_decim)) & ~isnan(ind_surf(ind_decim)));
                        tmp3= pk.layer(ii).ind_y_smooth(tmp2) - ind_surf(tmp2) + 1;
                        tmp3((tmp3 < 1) | (tmp3 > num_sample_trim)) ...
                            = NaN;
                        p_pksmoothdepth(ii) ...
                            = plot(block.dist_lin(tmp2(~isnan(tmp3))), (1e6 .* block.twtt(round(tmp3(~isnan(tmp3))))), 'g.', 'markersize', 12, 'visible', 'off');
                    else
                        tmp1= [tmp1 ii]; %#ok<AGROW>
                    end
                end
                if ~isempty(tmp1)
                    tmp2    = setdiff(1:pk.num_layer, tmp1);
                    for ii = tmp1
                        if (logical(p_pk(ii)) && ishandle(p_pk(ii)))
                            delete(p_pk(ii))
                        end
                        if (logical(p_pkdepth(ii)) && ishandle(p_pkdepth(ii)))
                            delete(p_pkdepth(ii))
                        end
                        if (logical(p_pksmooth(ii)) && ishandle(p_pksmooth(ii)))
                            delete(p_pksmooth(ii))
                        end
                        if (logical(p_pksmoothdepth(ii)) && ishandle(p_pksmoothdepth(ii)))
                            delete(p_pksmoothdepth(ii))
                        end
                    end
                    curr_layer ...
                            = interp1(1:length(tmp2), tmp2, curr_layer, 'nearest', 'extrap');
                    [pk.num_layer, pk.layer, smooth_done, p_pk, p_pkdepth, p_pkflat, p_pksmooth, p_pksmoothdepth, p_pksmoothflat, ind_y_flat_mean, ind_y_flat_smooth] ...
                            = deal(length(tmp2), pk.layer(tmp2), smooth_done(tmp2), p_pk(tmp2), p_pkdepth(tmp2), p_pkflat(tmp2), p_pksmooth(tmp2), p_pksmoothdepth(tmp2), p_pksmoothflat(tmp2), ind_y_flat_mean(tmp2, :), ind_y_flat_smooth(tmp2, :));
                    set(layer_list, 'string', [num2cell(1:pk.num_layer) 'surface' 'bed'], 'value', curr_layer)
                end
            end
            switch disp_type
                case 'twtt'
                    plot_twtt
                case '~depth'
                    plot_depth
                case 'phase'
                    plot_phase_diff
                case 'ARESP'
                    plot_aresp
            end
        end
        set(status_box, 'string', ['Decimation set to 1/' num2str(decim) ' samples.'])
        set(decim_edit, 'string', num2str(decim))
    end

%% Adjust length of each chunk

    function adj_length_chunk(source, eventdata)
        length_chunk        = abs(round(str2double(get(length_chunk_edit, 'string'))));
        % make the horizontal chunks of data to pick (too big is unwieldy)
        num_chunk           = floor((block.dist_lin(end) - block.dist_lin(1)) ./ length_chunk);
        dist_chunk          = block.dist_lin(1):length_chunk:block.dist_lin(end);
        if (dist_chunk(end) ~= block.dist_lin(end))
            dist_chunk(end) = block.dist_lin(end);
        end
        set(chunk_list, 'string', [num2cell(1:num_chunk) 'full'], 'value', (num_chunk + 1))
        set(status_box, 'string', ['Display chunk length adjusted to ' num2str(length_chunk) ' km.'])
        set(length_chunk_edit, 'string', num2str(length_chunk))
    end

%% Adjust assumed center frequency

    function adj_freq(source, eventdata)
        pk.freq             = 1e6 * abs(str2double(get(freq_edit, 'string')));
        set(status_box, 'string', ['Radar center frequency assumed to be ' num2str(pk.freq / 1e6) ' MHz.'])
        set(freq_edit, 'string', num2str(1e-6 * pk.freq))
        if phase_done
            prop_phase
            if (any(p_phase) && any(ishandle(p_phase)))
                delete(p_phase(logical(p_phase) & ishandle(p_phase)))
            end
            if (any(p_phasedepth) && any(ishandle(p_phasedepth)))
                delete(p_phasedepth(logical(p_phasedepth) & ishandle(p_phasedepth)))
            end
            axes(ax_radar) %#ok<*LAXES>
            [p_phase, p_phasedepth] ...
                            = deal(zeros(1, pk.num_phase));
            for ii = 1:pk.num_phase
                p_phase(ii) = plot(block.dist_lin(ind_decim), (1e6 .* block.twtt(round(ind_y_phase(ii, ind_decim)))), 'b', 'linewidth', 1, 'visible', 'off');
                p_phasedepth(ii) ...
                            = plot(block.dist_lin(ind_decim), (1e6 .* block.twtt(round(ind_y_phase(ii, ind_decim) - ind_surf(ind_decim) + 1))), 'b', 'linewidth', 1, 'visible', 'off');
            end
            if keep_phase_done
                set([p_phase(pk.ind_keep_phase) p_phasedepth(pk.ind_keep_phase)], 'color', 'w', 'linewidth', 2, 'visible', 'off')
            end
            show_phase
            set(status_box, 'string', ['Adjusted phase-tracked layers to ' num2str(1e-6 * pk.freq) ' MHz'])
        end
    end

%% Adjust number of indices over which to horizontally average flattened data

    function adj_decim_flat(source, eventdata)
        decim_flat          = abs(round(str2double(get(decim_flat_edit, 'string'))));
        set(status_box, 'string', ['Flattened averaging length adjusted to ' num2str(decim_flat) ' samples.'])
        set(decim_flat_edit, 'string', num2str(decim_flat))
        if flat_done
            mean_flat
        end
    end

%% Adjust number of indices above/below layer to search for peak

    function adj_num_win(source, eventdata)
        if ~round(str2double(get(num_win_edit, 'string')))
            set(num_win_edit, 'string', num2str(pk.num_win))
            set(status_box, 'string', 'Search window must be non-zero.')
            return
        end
        pk.num_win          = abs(round(str2double(get(num_win_edit, 'string'))));
        set(num_win_edit, 'string', num2str(pk.num_win))
        set(status_box, 'string', ['Vertical search window adjusted to +/- ' num2str(pk.num_win) ' samples.'])
    end

%% Adjust layer smoothing length

    function adj_length_smooth(source, eventdata)
        if ~round(10 * str2double(get(length_smooth_edit, 'string')))
            set(length_smooth_edit, 'string', num2str(pk.length_smooth))
            set(status_box, 'string', 'Smoothing length must be non-zero.')
            return
        end
        pk.length_smooth    = abs(round(10 * str2double(get(length_smooth_edit, 'string')))) / 10;
        smooth_done         = false(1, pk.num_layer); % now no layers have this smoothing length so all must be re-smoothed
        set(status_box, 'string', ['Layer smoothing length adjusted to ' num2str(pk.length_smooth) ' km.'])
        set(length_smooth_edit, 'string', num2str(pk.length_smooth))
        pause(0.1)
        pk_smooth
    end

%% Adjust matching range

    function adj_twtt_match(source, eventdata)
        pk.twtt_match       = 1e-6 * abs(str2double(get(twtt_match_edit, 'string')));
        set(status_box, 'string', ['Matching traveltime threshold set to ' num2str(1e6 * pk.twtt_match) ' us.'])
        set(twtt_match_edit, 'string', num2str(1e6 * pk.twtt_match))
        match_done          = false;
    end

%% Change colormap

    function change_cmap(source, eventdata)
        axes(ax_radar)
        colormap(cmaps{get(cmap_list, 'value')})
    end

%% Pop out figure

    function pop_fig(source, eventdata)
        if ~load_done
            set(status_box, 'string', 'Data must be loaded prior to popping out a figure.')
            return
        end
        % make a simple figure that also gets saved
        set(0, 'DefaultFigureWindowStyle', 'default')
        figure('position', [100 100 1200 800]);
        axis ij tight
        axis([dist_min dist_max (1e6 .* [twtt_min twtt_max])])
        hold on
        colormap(cmaps{get(cmap_list, 'value')})
        switch disp_type
            case 'twtt'
                imagesc(block.dist_lin(ind_decim), (1e6 .* block.twtt), amp_mean, [db_min db_max])
                caxis([db_min db_max])
                if get(surfbed_check, 'value')
                    if surf_avail
                        tmp1= plot(block.dist_lin(ind_decim), (1e6 .* block.twtt_surf(ind_decim)), 'm.', 'markersize', 12);
                    end
                    if bed_avail
                        tmp1= plot(block.dist_lin(ind_decim), (1e6 .* block.twtt_bed(ind_decim)), 'm.', 'markersize', 12);
                        if any(isnan(ind_bed))
                            set(tmp1, 'marker', '.', 'linestyle', 'none', 'markersize', 12)
                        end
                    end
                end
                if get(ref_check, 'value')
                    switch ref_start_or_end
                        case 'start'
                            for ii = 1:pk_ref.num_layer
                                plot(block.dist_lin(1:decim:block.ind_overlap(1)), (1e6 .* pk_ref.layer(ii).twtt_smooth(pk_ref.ind_overlap(2):decim:end)), 'y.', 'markersize', 12)
                            end
                        case 'end'
                            for ii = 1:pk_ref.num_layer
                                plot(block.dist_lin(block.ind_overlap(2):decim:end), (1e6 .* pk_ref.layer(ii).twtt_smooth(1:decim:pk_ref.ind_overlap(1))), 'y.', 'markersize', 12)
                            end
                    end
                end
                if get(phase_check, 'value')
                    for ii = 1:pk.num_phase
                        tmp1= plot(block.dist_lin(ind_decim), (1e6 .* block.twtt(round(ind_y_phase(ii, ind_decim)))), 'b', 'linewidth', 1);
                        if any(ii == pk.num_keep_phase)
                            set(tmp1, 'linewidth', 2, 'color', 'w')
                        end
                    end
                    plot((ones(1, pk.num_phase) .* block.dist_lin(pk.ind_x_start_phase)), (1e6 .* block.twtt(pk.ind_y_start_phase)), 'm.', 'markersize', 12)
                end
                if get(aresp_check, 'value')
                    for ii = 1:pk.num_aresp
                        tmp1= plot(block.dist_lin(ind_decim), (1e6 .* block.twtt(round(ind_y_aresp(ii, ind_decim)))), 'c', 'linewidth', 1);
                        if any(ii == pk.num_keep_aresp)
                            set(tmp1, 'linewidth', 2, 'color', 'w')
                        end
                    end
                    plot((ones(1, pk.num_aresp) .* block.dist_lin(pk.ind_x_start_aresp)), (1e6 .* block.twtt(pk.ind_y_start_aresp)), 'm.', 'markersize', 12)
                end
                if get(man_check, 'value')
                    for ii = 1:pk.num_man
                        plot(block.dist_lin(ind_decim), (1e6 .* block.twtt(round(pk.ind_y_man(ii, ind_decim)))), 'w', 'linewidth', 1)
                    end
                end
                if get(pk_check, 'value')
                    for ii = 1:pk.num_layer
                        if ~isempty(find(~isnan(pk.layer(ii).ind_y(ind_decim)), 1))
                            tmp1 = plot(block.dist_lin(ind_decim(~isnan(pk.layer(ii).ind_y(ind_decim)))), (1e6 .* block.twtt(round(pk.layer(ii).ind_y(ind_decim(~isnan(pk.layer(ii).ind_y(ind_decim))))))), 'r.', 'markersize', 12);
                            if (ii == curr_layer)
                                set(tmp1, 'markersize', 24)
                            end
                        end
                    end
                end
                if get(smooth_check, 'value')
                    for ii = 1:pk.num_layer
                        if ~isempty(find(~isnan(pk.layer(ii).ind_y_smooth(ind_decim)), 1))
                            tmp1    = plot(block.dist_lin(ind_decim(~isnan(pk.layer(ii).ind_y_smooth(ind_decim)))), ...
                                           (1e6 .* block.twtt(round(pk.layer(ii).ind_y_smooth(ind_decim(~isnan(pk.layer(ii).ind_y_smooth(ind_decim))))))), 'g.', 'markersize', 12);
                        end
                        if (ii == curr_layer)
                            set(tmp1, 'markersize', 24)
                        end
                    end
                end
                text((dist_max + (0.015 * (dist_max - dist_min))), (1e6 * (twtt_min - (0.04 * (twtt_max - twtt_min)))), '(dB)', 'color', 'k', 'fontsize', 20)
            case '~depth'
                imagesc(block.dist_lin(ind_decim), (1e6 .* block.twtt), amp_depth, [db_min db_max])
                caxis([db_min db_max])
                if get(surfbed_check, 'value')
                    if bed_avail
                        tmp1= ind_decim(~isnan(ind_bed) & ~isnan(ind_surf));
                        tmp2= plot(block.dist_lin(tmp1), (1e6 .* (block.twtt_bed(tmp1) - block.twtt_surf(tmp1))), 'm.', 'markersize', 12);
                    end
                end
                if get(ref_check, 'value')
                    switch ref_start_or_end
                        case 'start'
                            for ii = 1:pk_ref.num_layer
                                plot(block.dist_lin(1:decim:block.ind_overlap(1)), (1e6 .* (pk_ref.layer(ii).twtt_smooth(pk_ref.ind_overlap(2):decim:end) - pk_ref.twtt_surf(pk_ref.ind_overlap(2):decim:end))), 'y.', 'markersize', 12)
                            end
                        case 'end'
                            for ii = 1:pk_ref.num_layer
                                plot(block.dist_lin(block.ind_overlap(2):decim:end), (1e6 .* (pk_ref.layer(ii).twtt_smooth(1:decim:pk_ref.ind_overlap(1)) - pk_ref.twtt_surf(1:decim:pk_ref.ind_overlap(1)))), 'y.', 'markersize', 12)
                            end
                    end
                end
                if get(phase_check, 'value')
                    for ii = 1:pk.num_phase
                        tmp1= plot(block.dist_lin(ind_decim), (1e6 .* block.twtt(round(ind_y_phase(ii, ind_decim) - ind_surf(ind_decim) + 1))), 'b', 'linewidth', 1);
                        if any(ii == pk.num_keep_phase)
                            set(tmp1, 'linewidth', 2, 'color', 'w')
                        end
                    end
                    plot((ones(1, pk.num_phase) .* block.dist_lin(pk.ind_x_start_phase)), (1e6 .* block.twtt(pk.ind_y_start_phase - ind_surf(pk.ind_x_start_phase) + 1)), 'm.', 'markersize', 12)
                end
                if get(aresp_check, 'value')
                    for ii = 1:pk.num_aresp
                        tmp1= plot(block.dist_lin(ind_decim), (1e6 .* block.twtt(round(ind_y_aresp(ii, ind_decim) - ind_surf(ind_decim) + 1))), 'c', 'linewidth', 1);
                        if any(ii == pk.num_keep_aresp)
                            set(tmp1, 'linewidth', 2, 'color', 'w')
                        end
                    end
                    plot((ones(1, pk.num_aresp) .* block.dist_lin(pk.ind_x_start_aresp)), (1e6 .* block.twtt(pk.ind_y_start_aresp - ind_surf(pk.ind_x_start_aresp) + 1)), 'm.', 'markersize', 12)
                end
                if get(man_check, 'value')
                    for ii = 1:pk.num_man
                        plot(block.dist_lin(ind_decim), (1e6 .* block.twtt(round(pk.ind_y_man(ii, ind_decim) - ind_surf(ind_decim) + 1))), 'w', 'linewidth', 1)
                    end
                end
                if get(pk_check, 'value')
                    for ii = 1:pk.num_layer
                        if ~isempty(find(~isnan(pk.layer(ii).ind_y(ind_decim)), 1))
                            tmp1 = ind_decim(~isnan(pk.layer(ii).ind_y(ind_decim)) & ~isnan(ind_surf(ind_decim)));
                            tmp2 = plot(block.dist_lin(tmp1), (1e6 .* block.twtt(round(pk.layer(ii).ind_y(tmp1) - ind_surf(tmp1) + 1))), 'r.', 'markersize', 12);
                            if (ii == curr_layer)
                                set(tmp2, 'markersize', 24)
                            end
                        end
                    end
                end
                if get(smooth_check, 'value')
                    for ii = 1:pk.num_layer
                        if ~isempty(find(~isnan(pk.layer(ii).ind_y_smooth(ind_decim)), 1))
                            tmp1    = ind_decim(~isnan(pk.layer(ii).ind_y_smooth(ind_decim)) & ~isnan(ind_surf(ind_decim)));                            
                            tmp2    = plot(block.dist_lin(tmp1), (1e6 .* block.twtt(round(pk.layer(ii).ind_y_smooth(tmp1) - ind_surf(tmp1) + 1))), 'g.', 'markersize', 12);
                        end
                        if (ii == curr_layer)
                            set(tmp2, 'markersize', 24)
                        end
                    end
                end
                text((dist_max + (0.015 * (dist_max - dist_min))), (1e6 * (twtt_min - (0.04 * (twtt_max - twtt_min)))), '(dB)', 'color', 'k', 'fontsize', 20)
            case 'phase'
                imagesc(block.dist_lin(1:decim:size(block.phase_diff_filt, 2)), (1e6 .* block.twtt), block.phase_diff_filt(:, 1:decim:end), [phase_diff_min phase_diff_max])
                if get(surfbed_check, 'value')
                    if surf_avail
                        plot(block.dist_lin(ind_decim), (1e6 .* block.twtt_surf(ind_decim)), 'm.', 'markersize', 12)
                    end
                    if bed_avail
                        plot(block.dist_lin(ind_decim), (1e6 .* block.twtt_bed(ind_decim)), 'm.', 'markersize', 12)
                    end
                end
                text((dist_max + (0.015 * (dist_max - dist_min))), (1e6 * (twtt_min - (0.04 * (twtt_max - twtt_min)))), '(rad)', 'color', 'k', 'fontsize', 20)
            case 'ARESP'
                imagesc(block.dist_lin(1:decim:size(block.slope_aresp, 2)), (1e6 .* block.twtt), block.slope_aresp(:, 1:decim:end), [aresp_min aresp_max])
                if get(surfbed_check, 'value')
                    if surf_avail
                        plot(block.dist_lin(ind_decim), (1e6 .* block.twtt_surf(ind_decim)), 'm.', 'markersize', 12)
                    end
                    if bed_avail
                        plot(block.dist_lin(ind_decim), (1e6 .* block.twtt_bed(ind_decim)), 'm.', 'markersize', 12)
                    end
                end
                text((dist_max + (0.015 * (dist_max - dist_min))), (1e6 * (twtt_min - (0.04 * (twtt_max - twtt_min)))), '(tan \theta)', 'color', 'k', 'fontsize', 20)

            case 'flat'
                imagesc(block.dist_lin(ind_decim_flat), (1e6 .* block.twtt), amp_flat_mean, [db_min db_max])
                if get(surfbed_check, 'value')
                    if surf_avail
                        tmp1        = ind_surf_flat(ind_decim_flat);
                        plot(block.dist_lin(ind_decim_flat(~isnan(tmp1))), (1e6 .* block.twtt(tmp1(~isnan(tmp1)))), 'g.', 'markersize', 12)
                    end
                    if bed_avail
                        tmp1        = ind_bed_flat(ind_decim_flat);
                        plot(block.dist_lin(ind_decim_flat(~isnan(tmp1))), (1e6 .* block.twtt(tmp1(~isnan(tmp1)))), 'g.', 'markersize', 12)
                    end
                end
                if get(pk_check, 'value')
                    for ii = 1:pk.num_layer
                        tmp1 = plot(block.dist_lin(ind_decim_flat(~isnan(ind_y_flat_mean(ii, :)))), (1e6 .* block.twtt(round(ind_y_flat_mean(ii, ~isnan(ind_y_flat_mean(ii, :)))))), 'r.', 'markersize', 12);
                        if (ii == curr_layer)
                            set(tmp1, 'markersize', 24)
                        end
                    end
                end
                if get(smooth_check, 'value')
                    for ii = 1:pk.num_layer
                        tmp1 = plot(block.dist_lin(ind_decim_flat(~isnan(ind_y_flat_smooth(ii, :)))), (1e6 .* block.twtt(round(ind_y_flat_smooth(ii, ~isnan(ind_y_flat_smooth(ii, :)))))), 'g.', 'markersize', 12);
                        if (ii == curr_layer)
                            set(tmp1, 'markersize', 24)
                        end
                    end
                end
                text((dist_max + (0.015 * (dist_max - dist_min))), (1e6 * (twtt_min - (0.04 * (twtt_max - twtt_min)))), '(dB)', 'color', 'k', 'fontsize', 20)
        end
        set(gca, 'fontsize', 20, 'layer', 'top')
        xlabel('Distance (km)')
        ylabel('Traveltime ({\mu}s)')
        title(file_data(1:(end - 4)), 'fontweight', 'bold', 'interpreter', 'none')
        colorbar('fontsize', 20)
        if get(grid_check, 'value')
            grid on
        end
        box on
        set(0, 'DefaultFigureWindowStyle', 'docked')
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
        if get(cbfix_check2, 'value')
            axes(ax_radar)
            tmp1            = zeros(2);
            tmp1(1, :)      = interp1(block.dist_lin(ind_decim), 1:num_decim, [dist_min dist_max], 'nearest', 'extrap');
            tmp1(2, :)      = interp1(block.twtt, 1:num_sample_trim, [twtt_min twtt_max], 'nearest', 'extrap');
            switch disp_type
                case 'twtt'
                    tmp1    = amp_mean(tmp1(2, 1):tmp1(2, 2), tmp1(1, 1):tmp1(1, 2));
                case '~depth'
                    tmp1    = amp_depth(tmp1(2, 1):tmp1(2, 2), tmp1(1, 1):tmp1(1, 2));
                case 'phase'
                    tmp1    = block.phase_diff_filt(tmp1(2, 1):tmp1(2, 2), tmp1(1, 1):decim:tmp1(1, 2));
                case 'ARESP'
                    tmp1    = block.slope_aresp(tmp1(2, 1):tmp1(2, 2), tmp1(1, 1):decim:tmp1(1, 2));
                case 'flat'
                    tmp1(1, :) ...
                            = interp1(block.dist_lin(ind_decim_flat), 1:num_decim_flat, [dist_min dist_max], 'nearest', 'extrap');
                    tmp1    = amp_flat_mean(tmp1(2, 1):tmp1(2, 2), tmp1(1, 1):tmp1(1, 2));
            end
            tmp2            = NaN(1, 2);
            [tmp2(1), tmp2(2)] ...
                            = deal(nanmean(tmp1(~isinf(tmp1))), nanstd(tmp1(~isinf(tmp1))));
            [tmp1, tmp4]    = deal(zeros(1, 2));
            switch disp_type
                case {'twtt' '~depth' 'flat'}
                    tmp4    = [db_min_ref db_max_ref];
                case 'phase'
                    tmp4    = [phase_diff_min_ref phase_diff_max_ref];
                case 'aresp'
                    tmp4    = [aresp_min_ref aresp_max_ref];
            end
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
            if (tmp4(1) < get(cb_min_slide, 'min'))
                tmp4(1)     = db_min_ref;
            end
            set(cb_min_slide, 'value', tmp4(1))
            if (tmp4(2) > get(cb_max_slide, 'max'))
                tmp4(2)     = db_min_ref;                
            end
            set(cb_max_slide, 'value', tmp4(2))
            set(cb_min_edit, 'string', sprintf('%3.0f', tmp4(1)))
            set(cb_max_edit, 'string', sprintf('%3.0f', tmp4(2)))
            caxis(tmp4)
            switch disp_type
                case {'twtt' '~depth' 'flat'}
                    [db_min, db_max] ...
                            = deal(tmp4(1), tmp4(2));
                case 'phase'
                    [phase_diff_min, phase_diff_max] ...
                            = deal(tmp4(1), tmp4(2));
                case 'aresp'
                    [aresp_min, aresp_max] ...
                            = deal(tmp4(1), tmp4(2));
            end
        end
    end

%% Keyboard shortcuts for various functions

    function keypress(~, eventdata)
        switch eventdata.Key
            case '1'
                if get(ref_check, 'value')
                    set(ref_check, 'value', 0)
                else
                    set(ref_check, 'value', 1)
                end
                show_ref
            case '2'
                if get(phase_check, 'value')
                    set(phase_check, 'value', 0)
                else
                    set(phase_check, 'value', 1)
                end
                show_phase
            case '3'
                if get(aresp_check, 'value')
                    set(aresp_check, 'value', 0)
                else
                    set(aresp_check, 'value', 1)
                end
                show_aresp
            case '4'
                if get(man_check, 'value')
                    set(man_check, 'value', 0)
                else
                    set(man_check, 'value', 1)
                end
                show_man
            case '5'
                if get(pk_check, 'value')
                    set(pk_check, 'value', 0)
                else
                    set(pk_check, 'value', 1)
                end
                show_pk
            case '6'
                if get(smooth_check, 'value')
                    set(smooth_check, 'value', 0)
                else
                    set(smooth_check, 'value', 1)
                end
                show_smooth
            case '7'
                if get(surfbed_check, 'value')
                    set(surfbed_check, 'value', 0)
                else
                    set(surfbed_check, 'value', 1)
                end
                show_surfbed
            case 'a'
                pk_adj
            case 'b'
                if get(cbfix_check2, 'value')
                    set(cbfix_check2, 'value', 0)
                else
                    set(cbfix_check2, 'value', 1)
                end
                narrow_cb
            case 'c'
                track_phase
            case 'd'
                pk_del
            case 'e'
                reset_xy
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
                pk_match
            case 'i'
                load_pk
            case 'j'
                track_aresp
            case 'k'
                if phase_avail
                    pk_keep_phase
                elseif aresp_avail
                    pk_keep_aresp
                end
            case 'l'
                load_data
            case 'm'
                pk_merge
            case 'n'
                pk_next
            case 'o'
                pk_focus
            case 'p'
                pk_auto
            case 'q'
                pop_fig
            case 'r'
                load_ref
            case 's'
                pk_man
            case 't'
                trim_y
            case 'u'
                track_man
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
                if get(twttfix_check, 'value')
                    set(twttfix_check, 'value', 0)
                else
                    set(twttfix_check, 'value', 1)
                end
            case 'z'
                pk_surfbed
            case 'slash'
                switch ord_poly
                    case 2
                        ord_poly ...
                            = 3;
                        set(status_box, 'string', 'Flattening polynomial switched to 3rd order')
                    case 3
                        ord_poly ...
                            = 2;
                        set(status_box, 'string', 'Flattening polynomial switched to 2nd order.')
                end
            case 'backslash'
                switch pk.predict_or_pk
                    case 'predict'
                        if pk.num_layer
                            pk.predict_or_pk ...
                            = 'pk';
                            set(status_box, 'string', 'Flattening switched to using picked layers.')
                        end
                    case 'pk'
                        pk.predict_or_pk ...
                            = 'predict';
                        set(status_box, 'string', 'Flattening switched to using predicted layers.')
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
                            set(disp_group, 'selectedobject', disp_check(5))
                            disp_type = 'flat';
                            plot_flat
                        else
                            set(disp_group, 'selectedobject', disp_check(2))
                            disp_type = '~depth';
                            plot_depth
                        end
                    case '~depth'
                        if flat_done
                            set(disp_group, 'selectedobject', disp_check(5))
                            disp_type = 'flat';
                            plot_flat
                        else
                            set(disp_group, 'selectedobject', disp_check(1))
                            disp_type = 'twtt';
                            plot_twtt
                        end
                    otherwise
                        set(disp_group, 'selectedobject', disp_check(1))
                        disp_type = 'twtt';
                        plot_twtt
                end
        end
    end

%% Mouse wheel shortcut

    function wheel_zoom(~, eventdata)
        switch eventdata.VerticalScrollCount
            case -1
                zoom_in
            case 1
                zoom_out
        end
    end

%% Mouse click

    function mouse_click(source, eventdata)
        if ((logical(pk.num_layer) || surf_avail || bed_avail) && any(strcmp(disp_type, {'twtt' '~depth' 'flat'})))
            tmp1            = get(source, 'currentpoint');
            tmp2            = get(pkgui, 'position');
            tmp3            = get(ax_radar, 'position');
            tmp4            = [(tmp2(1) + (tmp2(3) * tmp3(1))) (tmp2(1) + (tmp2(3) * (tmp3(1) + tmp3(3)))); (tmp2(2) + (tmp2(4) * tmp3(2))) (tmp2(2) + (tmp2(4) * (tmp3(2) + tmp3(4))))];
            if ((tmp1(1) > (tmp4(1, 1))) && (tmp1(1) < (tmp4(1, 2))) && (tmp1(2) > (tmp4(2, 1))) && (tmp1(2) < (tmp4(2, 2))))
                tmp1        = [((tmp1(1) - tmp4(1, 1)) / diff(tmp4(1, :))) ((tmp4(2, 2) - tmp1(2)) / diff(tmp4(2, :)))];
                tmp2        = [get(ax_radar, 'xlim'); get(ax_radar, 'ylim')];
                [ind_x_pk, ind_y_pk] ...
                            = deal(((tmp1(1) * diff(tmp2(1, :))) + tmp2(1, 1)), ((tmp1(2) * diff(tmp2(2, :))) + tmp2(2, 1)));
                pk_select_gui
            end
        end
    end

%% Test something

    function misctest(source, eventdata)
        
    end

%%
end