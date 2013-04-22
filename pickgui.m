function pickgui
% PICKGUI Interactive radar-layer picker.
%   
%   PICKGUI loads a GUI for tracing layers in CReSIS radar data, either the
%   original data files or those that have been pre-processed by
%   RADBLOCKPROC. This GUI includes semi-automatic tracing of layers and
%   layer prediction using the horizontal phase gradient and/or ARESP, if
%   the layer slopes for those methods have been added to the blocks by
%   PHASEINTERP or ARESP, respectively. Layers from separate data blocks
%   can be merged later using MERGEGUI and compared between transects using
%   FENCEGUI.
%   
%   Refer to manual for operation (pickgui_man.docx).
%   
%   PICKGUI requires that the function SMOOTH_LOWESS be available within
%   the user's path. If original CReSIS data blocks are to be loaded, then
%   it also requires that either the function LL2PS be available within the
%   user's path or the Mapping Toolbox be licensed and available. If the
%   Parallel Computing Toolbox is licensed and available, then several
%   calculations related to data flattening will be parallelized.
% 
% Joe MacGregor (UTIG), Mark Fahnestock (UAF-GI)
% Last updated: 04/17/13

if ~exist('smooth_lowess', 'file')
    error('pickgui:smoothlowess', 'Function SMOOTH_LOWESS is not available within this user''s path.')
end

%% Intialize variables

[pk, pk_ref]                = deal(struct);
pk.layer                    = struct;
[pk.num_layer, pk.num_man, pk.ind_trim_start, pk.ind_y_man, pk.num_keep_phase, pk.num_keep_aresp, pk.num_phase, pk.num_aresp, pk.ind_keep_phase, pk.ind_keep_aresp, pk.ind_x_start_phase, pk.ind_x_start_aresp, ...
 pk.ind_y_start_phase, pk.ind_y_start_aresp, pk.ind_y_phase_max, pk.ind_y_aresp_max] ...
                            = deal(0);
pk.keep_or_flat             = 'keep';

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

% specularity default
[spec_min_ref, spec_max_ref]= deal(0, 1);
[spec_min, spec_max]        = deal(0, 0.2);

% ARESP default
[aresp_min_ref, aresp_max_ref] ...
                            = deal(0, 1);
[aresp_min, aresp_max]      = deal(0, 1);

% some default values
speed_vacuum                = 299792458;
permitt_ice                 = 3.15;
speed_ice                   = speed_vacuum / sqrt(permitt_ice);
decim                       = 5; % decimate radargram for display
length_chunk                = 10; % chunk length in km
int_track                   = 10; % number of indices (vertical) to separate phase-tracked layers
pk.num_ind_mean             = 5; % number of indices over which to average
pk.num_win                  = 1; % +/- number of vertical indices in window within which to search for peak/trough in flattened layers
pk.length_smooth            = 1; % km length over which layers will be smoothed
pk.freq                     = 195e6; % radar center frequency
pk.twtt_match               = 0.15e-6; % traveltime range about which to search for matching layers
[num_ind_mean_ref, num_win_ref, length_smooth_ref, freq_ref, twtt_match_ref] ...
                            = deal(pk.num_ind_mean, pk.num_win, pk.length_smooth, pk.freq, pk.twtt_match);
disp_type                   = 'amp.';
cmaps                       = {'bone' 'jet'}';
flat_switch                 = 'full';
ref_start_end               = 'start';
time_pause                  = 0.05;
if ispc
    time_pause              = 1; % windows can't do pauses very well
end
if license('checkout', 'distrib_computing_toolbox')
    if ~matlabpool('size')
        matlabpool open
    end
    parallel_check          = true;
else
    parallel_check          = false;
end

% simply allocate a bunch of variables
[aresp_avail, aresp_done, bed_avail, flat_done, keep_phase_done, keep_aresp_done, load_done, load_flat, match_done, mean_done, phase_avail, phase_done, pk_done, ref_done, smooth_done, spec_avail, ...
 surf_avail, trim_done]     = deal(false);
[amp_flat_mean, amp_mean, block, button, curr_chunk, curr_layer, dist_chunk, ii, ind_bed, ind_bed_flat, ind_surf, ind_decim, ind_surf_flat, ind_x_mean_old, ind_x_pk, ind_y_aresp, ind_y_flat, ...
    ind_y_mat, ind_y_phase, ind_y_pk, jj, num_chunk, num_mean, num_sample_trim, p_aresp, p_bed, p_bedflat, p_data, p_flat, pk_tmp, p_phase, p_pk, p_pkflat, p_pkflatmark, p_man, p_ref, p_refflat, ...
    p_pksmooth, p_pksmoothflat, p_startphase, p_startaresp, p_surf, p_surfflat, pkfig, rad_sample, tmp1, tmp2, tmp3, tmp4, tmp5] ...
                            = deal(0);
amp_flat                    = NaN;
[curr_type, file_data, file_pk, file_ref, path_data, path_pk] ...
                            = deal('');

%% Draw the GUI

set(0, 'DefaultFigureWindowStyle', 'docked')
if ispc % windows switch
    pkgui                   = figure('toolbar', 'figure', 'name', 'PICKGUI', 'position', [1920 940 1 1], 'menubar', 'none', 'keypressfcn', @keypress);
    ax_radar                = subplot('position', [0.065 0.06 1.42 0.81]);
    size_font               = 14;
    width_slide             = 0.01;
else
    pkgui                   = figure('toolbar', 'figure', 'name', 'PICKGUI', 'position', [1864 1100 1 1], 'menubar', 'none', 'keypressfcn', @keypress);
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
twtt_min_slide              = uicontrol(pkgui, 'style', 'slider', 'units', 'normalized', 'position', [0.005 0.50 width_slide 0.32], 'callback', @slide_twtt_min, 'min', 0, 'max', 1, 'value', twtt_max_ref, ...
                                               'sliderstep', [0.01 0.1]);
twtt_max_slide              = uicontrol(pkgui, 'style', 'slider', 'units', 'normalized', 'position', [0.005 0.07 width_slide 0.32], 'callback', @slide_twtt_max, 'min', 0, 'max', 1, 'value', twtt_min_ref, ...
                                               'sliderstep', [0.01 0.1]);
cb_min_slide                = uicontrol(pkgui, 'style', 'slider', 'units', 'normalized', 'position', [0.96 0.07 width_slide 0.32], 'callback', @slide_db_min, 'min', -150, 'max', 0, 'value', db_min_ref, ...
                                               'sliderstep', [0.01 0.1]);
cb_max_slide                = uicontrol(pkgui, 'style', 'slider', 'units', 'normalized', 'position', [0.96 0.50 width_slide 0.32], 'callback', @slide_db_max, 'min', -150, 'max', 0, 'value', db_max_ref, ...
                                               'sliderstep', [0.01 0.1]);
dist_min_slide              = uicontrol(pkgui, 'style', 'slider', 'units', 'normalized', 'position', [0.12 0.005 0.27 0.02], 'callback', @slide_dist_min, 'min', 0, 'max', 1, 'value', dist_min_ref, ...
                                               'sliderstep', [0.01 0.1]);
dist_max_slide              = uicontrol(pkgui, 'style', 'slider', 'units', 'normalized', 'position', [0.64 0.005 0.27 0.02], 'callback', @slide_dist_max, 'min', 0, 'max', 1, 'value', dist_max_ref, ...
                                               'sliderstep', [0.01 0.1]);
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
phase_push                  = uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Phase', 'units', 'normalized', 'position', [0.37 0.965 0.04 0.03], 'callback', @track_phase, 'fontsize', size_font, ...
                                               'foregroundcolor', 'm');
keep_push                   = uicontrol(pkgui, 'style', 'pushbutton', 'string', 'keep', 'units', 'normalized', 'position', [0.42 0.965 0.04 0.03], 'callback', @pk_keep_phase, 'fontsize', size_font, ...
                                               'foregroundcolor', 'm');
aresp_push                  = uicontrol(pkgui, 'style', 'pushbutton', 'string', 'ARESP', 'units', 'normalized', 'position', [0.37 0.925 0.04 0.03], 'callback', @track_aresp, 'fontsize', size_font, ...
                                               'foregroundcolor', 'm');
keep_aresp_push             = uicontrol(pkgui, 'style', 'pushbutton', 'string', 'keep', 'units', 'normalized', 'position', [0.42 0.925 0.04 0.03], 'callback', @pk_keep_aresp, 'fontsize', size_font, ...
                                               'foregroundcolor', 'm');
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Manual', 'units', 'normalized', 'position', [0.37 0.885 0.04 0.03], 'callback', @track_man, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Flatten', 'units', 'normalized', 'position', [0.42 0.885 0.04 0.03], 'callback', @flatten, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Mean flat data', 'units', 'normalized', 'position', [0.47 0.965 0.08 0.03], 'callback', @mean_flat, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'man.', 'units', 'normalized', 'position', [0.495 0.925 0.025 0.03], 'callback', @pk_man, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'flat', 'units', 'normalized', 'position', [0.5225 0.925 0.025 0.03], 'callback', @pk_flat, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Smooth layers', 'units', 'normalized', 'position', [0.47 0.885 0.08 0.03], 'callback', @pk_smooth, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Delete', 'units', 'normalized', 'position', [0.72 0.885 0.05 0.03], 'callback', @del_layer, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Adjust', 'units', 'normalized', 'position', [0.77 0.925 0.05 0.03], 'callback', @adj_layer, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Merge', 'units', 'normalized', 'position', [0.82 0.925 0.04 0.03], 'callback', @merge_layer, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Delete all', 'units', 'normalized', 'position', [0.77 0.885 0.05 0.03], 'callback', @del_all, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Focus', 'units', 'normalized', 'position', [0.72 0.925 0.05 0.03], 'callback', @focus_layer, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'choose', 'units', 'normalized', 'position', [0.68 0.885 0.035 0.03], 'callback', @choose_pk, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'next', 'units', 'normalized', 'position', [0.68 0.925 0.035 0.03], 'callback', @pk_next, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'last', 'units', 'normalized', 'position', [0.64 0.925 0.035 0.03], 'callback', @pk_last, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'test', 'units', 'normalized', 'position', [0.82 0.885 0.04 0.03], 'callback', @misctest, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Match', 'units', 'normalized', 'position', [0.925 0.925 0.04 0.03], 'callback', @pk_match, 'fontsize', size_font, 'foregroundcolor', 'm')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'Save', 'units', 'normalized', 'position', [0.965 0.925 0.03 0.03], 'callback', @pk_save, 'fontsize', size_font, 'foregroundcolor', 'g')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'reset x/y', 'units', 'normalized', 'position', [0.945 0.885 0.05 0.03], 'callback', @reset_xy, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.005 0.46 0.03 0.03], 'callback', @reset_twtt_min, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.005 0.03 0.03 0.03], 'callback', @reset_twtt_max, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.40 0.005 0.03 0.03], 'callback', @reset_dist_min, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.92 0.005 0.03 0.03], 'callback', @reset_dist_max, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.965 0.03 0.03 0.03], 'callback', @reset_db_min, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(pkgui, 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.965 0.46 0.03 0.03], 'callback', @reset_db_max, 'fontsize', size_font, 'foregroundcolor', 'r')

% fixed text annotations
a(1)                        = annotation('textbox', [0.13 0.965 0.03 0.03], 'string', 'N_{decimate}', 'fontsize', size_font, 'color', 'b', 'edgecolor', 'none');
a(2)                        = annotation('textbox', [0.21 0.965 0.03 0.03], 'string', 'Chunk', 'fontsize', size_font, 'color', 'b', 'edgecolor', 'none');
a(3)                        = annotation('textbox', [0.21 0.925 0.03 0.03], 'string', 'L_{chunk}', 'fontsize', size_font, 'color', 'b', 'edgecolor', 'none');
a(4)                        = annotation('textbox', [0.21 0.885 0.03 0.03], 'string', 'f_{center}', 'fontsize', size_font, 'color', 'b', 'edgecolor', 'none');
a(5)                        = annotation('textbox', [0.965 0.42 0.03 0.03], 'string', 'dB_{min}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(6)                        = annotation('textbox', [0.965 0.85 0.03 0.03], 'string', 'dB_{max}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(7)                        = annotation('textbox', [0.56 0.965 0.03 0.03], 'string', 'N_{mean}', 'fontsize', size_font, 'color', 'm', 'edgecolor', 'none');
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
a(20)                       = annotation('textbox', [0.47 0.925 0.03 0.03], 'string', 'Pick', 'fontsize', size_font, 'color', 'm', 'edgecolor', 'none');
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
decim_edit                  = uicontrol(pkgui, 'style', 'edit', 'string', num2str(decim), 'units', 'normalized', 'position', [0.175 0.965 0.03 0.03], 'fontsize', size_font, 'foregroundcolor', 'k', ...
                                               'backgroundcolor', 'w', 'callback', @adj_decim);
length_chunk_edit           = uicontrol(pkgui, 'style', 'edit', 'string', num2str(length_chunk), 'units', 'normalized', 'position', [0.25 0.925 0.04 0.03], 'fontsize', size_font, 'foregroundcolor', 'k', ...
                                               'backgroundcolor', 'w', 'callback', @adj_length_chunk);
freq_edit                   = uicontrol(pkgui, 'style', 'edit', 'string', num2str(1e-6 * pk.freq), 'units', 'normalized', 'position', [0.25 0.885 0.04 0.03], 'fontsize', size_font, 'foregroundcolor', 'k', ...
                                               'backgroundcolor', 'w', 'callback', @adj_freq);
num_ind_mean_edit           = uicontrol(pkgui, 'style', 'edit', 'string', num2str(pk.num_ind_mean), 'units', 'normalized', 'position', [0.595 0.965 0.04 0.03], 'fontsize', size_font, 'foregroundcolor', 'k', ...
                                               'backgroundcolor', 'w', 'callback', @adj_num_ind_mean);
num_win_edit                = uicontrol(pkgui, 'style', 'edit', 'string', num2str(pk.num_win), 'units', 'normalized', 'position', [0.595 0.925 0.04 0.03], 'fontsize', size_font, 'foregroundcolor', 'k', ...
                                               'backgroundcolor', 'w', 'callback', @adj_num_win);
length_smooth_edit          = uicontrol(pkgui, 'style', 'edit', 'string', num2str(pk.length_smooth), 'units', 'normalized', 'position', [0.595 0.885 0.04 0.03], 'fontsize', size_font, 'foregroundcolor', 'k', ...
                                               'backgroundcolor', 'w', 'callback', @adj_length_smooth);
twtt_match_edit             = uicontrol(pkgui, 'style', 'edit', 'string', num2str(1e6 * pk.twtt_match), 'units', 'normalized', 'position', [0.88 0.925 0.04 0.03], 'fontsize', size_font, 'foregroundcolor', 'k', ...
                                               'backgroundcolor', 'w', 'callback', @adj_twtt_match);
% menus
chunk_list                  = uicontrol(pkgui, 'style', 'popupmenu', 'string', 'N/A', 'value', 1, 'units', 'normalized', 'position', [0.25 0.955 0.045 0.04], 'fontsize', size_font, 'foregroundcolor', 'k', ...
                                               'callback', @plot_chunk);
cmap_list                   = uicontrol(pkgui, 'style', 'popupmenu', 'string', cmaps, 'value', 1, 'units', 'normalized', 'position', [0.89 0.865 0.05 0.05], 'callback', @change_cmap, 'fontsize', size_font);
layer_list                  = uicontrol(pkgui, 'style', 'popupmenu', 'string', 'N/A', 'value', 1, 'units', 'normalized', 'position', [0.64 0.875 0.0375 0.04], 'fontsize', size_font, 'foregroundcolor', 'k', ...
                                               'callback', @choose_layer);

% check boxes
twttfix_check               = uicontrol(pkgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.04 0.85 0.01 0.03], 'fontsize', size_font, 'value', 0);
distfix_check               = uicontrol(pkgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.02 0.005 0.01 0.03], 'fontsize', size_font, 'value', 0);
cbfix_check1                = uicontrol(pkgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.97 0.005 0.01 0.03], 'fontsize', size_font, 'value', 0);
cbfix_check2                = uicontrol(pkgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.985 0.005 0.01 0.03], 'callback', @narrow_cb, 'fontsize', size_font, 'value', 0);
ref_check                   = uicontrol(pkgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.36 0.92 0.01 0.03], 'callback', @show_ref, 'fontsize', size_font, 'value', 0);
man_check                   = uicontrol(pkgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.41 0.88 0.01 0.03], 'callback', @show_man, 'fontsize', size_font, 'value', 0);
phase_check                 = uicontrol(pkgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.41 0.96 0.01 0.03], 'callback', @show_phase, 'fontsize', size_font, 'value', 0);
aresp_check                 = uicontrol(pkgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.41 0.92 0.01 0.03], 'callback', @show_aresp, 'fontsize', size_font, 'value', 0);
flat_check                  = uicontrol(pkgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.46 0.88 0.01 0.03], 'callback', @show_flat, 'fontsize', size_font, 'value', 0);
mean_check                  = uicontrol(pkgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.55 0.96 0.01 0.03], 'callback', @show_mean, 'fontsize', size_font, 'value', 0);
pk_check                    = uicontrol(pkgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.55 0.92 0.01 0.03], 'callback', @show_pk, 'fontsize', size_font, 'value', 0);
smooth_check                = uicontrol(pkgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.55 0.88 0.01 0.03], 'callback', @show_smooth, 'fontsize', size_font, 'value', 0);
grid_check                  = uicontrol(pkgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.88 0.885 0.01 0.03], 'callback', @toggle_grid, 'fontsize', size_font, 'value', 0);

% display buttons
disp_group                  = uibuttongroup('position', [0.005 0.885 0.20 0.03], 'selectionchangefcn', @disp_radio);
uicontrol(pkgui, 'style', 'text', 'parent', disp_group, 'units', 'normalized', 'position', [0 0.6 0.9 0.3], 'fontsize', size_font)
disp_check(1)               = uicontrol(pkgui, 'style', 'radio', 'string', 'amp.', 'units', 'normalized', 'position', [0.01 0.1 0.3 0.8], 'parent', disp_group, 'fontsize', size_font, 'handlevisibility', 'off');
disp_check(2)               = uicontrol(pkgui, 'style', 'radio', 'string', 'phase', 'units', 'normalized', 'position', [0.2 0.1 0.5 0.8], 'parent', disp_group, 'fontsize', size_font, 'handlevisibility', 'off');
disp_check(3)               = uicontrol(pkgui, 'style', 'radio', 'string', 'ARESP', 'units', 'normalized', 'position', [0.41 0.1 0.5 0.8], 'parent', disp_group, 'fontsize', size_font, 'handlevisibility', 'off');
disp_check(4)               = uicontrol(pkgui, 'style', 'radio', 'string', 'spec.', 'units', 'normalized', 'position', [0.64 0.1 0.5 0.8], 'parent', disp_group, 'fontsize', size_font, 'handlevisibility', 'off');
disp_check(5)               = uicontrol(pkgui, 'style', 'radio', 'string', 'flat', 'units', 'normalized', 'position', [0.85 0.1 0.15 0.8], 'parent', disp_group, 'fontsize', size_font, 'handlevisibility', 'off');
set(disp_group, 'selectedobject', disp_check(1))

%% Clear plots

    function clear_plots(source, eventdata)
        if (any(p_aresp) & any(ishandle(p_aresp))) %#ok<AND2>
            delete(p_aresp(logical(p_aresp) & ishandle(p_aresp)))
        end
        if (logical(p_bedflat) && ishandle(p_bedflat))
            delete(p_bedflat)
        end
        if (any(p_flat) && any(ishandle(p_flat)))
            delete(p_flat(logical(p_flat) & ishandle(p_flat)))
        end
        if (any(p_phase) & any(ishandle(p_phase))) %#ok<AND2>
            delete(p_phase(logical(p_phase) & ishandle(p_phase)))
        end
        if (any(p_man(:)) && any(ishandle(p_man(:))))
            delete(p_man(logical(p_man(:)) & ishandle(p_man(:))))
        end
        if (any(p_pk) && any(ishandle(p_pk)))
            delete(p_pk(logical(p_pk) & ishandle(p_pk)))
        end
        if (any(p_pkflat) && any(ishandle(p_pkflat)))
            delete(p_pkflat(logical(p_pkflat) & ishandle(p_pkflat)))
        end
        if (any(p_pkflatmark) && any(ishandle(p_pkflatmark)))
            delete(p_pkflatmark(logical(p_pkflatmark) & ishandle(p_pkflatmark)))
        end
        if (any(p_pksmooth) && any(ishandle(p_pksmooth)))
            delete(p_pksmooth(logical(p_pksmooth) & ishandle(p_pksmooth)))
        end
        if (any(p_pksmoothflat) && any(ishandle(p_pksmoothflat)))
            delete(p_pksmoothflat(logical(p_pksmoothflat) & ishandle(p_pksmoothflat)))
        end
        if (logical(p_startphase) && ishandle(p_startphase))
            delete(p_startphase)
        end
        if (logical(p_startaresp) && ishandle(p_startaresp))
            delete(p_startaresp)
        end
        if (logical(p_surfflat) && ishandle(p_surfflat))
            delete(p_surfflat)
        end
        pause(time_pause)
        set([aresp_check phase_check flat_check man_check mean_check pk_check smooth_check], 'value', 0)
        set(disp_group, 'selectedobject', disp_check(1))
        set(layer_list, 'string', 'N/A', 'value', 1)
    end

%% Clear data and picks

    function clear_data(source, eventdata)
        [pk, pk_ref]        = deal(struct);
        pk.layer            = struct;
        [pk.num_layer, pk.ind_trim_start, pk.num_man, pk.ind_y_man, pk.num_keep_phase, pk.num_keep_aresp, pk.num_aresp, pk.num_phase, pk.ind_keep_phase, pk.ind_keep_aresp, ...
            pk.ind_x_start_phase, pk.ind_x_start_aresp, pk.ind_y_start_phase, pk.ind_y_start_aresp, pk.ind_y_phase_max, pk.ind_y_aresp_max] ...
                            = deal(0);
        pk.keep_or_flat     = 'keep';
        flat_switch         = 'full';
        [aresp_avail, aresp_done, bed_avail, pk_done, flat_done, keep_phase_done, keep_aresp_done, load_done, load_flat, match_done, mean_done, phase_avail, phase_done, ref_done, smooth_done, spec_avail, ...
            surf_avail, trim_done] ...
                            = deal(false);
        [amp_flat_mean, amp_mean, button, curr_chunk, curr_layer, ind_decim, dist_chunk, ii, ind_bed, ind_bed_flat, ind_surf, ind_surf_flat, ind_x_mean_old, ind_x_pk, ind_y_aresp, ind_y_flat, ...
            ind_y_mat, ind_y_phase, ind_y_pk, jj, num_chunk, num_mean, num_sample_trim, p_aresp, p_bed, p_bedflat, p_data, p_flat, pk_tmp, p_phase, p_man, p_pk, p_pkflat, p_pkflatmark, p_pksmooth, ...
            p_pksmoothflat, p_ref, p_refflat, p_startphase, p_startaresp, p_surf, p_surfflat, pkfig, rad_sample, tmp1, tmp2, tmp3, tmp4, tmp5] ...
                            = deal(0);
        amp_flat            = NaN;
        [curr_type, file_ref, file_pk, ref_start_end] ...
                            = deal('');
        [pk.num_ind_mean, pk.num_win, pk.length_smooth, pk.freq, pk.twtt_match] ...
                            = deal(num_ind_mean_ref, num_win_ref, length_smooth_ref, freq_ref, twtt_match_ref);
        set(num_ind_mean_edit, 'string', num2str(pk.num_ind_mean))
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
        
        if ~isempty(file_data)
            
            % attempt to load data
            set(status_box, 'string', 'Loading data...')
            pause(time_pause)
            tmp1            = load([path_data file_data]);
            try
                block       = tmp1.block;
                tmp1        = 0;
                set(file_box, 'string', file_data(1:(end - 4)))
            catch % if not block structure than possibly an original CReSIS data file
                if isfield(tmp1, 'Data')
                    if (~license('checkout', 'map_toolbox') && ~exist('ll2ps', 'file'))
                        error('pickgui:ll2ps', 'Function LL2PS is not available within this user''s path and Mapping Toolbox is not available.')
                    elseif license('checkout', 'map_toolbox')
                        wgs84               = almanac('earth', 'wgs84', 'meters');
                        ps_struct           = defaultm('ups');
                        [ps_struct.geoid, ps_struct.mapparallels, ps_struct.falsenorthing, ps_struct.falseeasting, ps_struct.origin] ...
                                            = deal(wgs84, 70, 0, 0, [90 -45 0]);
                    end
                    block                   = struct;
                    try
                        [block.amp, block.lat, block.lon, block.num_trace, block.twtt, block.twtt_surf, block.twtt_bed, block.elev_air, block.dt, block.num_sample, block.file_in, block.ind_overlap] ...
                                            = deal(double(tmp1.Data), tmp1.Latitude, tmp1.Longitude, length(tmp1.Latitude), tmp1.Time, tmp1.Surface, tmp1.Bottom, tmp1.Elevation, (tmp1.Time(2) - tmp1.Time(1)), ...
                                                   length(tmp1.Time), file_data(1:(end - 4)), NaN(1, 2));
                    catch
                        set(status_box, 'string', 'Selected file does not contain expected variables.')
                        return
                    end
                    block.amp               = single(block.amp);
                    if isrow(block.twtt)
                        block.twtt          = block.twtt';
                    end
                    if (license('test', 'map_toolbox') && (block.lat(1) > 0)) % in Greenland and can do better x/y coordinates
                        [block.x, block.y]  = mfwdtran(ps_struct, block.lat, block.lon);
                    else
                        [block.x, block.y]  = ll2ps(block.lat, block.lon, 70); % convert lat/lon to polar stereographic x/y (better for Antarctica then Greenland)
                    end
                    [block.x, block.y]      = deal((1e-3 .* block.x), (1e-3 .* block.y)); % m to km
                    block.dist              = cumsum([0 sqrt((diff(block.x) .^ 2) + (diff(block.y) .^ 2))]); % distance vector
                    block.dist_lin          = interp1([1 block.num_trace], block.dist([1 end]), 1:block.num_trace); % monotonic distance vector
                    tmp1                    = 0;
                    set(file_box, 'string', file_data(6:(end - 4)))
                else
                    set(status_box, 'string', 'Selected file does not contain a block structure.')
                    return
                end
            end
            
            if load_done
                % delete old handles
                if (logical(p_data) && ishandle(p_data))
                    delete(p_data)
                end
                if (logical(p_surf) && ishandle(p_surf))
                    delete(p_surf)
                end
                if (logical(p_bed) && ishandle(p_bed))
                    delete(p_bed)
                end
                if (any(p_ref) && any(ishandle(p_ref)))
                    delete(p_ref(logical(p_ref) & ishandle(p_ref)))
                end
                if (any(p_refflat) && any(ishandle(p_refflat)))
                    delete(p_refflat(logical(p_refflat) & ishandle(p_refflat)))
                end
                set(ref_check, 'value', 0)
                clear_plots
                clear_data
                set(twttfix_check, 'value', 0)
                load_done   = false;
                set(disp_group, 'selectedobject', disp_check(1))
                disp_type   = 'amp.';
                pause(time_pause)
            end
            
            if ~isfield(block, 'phase_diff_filt') % no phase to work with
                set(phase_push, 'visible', 'off')
                set(phase_check, 'visible', 'off')
                set(keep_push, 'visible', 'off')
                set(disp_check(2), 'visible', 'off')
                [decim, pk.num_ind_mean] ...
                            = deal(1, 5);
                set(decim_edit, 'string', num2str(decim))
                set(num_ind_mean_edit, 'string', num2str(pk.num_ind_mean))
                pk.ind_x_mean ...
                            = ceil((pk.num_ind_mean / 2) + 1):pk.num_ind_mean:(block.num_trace - ceil(pk.num_ind_mean / 2));
                num_mean    = length(pk.ind_x_mean);
                phase_avail = false;
            else
                set(phase_push, 'visible', 'on')
                set(phase_check, 'visible', 'on')
                set(keep_push, 'visible', 'on')
                set(disp_check(2), 'visible', 'on')
                phase_avail = true;
            end
            
            if ~isfield(block, 'slope_aresp') % no aresp
                set(aresp_push, 'visible', 'off')
                set(aresp_check, 'visible', 'off')
                set(keep_aresp_push, 'visible', 'off')
                set(disp_check(3), 'visible', 'off')
                aresp_avail = false;
            else
                set(aresp_push, 'visible', 'on')
                set(aresp_check, 'visible', 'on')
                set(keep_aresp_push, 'visible', 'on')
                set(disp_check(3), 'visible', 'on')
                aresp_avail = true;
            end
            
            if ~isfield(block, 'specularity') % no specularity
                set(disp_check(4), 'visible', 'off')
                spec_avail  = false;
            else
                [spec_min_ref, spec_max_ref] ...
                            = deal(min(block.specularity(:)), max(block.specularity(:)));
                set(disp_check(4), 'visible', 'on')
                spec_avail  = true;
            end
            
            % decimation vector
            if (decim > 1)
                ind_decim   = (1 + ceil(decim / 2)):decim:(block.num_trace - ceil(decim / 2));
            else
                ind_decim   = 1:block.num_trace;
            end
            
            % convert to dB
            block.amp(block.amp == 0) ...
                            = NaN;
            block.amp(isinf(block.amp)) ...
                            = NaN;
            block.amp       = 10 .* log10(abs(block.amp));
            num_sample_trim = block.num_sample;
            if (length(block.twtt) ~= block.num_sample) % fix for data that were preprocessed with blockproc before this same fix was implemented there
                block.twtt  = interp1(1:length(block.twtt), block.twtt, 1:block.num_sample, 'linear', 'extrap');
            end
            if isrow(block.twtt)
                block.twtt  = block.twtt';
            end
            if (decim > 1)
                amp_mean    = NaN(num_sample_trim, length(ind_decim), 'single');
                tmp1        = floor(decim / 2);
                for ii = 1:length(ind_decim)
                    amp_mean(:, ii) ...
                            = nanmean(block.amp(:, (ind_decim(ii) - tmp1):(ind_decim(ii) + tmp1)), 2);
                end
            else
                amp_mean    = single(block.amp);
            end
            
            % make chunks
            adj_length_chunk
            
            % assign traveltime and distance reference values/sliders based on data
            [twtt_min_ref, twtt_max_ref, dist_min_ref, dist_max_ref, twtt_min, twtt_max, dist_min, dist_max, db_min_ref, db_max_ref, db_min, db_max] ...
                            = deal(block.twtt(1), block.twtt(end), block.dist_lin(1), block.dist_lin(end), block.twtt(1), block.twtt(end), block.dist_lin(1), block.dist_lin(end), min(block.amp(:)), ...
                                   max(block.amp(:)), min(block.amp(:)), max(block.amp(:)));
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
            
            if isfield(block, 'twtt_surf')
                ind_surf    = NaN(1, block.num_trace);
                if any(~isnan(block.twtt_surf) & any(~isinf(block.twtt_surf)))
                    ind_surf(~isnan(block.twtt_surf) & ~isinf(block.twtt_surf)) ...
                            = interp1(block.twtt, 1:num_sample_trim, block.twtt_surf(~isnan(block.twtt_surf) & ~isinf(block.twtt_surf)), 'nearest', 'extrap');
                    surf_avail ...
                            = true;
                end
            end
            if isfield(block, 'twtt_bed')
                ind_bed     = NaN(1, block.num_trace);
                if any(~isnan(block.twtt_bed) & any(~isinf(block.twtt_bed)))
                    ind_bed(~isnan(block.twtt_bed) & ~isinf(block.twtt_bed)) ...
                            = interp1(block.twtt, 1:num_sample_trim, block.twtt_bed(~isnan(block.twtt_bed) & ~isinf(block.twtt_bed)), 'nearest', 'extrap');
                    bed_avail ...
                            = true;
                end
            end
            
            % plot data
            set(disp_group, 'selectedobject', disp_check(1))
            load_done       = true;           
            plot_db
            set(status_box, 'string', 'Block loaded.')
            
        else
            file_data       = tmp1;
            set(status_box, 'string', 'No data loaded.')
        end
        
    end

%% Trim extraneous data using current y axis (e.g., before surface reflection and after deepest bed reflection)

    function trim_y(source, eventdata)
        
        if ~load_done
            set(status_box, 'string', 'No data to trim.')
            return
        end
        if (trim_done && flat_done)
            set(status_box, 'string', 'Data that has already been flattened should not be trimmed twice.')
            return
        end
        
        set(status_box, 'string', 'Trimming data...')
        pause(time_pause)
        
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
        
        if phase_avail
            block.phase_diff_filt ...
                            = block.phase_diff_filt(tmp1, :);
        end
        if aresp_avail
            block.slope_aresp ...
                            = block.slope_aresp(tmp1, :);
        end
        if spec_avail
            block.specularity ...
                            = block.specularity(tmp1, :);
        end
        num_sample_trim     = length(tmp1);
        
        % if layers already exist, then trim them too
        if pk_done
            for ii = 1:pk.num_layer
                pk.layer(ii).ind_y ...
                            = pk.layer(ii).ind_y - pk.ind_trim_start;
            end
            if (load_flat && isfield(pk.layer, 'ind_y_flat_mean'))
                for ii = 1:pk.num_layer
                    pk.layer(ii).ind_y_flat_mean ...
                            = pk.layer(ii).ind_y_flat_mean - pk.ind_trim_start;
                end
            end
            if any(smooth_done)
                for ii = find(smooth_done)
                    pk.layer(ii).ind_y_smooth ...
                            = pk.layer(ii).ind_y_smooth - pk.ind_trim_start;
                end
                if (load_flat && isfield(pk.layer, 'ind_y_flat_smooth'))
                    for ii = find(smooth_done)
                        pk.layer(ii).ind_y_flat_smooth ...
                            = pk.layer(ii).ind_y_flat_smooth - pk.ind_trim_start;
                    end
                end
            end
        end
        
        set(status_box, 'string', ['Data trimmed off before ' num2str((1e6 * block.twtt(tmp1(1))), '%2.1f') ' us and after ' num2str((1e6 * block.twtt(tmp1(end))), '%2.1f') ' us.'])
        
        % save traveltime for the end
        block.twtt          = block.twtt(tmp1);
        [twtt_min_ref, twtt_max_ref] ...
                            = deal(block.twtt(1), block.twtt(end));
        set(twtt_min_slide, 'min', (1e6 * twtt_min_ref), 'max', (1e6 * twtt_max_ref), 'value', (1e6 * twtt_max_ref))
        set(twtt_max_slide, 'min', (1e6 * twtt_min_ref), 'max', (1e6 * twtt_max_ref), 'value', (1e6 * twtt_min_ref))
        
        if trim_done
            pk.ind_trim_start ...
                            = tmp2;
        else
            trim_done       = true;
        end
        
        update_twtt_range
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
        if ~strcmp(disp_type, 'amp.')
            set(disp_group, 'selectedobject', disp_check(1))
            disp_type       = 'amp.';
            plot_db
        end
        
        if ispc
            tmp1            = '..\..\pk\';
        else
            tmp1            = '../../pk/';
        end
        
        % look for picks file in expected location
        try %#ok<TRYNC>
            if (~isempty(path_data) && exist([path_data tmp1 file_data(1:11)], 'dir'))
                if ispc
                    if isempty(path_pk)
                        path_pk = [path_data(1:strfind(path_data, '\block')) 'pk\' file_data(1:11) '/'];
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
        
        if ~isempty(file_pk)
            
            % check against data file
            try %#ok<TRYNC>
                if ~strcmp(file_data(1:(end - 4)), file_pk(1:(end - 7)))
                    set(status_box, 'string', ['Selected picks file (' file_pk(1:(end - 4)) ') might not match data file. Continue loading? Y: yes; otherwise: no...'])
                    waitforbuttonpress
                    if ~strcmpi(get(pkgui, 'currentcharacter'), 'Y')
                        set(status_box, 'string', 'Loading of picks file cancelled by user.')
                        return
                    end
                end
            end
            
            % load picks file
            tmp1            = load([path_pk file_pk]);
            try
                pk          = tmp1.pk;
                tmp1        = 0;
            catch
                set(status_box, 'string', 'Selected file does not contained a pk structure.')
                return
            end
            
            if ~pk.num_layer
                set(status_box, 'string', 'Picks file has no picks.')
                return
            end
            
            clear_plots
            [aresp_done, flat_done, keep_aresp_done, keep_phase_done, load_flat, match_done, mean_done, phase_done, pk_done, smooth_done] ...
                            = deal(false);
            if isfield(pk, 'poly_flat')
                load_flat   = true;
            else % sort layers by index
                tmp1        = NaN(pk.num_layer, block.num_trace);
                tmp2        = NaN(1, pk.num_layer);
                for ii = 1:pk.num_layer
                    pk.layer(ii).ind_y ...
                            = interp1(block.twtt, 1:num_sample_trim, pk.layer(ii).twtt, 'nearest', 'extrap');
                    pk.layer(ii).ind_y_smooth ...
                            = interp1(block.twtt, 1:num_sample_trim, pk.layer(ii).twtt_smooth, 'nearest', 'extrap');
                    pk.layer(ii).type ...
                            = 'max';
                    tmp1(ii, :) ...
                            = pk.layer(ii).twtt;
                    tmp2(ii)= nanmean(pk.layer(ii).ind_y);
                end
                [~, tmp3]   = sort(tmp2);
                pk.layer    = pk.layer(tmp3);
                if (length(find(all(~isnan(tmp1')))) > 3)
                    pk.keep_or_flat ...
                            = 'flat';
                else
                    pk.keep_or_flat ...
                            = 'keep';
                end
            end
            num_mean        = length(pk.ind_x_mean);
            
            set(freq_edit, 'string', num2str(1e-6 * pk.freq))
            set(num_ind_mean_edit, 'string', num2str(pk.num_ind_mean))
            set(num_win_edit, 'string', num2str(pk.num_win))
            set(length_smooth_edit, 'string', num2str(pk.length_smooth))
            set(twtt_match_edit, 'string', num2str(1e6 * pk.twtt_match))
            set(layer_list, 'string', num2cell(1:pk.num_layer), 'value', 1)
            
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
            pk_done         = true;
            smooth_done     = true(1, pk.num_layer);
            trim_y
            pause(time_pause)
            
            % change old variable names
            if isfield(pk, 'num_layer_keep')
                [pk.ind_keep_phase, pk.num_keep_phase] ...
                            = deal(pk.ind_layer_keep, pk.num_layer_keep);
                pk          = rmfield(pk, {'ind_layer_keep' 'num_layer_keep'});
                if isfield(pk, 'ind_x_start')
                    pk.ind_x_start_phase ...
                            = pk.ind_x_start;
                    pk      = rmfield(pk, 'ind_x_start');
                else
                    pk.ind_x_start_phase ...
                            = 1;
                end
                if isfield(pk, 'ind_y_start')
                    pk.ind_y_start_phase ...
                            = pk.ind_y_start;
                    pk      = rmfield(pk, 'ind_y_start');
                end
            end
            
            % load and plot phase tracking
            if (phase_avail && pk.num_keep_phase)
                prop_phase
                p_phase     = zeros(1, pk.num_phase);
                for ii = 1:pk.num_phase
                    p_phase(ii) = plot(block.dist_lin(ind_decim), (1e6 .* block.twtt(round(ind_y_phase(ii, ind_decim)))), 'b', 'linewidth', 1, 'visible', 'off'); % plot the phase-tracked layers
                end
                p_startphase= plot((ones(1, pk.num_phase) .* block.dist_lin(pk.ind_x_start_phase)), (1e6 .* block.twtt(pk.ind_y_start_phase)), 'm.', 'markersize', 12, 'visible', 'off'); 
                phase_done  = true;
                set(phase_check, 'value', 1)
                if pk.num_keep_phase
                    for ii = 1:pk.num_keep_phase
                        set(p_phase(pk.ind_keep_phase(ii)), 'linewidth', 2, 'color', 'w') % increase linewidth and change color
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
                [pk.ind_layer_keep_aresp, pk.num_layer_keep_aresp] ...
                            = deal(pk.ind_layer_keep_aresp, pk.num_layer_keep_aresp);
                pk          = rmfield(pk, {'ind_layer_keep_aresp' 'num_layer_keep_aresp'});
            end
            
            if (aresp_avail && isfield(pk, 'num_keep_aresp'))
                if pk.num_keep_aresp
                    prop_aresp
                    p_aresp = zeros(1, pk.num_aresp);
                    for ii = 1:pk.num_aresp
                        p_aresp(ii) = plot(block.dist_lin(ind_decim), (1e6 .* block.twtt(round(ind_y_aresp(ii, ind_decim)))), 'b', 'linewidth', 1, 'visible', 'off'); % plot the phase-tracked layers
                    end
                    p_startaresp ...
                            = plot((ones(1, pk.num_aresp) .* block.dist_lin(pk.ind_x_start_aresp)), (1e6 .* block.twtt(pk.ind_y_start_aresp)), 'm.', 'markersize', 12, 'visible', 'off'); % plot the starting y indices
                    aresp_done ...
                            = true;
                    set(aresp_check, 'value', 1)
                    if pk.num_keep_aresp
                        for ii = 1:pk.num_keep_aresp
                            set(p_aresp(pk.ind_keep_aresp(ii)), 'linewidth', 2, 'color', 'w') % increase linewidth and change color
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
                    set(status_box, 'string', 'Cannot display manually picked layers due to trimming issue...');
                    [pk.num_man, pk.ind_y_man] ...
                            = deal(0);
                    pause(time_pause)
                else
                    p_man   = zeros(pk.num_man, 2);
                    for ii = 1:pk.num_man
                        p_man(ii, 1) ...
                            = plot(block.dist_lin(ind_decim), (1e6 .* block.twtt(round(pk.ind_y_man(ii, ind_decim)))), 'linewidth', 2, 'color', [0.85 0.85 0.85], 'visible', 'off'); % max-picking spline
                        p_man(ii, 2) ...
                            = plot(block.dist_lin(ind_decim), (1e6 .* block.twtt(round(pk.ind_y_man(ii, ind_decim)))), 'linewidth', 2, 'color', [0.85 0.85 0.85], 'visible', 'off'); % max-picking spline
                    end
                    set(man_check, 'value', 1)
                end
            end
            
            % plot layers in flat space
            [p_flat, p_pk, p_pkflat, p_pkflatmark, p_pksmooth, p_pksmoothflat] ...
                            = deal(zeros(1, pk.num_layer));
            tmp1            = [];
            for ii = 1:pk.num_layer
                if ~isempty(ind_decim(~isnan(pk.layer(ii).ind_y(ind_decim))))
                    p_pk(ii)= plot(block.dist_lin(ind_decim(~isnan(pk.layer(ii).ind_y(ind_decim)))), (1e6 .* block.twtt(round(pk.layer(ii).ind_y(ind_decim(~isnan(pk.layer(ii).ind_y(ind_decim))))))), ...
                                   'r.', 'markersize', 12, 'visible', 'off');
                else
                    tmp1    = [tmp1 ii]; %#ok<AGROW>
                    continue
                end
                if ~isempty(ind_decim(~isnan(pk.layer(ii).ind_y_smooth(ind_decim))))
                    p_pksmooth(ii) = plot(block.dist_lin(ind_decim(~isnan(pk.layer(ii).ind_y_smooth(ind_decim)))), ...
                                        (1e6 .* block.twtt(round(pk.layer(ii).ind_y_smooth(ind_decim(~isnan(pk.layer(ii).ind_y_smooth(ind_decim))))))), 'g.', 'markersize', 12, 'visible', 'off');
                else
                    tmp1    = [tmp1 ii]; %#ok<AGROW>
                end
            end
            
            % get rid of empty layers
            if ~isempty(tmp1)
                tmp2        = setdiff(1:pk.num_layer, tmp1);
                delete(p_pk(tmp1))
                for ii = tmp1
                    if (logical(p_pksmooth(ii)) && ishandle(p_pksmooth(ii)))
                        delete(p_pksmooth(ii))
                    end
                    if (logical(p_pksmoothflat(ii)) && ishandle(p_pksmoothflat(ii)))
                        delete(p_pksmoothflat(tmp1(ii)))
                    end
                end
                curr_layer  = interp1(1:length(tmp2), tmp2, curr_layer, 'nearest', 'extrap');
                [pk.num_layer, pk.layer, p_pkflatmark, p_flat, p_pkflat, p_pk, smooth_done, p_pksmooth, p_pksmoothflat] ...
                            = deal(length(tmp2), pk.layer(tmp2), p_pkflatmark(tmp2), p_flat(tmp2), p_pkflat(tmp2), p_pk(tmp2), smooth_done(tmp2), p_pksmooth(tmp2), p_pksmoothflat(tmp2));
                set(layer_list, 'string', num2cell(1:pk.num_layer), 'value', curr_layer)
            end
            
            if load_flat % reflatten if possible
                
                set(status_box, 'string', 'Flattening amplitude data...')
                
                % load polynomial flattening
                ind_y_mat   = single((1:size(block.amp, 1))');
                ind_y_mat   = ind_y_mat(:, ones(1, block.num_trace)); % matrix of y indices
                if (size(pk.poly_flat, 1) == 3)
                    ind_y_flat = ((ind_y_mat .^ 2) .* pk.poly_flat(ones(num_sample_trim, 1), :)) + (ind_y_mat .* (pk.poly_flat((2 .* ones(num_sample_trim, 1)), :))) + pk.poly_flat((3 .* ones(num_sample_trim, 1)), :);
                else
                    ind_y_flat = ((ind_y_mat .^ 3) .* pk.poly_flat(ones(num_sample_trim, 1), :)) + ((ind_y_mat .^ 2) .* pk.poly_flat((2 .* ones(num_sample_trim, 1)), :)) + ...
                                 (ind_y_mat .* (pk.poly_flat((3 .* ones(num_sample_trim, 1)), :))) + pk.poly_flat((4 .* ones(num_sample_trim, 1)), :);
                end
                ind_y_mat   = 0;
                ind_y_flat(ind_y_flat < 1) ...
                            = 1;
                ind_y_flat(ind_y_flat > num_sample_trim) ...
                            = num_sample_trim;
                amp_flat    = NaN(size(block.amp, 1), block.num_trace, 'single');
                if parallel_check
                    tmp1    = block.amp;
                    parfor ii = 1:block.num_trace
                        amp_flat(:, ii) ...
                            = interp1(tmp1(:, ii), ind_y_flat(:, ii));
                    end
                    tmp1    = 0;
                else
                    for ii = 1:block.num_trace
                        amp_flat(:, ii) ...
                            = interp1(block.amp(:, ii), ind_y_flat(:, ii));
                    end
                end
                flat_done   = true;
                
                if ~isfield(pk.layer, 'ind_y_flat_mean') % hack for old merge files or updated ones
                    warning('off', 'MATLAB:interp1:NaNinY')
                    tmp1    = NaN(pk.num_layer, block.num_trace);
                    tmp2    = NaN(pk.num_layer, num_mean);
                    for ii = 1:pk.num_layer
                        tmp1(ii, :) = pk.layer(ii).ind_y;
                    end
                    for ii = 1:num_mean
                        [~, tmp3] ...
                            = unique(ind_y_flat(:, pk.ind_x_mean(ii)));
                        if any(~isnan(tmp1(:, pk.ind_x_mean(ii))))
                            tmp2(~isnan(tmp1(:, pk.ind_x_mean(ii))), ii) ...
                            = interp1(ind_y_flat(tmp3, pk.ind_x_mean(ii)), tmp3, tmp1(~isnan(tmp1(:, pk.ind_x_mean(ii))), pk.ind_x_mean(ii)), 'nearest', 'extrap');
                        end
                    end
                    
                    for ii = 1:pk.num_layer
                        pk.layer(ii).ind_y_flat_mean ...
                            = tmp2(ii, :); % get back to the pk structure
                    end
                    
                    [tmp1, tmp2] = deal(NaN(pk.num_layer, block.num_trace));
                    for ii = 1:pk.num_layer
                        tmp1(ii, :) = pk.layer(ii).ind_y_smooth;
                    end
                    for ii = 1:block.num_trace
                        [~, tmp3] ...
                            = unique(ind_y_flat(:, ii));
                        if any(~isnan(tmp1(:, ii)))
                            tmp2(~isnan(tmp1(:, ii)), ii) ...
                            = interp1(ind_y_flat(tmp3, ii), tmp3, tmp1(~isnan(tmp1(:, ii)), ii), 'nearest', 'extrap');
                        end
                    end
                    for ii = 1:pk.num_layer
                        pk.layer(ii).ind_y_flat_smooth ...
                            = tmp2(ii, :); % get back to the pk structure
                    end
                    warning('on', 'MATLAB:interp1:NaNinY')
                end
                
                % flatten surface pick
                if surf_avail
                    ind_surf_flat ...
                            = NaN(1, block.num_trace);
                    if parallel_check
                        parfor ii = 1:block.num_trace
                            [~, tmp2] = unique(ind_y_flat(:, ii));
                            if (length(tmp2) > 1)
                                ind_surf_flat(ii) = interp1(ind_y_flat(tmp2, ii), tmp2, ind_surf(ii), 'nearest', 'extrap'); %#ok<PFBNS>
                            end
                        end
                    else
                        for ii = 1:block.num_trace
                            [~, tmp2] = unique(ind_y_flat(:, ii));
                            if (length(tmp2) > 1)
                                ind_surf_flat(ii) = interp1(ind_y_flat(tmp2, ii), tmp2, ind_surf(ii), 'nearest', 'extrap');
                            end
                        end
                    end
                    ind_surf_flat ...
                            = round(ind_surf_flat);
                    ind_surf_flat((ind_surf_flat < 1) | (ind_surf_flat > num_sample_trim)) ...
                            = NaN;
                end
                
                % flatten bed pick
                if bed_avail
                    ind_bed_flat= NaN(1, block.num_trace);
                    if parallel_check
                        parfor ii = 1:block.num_trace
                            [~, tmp2] = unique(ind_y_flat(:, ii));
                            if (length(tmp2) > 1)
                                ind_bed_flat(ii) = interp1(ind_y_flat(tmp2, ii), tmp2, ind_bed(ii), 'nearest', 'extrap'); %#ok<PFBNS>
                            end
                        end
                    else
                        for ii = 1:block.num_trace
                            [~, tmp2] = unique(ind_y_flat(:, ii));
                            if (length(tmp2) > 1)
                                ind_bed_flat(ii) = interp1(ind_y_flat(tmp2, ii), tmp2, ind_bed(ii), 'nearest', 'extrap');
                            end
                        end
                    end
                    ind_bed_flat ...
                            = round(ind_bed_flat);
                    ind_bed_flat((ind_bed_flat < 1) | (ind_bed_flat > num_sample_trim)) ...
                            = NaN;
                end
                
                mean_done   = true;
                mean_flat
                pause(time_pause)
                
                tmp1        = NaN(pk.num_layer, block.num_trace);
                for ii = 1:pk.num_layer
                    tmp1(ii, :) ...
                            = pk.layer(ii).ind_y_smooth;
                end
                ind_x_pk    = sum(isnan(tmp1));
                ind_x_pk    = find((ind_x_pk == min(ind_x_pk)), 1); % find the first trace in the record that has the maximum number of layers
                tmp1        = interp1(pk.ind_x_mean, 1:num_mean, ind_x_pk, 'nearest', 'extrap');
                
                % plot layers
                for ii = 1:pk.num_layer
                    if ~isnan(pk.layer(ii).ind_y_smooth(ind_x_pk))
                        p_flat(ii) ...
                            = plot(block.dist_lin([1 block.num_trace]), (1e6 .* block.twtt(repmat(round(ind_y_flat(round(pk.layer(ii).ind_y_smooth(ind_x_pk)), ind_x_pk)), 1, 2))), 'w:', 'linewidth', 2, ...
                                   'visible', 'off');
                    else
                        p_flat(ii) ...
                            = plot(block.dist_lin([1 block.num_trace]), (1e6 .* block.twtt(repmat(round(ind_y_flat(round(pk.layer(ii).ind_y_smooth(find(~isnan(pk.layer(ii).ind_y_smooth), 1, 'first'))), ...
                                   find(~isnan(pk.layer(ii).ind_y_smooth), 1, 'first'))), 1, 2))), 'w:', 'linewidth', 2, 'visible', 'off');
                    end
                    p_pkflat(ii) ...
                            = plot(block.dist_lin(pk.ind_x_mean(~isnan(pk.layer(ii).ind_y_flat_mean))), (1e6 .* block.twtt(round(pk.layer(ii).ind_y_flat_mean(~isnan(pk.layer(ii).ind_y_flat_mean))))), ...
                                   'r.', 'markersize', 12, 'visible', 'off');
                    if ~isnan(pk.layer(ii).ind_y_flat_mean(tmp1))
                        p_pkflatmark(ii) ...
                            = plot(block.dist_lin(pk.ind_x_mean(tmp1)), (1e6 * block.twtt(pk.layer(ii).ind_y_flat_mean(tmp1))), 'ko', 'markersize', 12, 'markerfacecolor', 'y', 'visible', 'off');
                    else
                        p_pkflatmark(ii) ...
                            = plot(block.dist_lin(pk.ind_x_mean(find(~isnan(pk.layer(ii).ind_y_flat_mean), 1, 'first'))), ...
                                   (1e6 * block.twtt(pk.layer(ii).ind_y_flat_mean(find(~isnan(pk.layer(ii).ind_y_flat_mean), 1, 'first')))), 'ko', 'markersize', 12, 'markerfacecolor', 'y', 'visible', 'off');
                    end
                    p_pksmoothflat(ii) ...
                            = plot(block.dist_lin(ind_decim(~isnan(pk.layer(ii).ind_y_flat_smooth(ind_decim)))), ...
                                   (1e6 .* block.twtt(round(pk.layer(ii).ind_y_flat_smooth(ind_decim(~isnan(pk.layer(ii).ind_y_flat_smooth(ind_decim))))))), 'g.', 'markersize', 12, 'visible', 'off');
                end
                
                % flatten reference picks if those were already loaded
                if ref_done
                    p_refflat ...
                            = zeros(1, pk_ref.num_layer);
                    switch ref_start_end
                        case 'start'
                            for ii = 1:pk_ref.num_layer
                                try
                                    p_refflat(ii) = plot(block.dist_lin([1 block.ind_overlap(1)]), ...
                                        (1e6 .* block.twtt(repmat(round(ind_y_flat(round(pk_ref.layer(ii).ind_y_smooth(find(~isnan(pk_ref.layer(ii).ind_y_smooth), 1, 'last')) - pk.ind_trim_start + 1), ...
                                        (pk_ref.num_trace - pk_ref.ind_overlap(2) + 1))), 1, 2))), 'y:', 'linewidth', 2, 'visible', 'off');
                                catch
                                    continue
                                end
                            end
                        case 'end'
                            for ii = 1:pk_ref.num_layer
                                try
                                    p_refflat(ii) = plot(block.dist_lin([block.ind_overlap(2) end]), ...
                                                     (1e6 .* block.twtt(repmat(round(ind_y_flat(round(pk_ref.layer(ii).ind_y_smooth(find(~isnan(pk_ref.layer(ii).ind_y_smooth), 1, 'first')) - ...
                                                     pk.ind_trim_start + 1), (pk_ref.num_trace - pk_ref.ind_overlap(2) + 1))), 1, 2))), 'y:', 'linewidth', 2, 'visible', 'off');
                                catch
                                    continue
                                end
                            end
                    end
                end
                match_done  = true;
            end
            
            set([pk_check smooth_check], 'value', 1)
            show_phase
            show_aresp
            show_ref
            show_man
            show_smooth
            if load_flat
                set(status_box, 'string', ['Picks loaded from ' file_pk(1:(end - 4)) '.'])
            else
                set(status_box, 'string', ['Picks loaded from ' file_pk(1:(end - 4)) ' (no flattening).'])
            end
            
        else
            file_pk         = tmp1;
            set(status_box, 'string', 'No picks loaded.')
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
        
        ref_start_end       = '';
        
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
                    ref_start_end ...
                            = 'start';
                elseif (~exist([path_pk file_pk(1:(end - 9)) tmp1 '_pk.mat'], 'file') && exist([path_pk file_pk(1:(end - 9)) tmp2 '_pk.mat'], 'file'))
                    set(status_box, 'string', 'Only right-hand-side reference picks available. Load these? Y: yes; otherwise: no...')
                    waitforbuttonpress
                    if strcmpi(get(pkgui, 'currentcharacter'), 'Y')
                        file_ref ...
                            = [file_pk(1:(end - 9)) tmp2 '_pk.mat'];
                        ref_start_end ...
                            = 'end';
                    else
                        [file_ref, ref_start_end] ...
                            = deal('');
                    end
                elseif (exist([path_pk file_pk(1:(end - 9)) tmp1 '_pk.mat'], 'file') && exist([path_pk file_pk(1:(end - 9)) tmp2 '_pk.mat'], 'file'))
                    set(status_box, 'string', 'Left (L) or right (R) reference picks?')
                    waitforbuttonpress
                    if strcmpi(get(pkgui, 'currentcharacter'), 'L')
                        file_ref= [file_pk(1:(end - 9)) tmp1 '_pk.mat'];
                        ref_start_end ...
                                = 'start';
                    elseif strcmpi(get(pkgui, 'currentcharacter'), 'R')
                        file_ref= [file_pk(1:(end - 9)) tmp2 '_pk.mat'];
                        ref_start_end ...
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
        
        if ~isempty(file_ref)
            
            if strcmp(file_ref, file_pk)
                set(status_box, 'string', 'Reference picks file should not be same as current block''s pick file.')
                return
            end
            
            % load reference picks file
            tmp1            = load([path_pk file_ref]);
            try
                pk_ref      = tmp1.pk;
                tmp1        = 0;
            catch
                set(status_box, 'string', 'Selected file does not contain a pk structure.')
                return
            end
            
            if (any(p_ref) && any(ishandle(p_ref)))
                delete(p_ref(logical(p_ref) & ishandle(p_ref)))
            end
            if (any(p_refflat) && any(ishandle(p_refflat)))
                delete(p_refflat(logical(p_refflat) & ishandle(p_refflat)))
            end
            
            ref_done        = false;
            set(status_box, 'string', 'Loading reference picks...')
            
            % plot reference picks
            p_ref           = zeros(1, pk_ref.num_layer);
            
            if ~isempty(ref_start_end)
                switch ref_start_end
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
                    ref_start_end = 'start'; % reference picks from start of transect
                    for ii = 1:pk_ref.num_layer
                        p_ref(ii) = plot(block.dist_lin(1:decim:block.ind_overlap(1)), (1e6 .* pk_ref.layer(ii).twtt_smooth(pk_ref.ind_overlap(2):decim:end)), 'y.', 'markersize', 12, 'visible', 'off');
                    end
                end
            elseif ~isnan(pk_ref.ind_overlap(1))
                if (abs(pk_ref.dist(pk_ref.ind_overlap(1)) - block.dist(end)) < 0.01)
                    ref_start_end = 'end';
                    for ii = 1:pk_ref.num_layer
                        p_ref(ii) = plot(block.dist_lin(block.ind_overlap(2):decim:end), (1e6 .* pk_ref.layer(ii).twtt_smooth(1:decim:pk_ref.ind_overlap(1))), 'y.', 'markersize', 12, 'visible', 'off');
                    end
                end
            else
                set(status_box, 'string', 'This picks file does not overlap with the current block.')
                pk_ref      = struct;
                return
            end
            
            % plot flattened reference picks if flattening already done
            if flat_done
                p_refflat   = zeros(1, pk_ref.num_layer);
                switch ref_start_end
                    case 'start'
                        for ii = 1:pk_ref.num_layer
                            try
                                p_refflat(ii) = plot(block.dist_lin([1 block.ind_overlap(1)]), ...
                                    (1e6 .* block.twtt(repmat(round(ind_y_flat(round(pk_ref.layer(ii).ind_y_smooth(find(~isnan(pk_ref.layer(ii).ind_y_smooth), 1, 'last')) - pk.ind_trim_start + 1), ...
                                    (pk_ref.num_trace - pk_ref.ind_overlap(2) + 1))), 1, 2))), 'y:', 'linewidth', 2, 'visible', 'off');
                            catch
                                continue
                            end
                        end
                    case 'end'
                        for ii = 1:pk_ref.num_layer
                            try
                                p_refflat(ii) = plot(block.dist_lin([block.ind_overlap(2) end]), ...
                                                 (1e6 .* block.twtt(repmat(round(ind_y_flat(round(pk_ref.layer(ii).ind_y_smooth(find(~isnan(pk_ref.layer(ii).ind_y_smooth), 1, 'first')) - pk.ind_trim_start + 1), ...
                                                 (pk_ref.num_trace - pk_ref.ind_overlap(2) + 1))), 1, 2))), 'y:', 'linewidth', 2, 'visible', 'off');
                            catch
                                continue
                            end
                        end
                end
            end
            
            set(ref_check, 'value', 1)
            ref_done        = true;
            show_ref
            set(status_box, 'string', ['Reference picks loaded from ' file_ref(1:(end - 4)) '.'])
            
        else
            set(status_box, 'string', 'No reference picks loaded.')
        end
    end

%% Automatically track layers using the horizontal phase gradient

    function track_phase(source, eventdata)

        if ~load_done
            set(status_box, 'string', 'Load data first.')
            return
        end
        if ~trim_done
            set(status_box, 'string', 'Trim excess data first.')
            return
        end
        if ~phase_avail
            set(status_box, 'string', 'Horizontal phase gradient not available.')
            return
        end
        
        if (logical(p_startphase) && ishandle(p_startphase));
            delete(p_startphase)
        end
        if (any(p_phase) && any(ishandle(p_phase)))
            delete(p_phase(logical(p_phase) & ishandle(p_phase)))
        end
        
        if keep_phase_done
            [pk.num_keep_phase, pk.ind_keep_phase] ...
                            = deal(0);
        end
        [phase_done, keep_phase_done] ...
                            = deal(false);
        
        set(phase_check, 'value', 0)
        set(status_box, 'string', 'Choose trace to propagate layers from...(Q: cancel)')
        tic
        ii = 0;
        
        set(pkgui, 'keypressfcn', [])
        while true
            
            [ind_x_pk, ind_y_pk, button] ...
                            = ginput(1); % get x/y indices to start tracing from (e.g., thickest ice)
            if ~ii
                ii          = 1;
                continue
            end
            
            if strcmpi(char(button), 'Q')
                
                set(status_box, 'string', 'Cancelled tracking phase.')
                set(pkgui, 'keypressfcn', @keypress)
                return
                
            elseif (button == 1) % left-click
                
                set(status_box, 'string', 'Propagating phase...')
                
                pk.ind_x_start_phase ...
                            = interp1(block.dist_lin(ind_decim), ind_decim, ind_x_pk, 'nearest', 'extrap'); % convert to index
                pk.ind_y_phase_max = interp1(block.twtt, 1:num_sample_trim, (1e-6 * ind_y_pk), 'nearest', 'extrap'); % maximum index to pick, not much in data below here for layer tracing; used in flattening section
                if (block.twtt_surf(pk.ind_x_start_phase) > twtt_min)
                    pk.ind_y_start_phase = (interp1(block.twtt, 1:num_sample_trim, block.twtt_surf(pk.ind_x_start_phase), 'nearest', 'extrap') + int_track):int_track:pk.ind_y_phase_max; ...
                        % y indices from which to propagate phase-tracked layers
                else
                    pk.ind_y_start_phase = (1 + int_track):int_track:pk.ind_y_phase_max; % y indices from which to propagate first attempt at phase-tracked layers
                end
                pk.num_phase    = length(pk.ind_y_start_phase); % number of y indices to test in first attempt
                
                prop_phase % propagate phase in a separate sub-function
                
                axes(ax_radar) %#ok<*LAXES>
                p_phase         = zeros(1, pk.num_phase);
                for ii = 1:pk.num_phase
                    p_phase(ii) = plot(block.dist_lin(ind_decim), (1e6 .* block.twtt(round(ind_y_phase(ii, ind_decim)))), 'b', 'linewidth', 1); % plot the phase-tracked layers
                end
                p_startphase    = plot((ones(1, pk.num_phase) .* block.dist_lin(pk.ind_x_start_phase)), (1e6 .* block.twtt(pk.ind_y_start_phase)), 'm.', 'markersize', 12); % plot the starting y indices
                set(phase_check, 'value', 1)
                phase_done      = true;
                set(status_box, 'string', ['Traced ' num2str(pk.num_phase) ' layers starting from ' num2str(block.dist(pk.ind_x_start_phase), '%3.1f') ' km in ' num2str(toc, '%.0f') ' s.'])
                break
                
            end
        end
        set(pkgui, 'keypressfcn', @keypress)
    end

%% Propagate phase starts

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
        
        tmp1                 = 0;
        
    end

%% Pick phase-propagated layers to keep for polynomial fitting

    function pk_keep_phase(source, eventdata)
        
        if ~phase_done
            set(status_box, 'string', 'No phase-tracked layers yet.')
            return
        end
        if ~strcmp(disp_type, 'amp.')
            set(disp_group, 'selectedobject', disp_check(1))
            plot_db
        end
        if ~get(phase_check, 'value')
            set(phase_check, 'value', 1)
            show_phase
            keep_phase_done = false;
        end
        
        set(status_box, 'string', 'Pick best phase-traced layers for flattening...')
        pause(time_pause)
        set(status_box, 'string', 'Return = keep; d = unkeep; q = stop; otherwise = pick again.')
        
        axes(ax_radar)
        ii                  = 0;
        
        set(pkgui, 'keypressfcn', [])
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
                    ind_y_pk= interp1(block.twtt, 1:num_sample_trim, (1e-6 * ind_y_pk), 'nearest', 'extrap');
                    tmp1    = pk.ind_keep_phase(~isnan(ind_y_phase(pk.ind_keep_phase, ind_x_pk)));
                    if (length(tmp1) > 1)
                        tmp1= interp1(ind_y_phase(tmp1, ind_x_pk), tmp1, ind_y_pk, 'nearest', 'extrap');
                    end
                    set(p_phase(tmp1), 'color', 'b', 'linewidth', 1)
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
                    ind_y_pk= interp1(tmp1, tmp2(tmp3), interp1(block.twtt, 1:num_sample_trim, (1e-6 * ind_y_pk), 'nearest', 'extrap'), 'nearest', 'extrap'); % index of nearest layer
                else
                    ind_y_pk= tmp2;
                end
                if any(ind_y_pk == pk.ind_keep_phase)
                    set(status_box, 'string', 'Chosen layer already kept. Pick again...')
                    continue
                end
                set(p_phase(ind_y_pk), 'color', 'w') % change layer color to white
                
                waitforbuttonpress
                
                if (double(get(pkgui, 'currentcharacter')) == 13)
                    set(p_phase(ind_y_pk), 'linewidth', 2) % increase linewidth
                    pk.ind_keep_phase(pk.num_keep_phase + 1) = ind_y_pk; % preserve layer number
                    pk.num_keep_phase = pk.num_keep_phase + 1;
                else
                    set(p_phase(ind_y_pk), 'color', 'b') % return to previous color
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
        
        set(pkgui, 'keypressfcn', @keypress)
    end

%% Automatically track layers using ARESP

    function track_aresp(source, eventdata)

        if ~load_done
            set(status_box, 'string', 'Load data first.')
            return
        end
        if ~trim_done
            set(status_box, 'string', 'Trim excess data first.')
            return
        end
        if ~aresp_avail
            set(status_box, 'string', 'ARESP not available.')
            return
        end
        
        if (logical(p_startaresp) && ishandle(p_startaresp))
            delete(p_startaresp)
        end
        if (any(p_aresp) && any(ishandle(p_aresp)))
            delete(p_aresp(logical(p_aresp) & ishandle(p_aresp)))
        end
        
        if keep_aresp_done
            [pk.num_keep_aresp, pk.ind_keep_aresp] ...
                            = deal(0);
        end
        [aresp_done, keep_aresp_done] ...
                            = deal(false);
        
        set(aresp_check, 'value', 0)
        set(status_box, 'string', 'Choose trace to propagate layers from...(Q: cancel)')
        tic
        ii = 0;
        
        set(pkgui, 'keypressfcn', [])
        while true
            
            [ind_x_pk, ind_y_pk, button] ...
                            = ginput(1); % get x index to start tracing from (e.g., thickest ice)
            if ~ii
                ii          = 1;
                continue
            end
            
            if strcmpi(char(button), 'Q')
                
                set(status_box, 'string', 'Cancelled tracking ARESP.')
                set(pkgui, 'keypressfcn', @keypress)
                return
                
            elseif (button == 1) % left-click
                
                set(status_box, 'string', 'Tracking ARESP slopes...')
                
                pk.ind_x_start_aresp ...
                            = interp1(block.dist_lin(ind_decim), ind_decim, ind_x_pk, 'nearest', 'extrap'); % convert to index
                pk.ind_y_aresp_max ...
                            = interp1(block.twtt, 1:num_sample_trim, (1e-6 * ind_y_pk), 'nearest', 'extrap'); % maximum index to pick, not much in data below here for layer tracing; used in flattening section
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
                
                axes(ax_radar) %#ok<*LAXES>
                p_aresp     = zeros(1, pk.num_aresp);
                for ii = 1:pk.num_aresp
                    p_aresp(ii) = plot(block.dist_lin(ind_decim(~isnan(ind_y_aresp(ii, ind_decim)))), (1e6 .* block.twtt(round(ind_y_aresp(ii, ind_decim(~isnan(ind_y_aresp(ii, ind_decim))))))), 'c', 'linewidth', 1);
                end
                p_startaresp= plot((ones(1, pk.num_aresp) .* block.dist_lin(pk.ind_x_start_aresp)), (1e6 .* block.twtt(pk.ind_y_start_aresp)), 'r.', 'markersize', 12); % plot the starting y indices
                set(aresp_check, 'value', 1)
                aresp_done  = true;
                set(status_box, 'string', ['Traced ' num2str(pk.num_aresp) ' layers starting from ' num2str(block.dist(pk.ind_x_start_aresp), '%3.1f') ' km in ' num2str(toc, '%.0f') ' s.'])
                break
            end
        end
        set(pkgui, 'keypressfcn', @keypress)
    end

%% Propagate ARESP starts

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

%% Pick ARESP-propagated layers to keep for polynomial fitting

    function pk_keep_aresp(source, eventdata)
        
        if ~aresp_done
            set(status_box, 'string', 'No ARESP-tracked layers yet.')
            return
        end
        if ~strcmp(disp_type, 'amp.')
            set(disp_group, 'selectedobject', disp_check(1))
            plot_db
        end
        if ~get(aresp_check, 'value')
            set(aresp_check, 'value', 1)
            show_aresp
            keep_aresp_done = false;
        end
        
        set(status_box, 'string', 'Pick best ARESP-traced layers for flattening...')
        pause(time_pause)
        set(status_box, 'string', 'Return = keep; d = unkeep; q = stop; otherwise = pick again.')
        
        axes(ax_radar)
        ii                  = 0;
        
        set(pkgui, 'keypressfcn', [])
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
                    ind_y_pk= interp1(block.twtt, 1:num_sample_trim, (1e-6 * ind_y_pk), 'nearest', 'extrap');
                    tmp1    = pk.ind_keep_aresp(~isnan(ind_y_aresp(pk.ind_keep_aresp, ind_x_pk)));
                    if (length(tmp1) > 1)
                        tmp1= interp1(ind_y_aresp(tmp1, ind_x_pk), tmp1, ind_y_pk, 'nearest', 'extrap');
                    end
                    set(p_aresp(tmp1), 'color', 'c', 'linewidth', 1)
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
                    ind_y_pk= interp1(tmp1, tmp2(tmp3), interp1(block.twtt, 1:num_sample_trim, (1e-6 * ind_y_pk), 'nearest', 'extrap'), 'nearest', 'extrap'); % index of nearest layer
                else
                    ind_y_pk= tmp2;
                end
                if any(ind_y_pk == pk.ind_keep_aresp)
                    set(status_box, 'string', 'Chosen layer already kept. Pick again...')
                    continue
                end
                set(p_aresp(ind_y_pk), 'color', 'w') % change layer color to white
                
                waitforbuttonpress
                
                if (double(get(pkgui, 'currentcharacter')) == 13)
                    set(p_aresp(ind_y_pk), 'linewidth', 2) % increase linewidth
                    pk.ind_keep_aresp(pk.num_keep_aresp + 1) ...
                            = ind_y_pk; % preserve layer number
                    pk.num_keep_aresp ...
                            = pk.num_keep_aresp + 1;
                else
                    set(p_aresp(ind_y_pk), 'color', 'c') % return to previous color
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
        
        set(pkgui, 'keypressfcn', @keypress)
    end

%% Manually pick a layer to use for prediction

    function track_man(source, eventdata)
        
        if ~trim_done
            set(status_box, 'string', 'Trim excess data first.')
            return
        end
        if ~strcmp(disp_type, 'amp.')
            set(disp_group, 'selectedobject', disp_check(1))
            disp_type       = 'amp.';
            plot_db
        end
        
        axes(ax_radar)
        ii                  = 0;
        [ind_x_pk, ind_y_pk, button] ...
                            = deal([]);
        if ~p_man
            p_man           = zeros(0, 2);
        end
        
        tmp4                = 0;
        
        set(pkgui, 'keypressfcn', [])
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
                ind_y_pk(ii)= interp1(block.twtt, 1:num_sample_trim, (1e-6 .* ind_y_pk(ii)), 'nearest', 'extrap');
                if (ii > 1)
                    delete(p_man((pk.num_man + 1), 1))
                end
                p_man((pk.num_man + 1), 1) ...
                            = plot(block.dist_lin(ind_x_pk), (1e6 .* block.twtt(ind_y_pk)), 'wx', 'markersize', 12); % original picks
                ii          = ii + 1;
                tmp4        = button;
                
            elseif (strcmpi(char(button), 'U') && (ii > 1))
                
                [ind_x_pk, ind_y_pk] ...
                            = deal(ind_x_pk(1:(end - 1)), ind_y_pk(1:(end - 1)));
                delete(p_man((pk.num_man + 1), 1))
                if (ii > 2)
                    p_man((pk.num_man + 1), 1) ...
                            = plot(block.dist_lin(ind_x_pk), (1e6 .* block.twtt(ind_y_pk)), 'wx', 'markersize', 12); % original picks
                end
                ii          = ii - 1;
                tmp4        = button;
                
            elseif strcmpi(char(button), 'D') % delete a manual layer
                
                if (tmp4 == 1)
                    set(status_box, 'string', 'Cannot delete a manual layer in the middle of picking one...')
                    pause(time_pause)
                    continue
                end
                if pk.num_man
                    tmp1    = interp1(block.dist_lin(ind_decim), ind_decim, tmp1, 'nearest', 'extrap'); % raw picks must be indices, not dimensionalized vectors (horizontal)
                    tmp2    = interp1(block.twtt, 1:num_sample_trim, (1e-6 .* tmp2), 'nearest', 'extrap');
                    tmp3    = NaN(1, pk.num_man);
                    for ii = 1:pk.num_man
                        tmp3(ii) ...
                            = pk.ind_y_man(ii, tmp1);
                    end
                    if (length(tmp3) > 1)
                        tmp3= interp1(tmp3, 1:pk.num_man, tmp2, 'nearest', 'extrap');
                    else
                        tmp3= 1;
                    end
                    pk.ind_y_man ...
                            = pk.ind_y_man(setdiff(1:pk.num_man, tmp3), :);
                    delete(p_man(tmp3, :))
                    p_man   = p_man(setdiff(1:pk.num_man, tmp3), :);
                    pk.num_man ...
                            = pk.num_man - 1;
                    set(status_box, 'string', ['Deleted manual layer #' num2str(tmp3) '.'])
                    pause(time_pause)
                end
                
            elseif (double(get(pkgui, 'currentcharacter')) == 13)
                
                if (length(ind_x_pk) < 3)
                    if (any(ishandle(p_man(end, :))) && any(p_man(end, :)))
                        delete(p_man(end, (ishandle(p_man(end, :)) & logical(p_man(end, :)))))
                    end
                    p_man   = p_man(1:(end - 1), :);
                    set(status_box, 'string', 'Not enough picked points to make a manual layer. Start over.')
                    pause(time_pause)
                    ii      = 1;
                    continue
                end
                pk.num_man  = pk.num_man + 1;
                
                if ~issorted(ind_x_pk) % resort picks
                    [ind_x_pk, tmp1]= sort(ind_x_pk);
                    ind_y_pk= ind_y_pk(tmp1);
                end
                if (length(unique(ind_x_pk)) < length(ind_x_pk)) % don't keep picks that are accidentally at the same horizontal index, otherwise spline will fail
                    [ind_x_pk, tmp1]= unique(ind_x_pk);
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
                
                pk.ind_y_man(pk.num_man, 1:block.num_trace) ...
                            = spline(ind_x_pk, ind_y_pk, 1:block.num_trace); % interpolate spline through picks
                
                p_man(pk.num_man, 2)= plot(block.dist_lin(ind_decim), (1e6 .* block.twtt(round(pk.ind_y_man(pk.num_man, ind_decim)))), 'linewidth', 2, 'color', [0.85 0.85 0.85]); % max-picking spline
                
                set(man_check, 'value', 1)
                show_man
                ii          = 1;
                [ind_x_pk, ind_y_pk] ...
                            = deal([]);
                set(status_box, 'string', 'Manual layer successfully picked. Now on to another...')
                pause(time_pause)
                tmp4        = button;
                
            elseif strcmpi(char(button), 'Q')
                
                if (size(p_man, 1) > pk.num_man)
                    if (any(p_man(end, :)) && any(ishandle(p_man(end, :))))
                        delete(p_man(end, (logical(p_man(end, :)) & ishandle(p_man(end, :)))))
                    end
                    p_man   = p_man(1:(end - 1), :);
                end
                set(status_box, 'string', 'Ended manual picking.')
                set(pkgui, 'keypressfcn', @keypress)
                return                
            end
        end
    end

%% Flatten amplitudes using layers

    function flatten(source, eventdata)
        
        if pk.num_layer
            tmp1            = NaN(pk.num_layer, block.num_trace);
            for ii = 1:pk.num_layer
                tmp1(ii, :) = pk.layer(ii).ind_y_smooth;
            end
        else
            tmp1            = NaN;
        end
        
        if (strcmp(pk.keep_or_flat, 'keep') && ((pk.num_keep_phase + pk.num_man + pk.num_keep_aresp + length(find(all(~isnan(tmp1'))))) < 3) && ~load_flat)
            set(status_box, 'string', 'Not enough approximate layers to flatten (need 3+ including surface).')
            return
        end
        if strcmp(pk.keep_or_flat, 'flat')
            if ~all(smooth_done)
                set(status_box, 'string', 'Not all layers have been smoothed prior to flattening.')
                return
            end
            if (pk.num_layer < 3)
                set(status_box, 'string', 'Not enough layers to attempt flattening.')
                return
            end
        end
        
        if (any(p_flat) && any(ishandle(p_flat)))
            delete(p_flat(logical(p_flat) & ishandle(p_flat)))
        end
        if (logical(p_surfflat) && ishandle(p_surfflat))
            delete(p_surfflat)
        end
        
        [flat_done, pk_done]= deal(false);
        
        set(status_box, 'string', 'Starting polynomial fits...')
        pause(time_pause)
        
        tic
        
        pk.poly_flat       = NaN(4, block.num_trace, 'single');
        % 2nd-order polynomial fits between kept layers at pk.ind_x_start (x) and kept layers at each trace (y)
        
        if parallel_check
            pctRunOnAll warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale')
            pctRunOnAll warning('off', 'MATLAB:polyfit:PolyNotUnique')
        else
            warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale')
            warning('off', 'MATLAB:polyfit:PolyNotUnique')
        end
        
        switch pk.keep_or_flat
            
            case 'keep'
                
                % compile layers to use in flattening, then sort them all together
                if surf_avail
                    tmp2    = ind_surf;
                else
                    tmp2    = [];
                end
                if pk.num_keep_phase
                    tmp2    = [tmp2; ind_y_phase(pk.ind_keep_phase, :)];
                end
                if pk.num_keep_aresp
                    tmp2    = [tmp2; ind_y_aresp(pk.ind_keep_aresp, :)];
                end
                if pk.num_man
                    tmp2    = [tmp2; pk.ind_y_man];
                end
                if any(pk_tmp) % earlier layers picked in flattened space
                    tmp2    = [tmp2; pk_tmp];
                    pk_tmp  = 0;
                elseif (~load_flat && pk.num_layer)
                    if ~isempty(find(all(~isnan(tmp1')), 1))
                        tmp2 = [tmp2; tmp1(all(~isnan(tmp1')), :)];
                    end
                end
                
                [~, tmp1]   = sort(nanmean(tmp2, 2));
                tmp2        = tmp2(tmp1, :); % all layers sorted by mean value
                
                % now keep the maximum number of layers
                ind_x_pk    = sum(isnan(tmp2));
                ind_x_pk    = find((ind_x_pk == min(ind_x_pk)), 1);
                ind_y_pk    = tmp2(:, ind_x_pk);
                tmp2        = tmp2(~isnan(ind_y_pk), :);
                tmp1        = tmp2(:, ind_x_pk);
                tmp4        = pk.poly_flat;
                
                % now polyfit
                if parallel_check
                    parfor ii = 1:block.num_trace
                        tmp4(:, ii) ...
                            = polyfit(tmp1(~isnan(tmp2(:, ii))), tmp2(~isnan(tmp2(:, ii)), ii), 3)'; %#ok<PFBNS>
                    end
                else
                    for ii = 1:block.num_trace
                        tmp4(:, ii) ...
                            = polyfit(tmp1(~isnan(tmp2(:, ii))), tmp2(~isnan(tmp2(:, ii)), ii), 3)';
                    end
                end
                pk.poly_flat= tmp4;
                tmp4        = 0;
                
                tmp1        = round((2 * pk.length_smooth) / nanmean(diff(block.dist)));
                if parallel_check
                    tmp2    = pk.poly_flat;
                    parfor ii = 1:4
                        tmp2(ii, :) = smooth_lowess(tmp2(ii, :), tmp1);
                    end
                    pk.poly_flat = tmp2;
                    tmp2 = 0;
                else
                    for ii = 1:4
                        pk.poly_flat(ii, :) = smooth_lowess(pk.poly_flat(ii, :), tmp1);
                    end
                end
                
                
            case 'flat' % more complicated if re-doing based on layers picked in flattened projection
                
                tmp2        = NaN(pk.num_layer, block.num_trace);
                for ii = 1:pk.num_layer
                    tmp2(ii, :) = pk.layer(ii).ind_y_smooth;
                end
                if surf_avail
                    tmp2    = [ind_surf; tmp2];
                end
                
                ind_x_pk    = sum(isnan(tmp2));
                ind_x_pk    = find((ind_x_pk == min(ind_x_pk)), 1); % find the first trace in the record that has the maximum number of layers, i.e., minimum number of NaNs
                ind_y_pk    = tmp2(:, ind_x_pk);
                tmp3        = interp1(pk.ind_x_mean, 1:num_mean, ind_x_pk, 'nearest', 'extrap'); % index in pk.ind_x_mean where ind_x_pk is
                
                % can't work with NaN's
                pk.num_layer= length(find(~isnan(ind_y_pk))); % reduce number of layers because they don't exist at overlapping points
                tmp2        = tmp2(~isnan(ind_y_pk), :);
                ind_y_pk    = ind_y_pk(~isnan(ind_y_pk));
                
                % polynomial fitting, including surface
                if (length(find(all(~isnan(tmp2')))) < 4) % not enough best flattened layers for polyfit at some traces
                    pk_tmp  = tmp2(2:end, :);
                    pk.keep_or_flat ...
                            = 'keep';
                    set(status_box, 'string', 'Not enough flat layers. Re-run flattening with combined approach.')
                    return
                else
                    tmp4    = pk.poly_flat;
                    tmp1    = tmp2(:, ind_x_pk);
                    if parallel_check
                        parfor ii = 1:block.num_trace
                            tmp4(:, ii) ...
                                = polyfit(tmp1(~isnan(tmp2(:, ii))), tmp2(~isnan(tmp2(:, ii)), ii), 3)'; %#ok<PFBNS>
                        end
                    else
                        for ii = 1:block.num_trace
                            tmp4(:, ii) ...
                                = polyfit(tmp1(~isnan(tmp2(:, ii))), tmp2(~isnan(tmp2(:, ii)), ii), 3)';
                        end
                    end
                end
                pk.poly_flat= tmp4;
                if surf_avail
                    tmp2    = tmp2(2:end, :);
                end
                tmp4        = 0;
                
                % smooth out the polynomials if there are any gaps in the layers, which create awkward steps in the flattened space
                if any(isnan(tmp2(:)))
                    tmp1     = round((2 * pk.length_smooth) / nanmean(diff(block.dist)));
                    if parallel_check
                        tmp4 = pk.poly_flat;
                        parfor ii = 1:4
                            tmp4(ii, :) = smooth_lowess(tmp4(ii, :), tmp1);
                        end
                        pk.poly_flat = tmp4;
                        tmp4 = 0;
                    else
                        for ii = 1:4
                            pk.poly_flat(ii, :) = smooth_lowess(pk.poly_flat(ii, :), tmp1);
                        end
                    end
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
        ind_y_flat          = ((ind_y_mat .^ 3) .* pk.poly_flat(ones(num_sample_trim, 1), :)) + ((ind_y_mat .^ 2) .* pk.poly_flat((2 .* ones(num_sample_trim, 1)), :)) + ...
                              (ind_y_mat .* (pk.poly_flat((3 .* ones(num_sample_trim, 1)), :))) + pk.poly_flat((4 .* ones(num_sample_trim, 1)), :);
        ind_y_mat           = 0;
        ind_y_flat(ind_y_flat < 1) ...
                            = 1;
        ind_y_flat(ind_y_flat > num_sample_trim) ...
                            = num_sample_trim;
        set(status_box, 'string', ['Fit layers and calculated remapping in ' num2str(toc, '%.0f') ' s...'])
        pause(time_pause)
        
        tic
        
        % flattened radargram based on the fits to the kept layers
        amp_flat            = NaN(size(block.amp, 1), block.num_trace, 'single');
        if parallel_check
            pctRunOnAll warning('off', 'MATLAB:interp1:NaNinY')
            tmp1            = block.amp;
            parfor ii = 1:block.num_trace
                amp_flat(:, ii) = interp1(tmp1(:, ii), ind_y_flat(:, ii));
            end
            tmp1            = 0;
            pctRunOnAll warning('on', 'MATLAB:interp1:NaNinY')
        else
            warning('off', 'MATLAB:interp1:NaNinY')
            for ii = 1:block.num_trace
                amp_flat(:, ii) = interp1(block.amp(:, ii), ind_y_flat(:, ii));
            end
            warning('on', 'MATLAB:interp1:NaNinY')
        end
        
        set(status_box, 'string', ['Interpolated flattened layers in ' num2str(toc, '%.0f') ' s...'])
        
        % flattened radargram and the kept layers (horizontal when flattened)
        flat_switch         = 'full';
        flat_done           = true;
        smooth_done         = false(1, pk.num_layer);
        
        ind_surf_flat       = NaN(1, block.num_trace);
        if surf_avail % flatten surface pick
            if parallel_check
                parfor ii = 1:block.num_trace
                    [~, tmp4] = unique(ind_y_flat(:, ii));
                    if (length(tmp4) > 1)
                        ind_surf_flat(ii) = interp1(ind_y_flat(tmp4, ii), tmp4, ind_surf(ii), 'nearest', 'extrap'); %#ok<PFBNS>
                    end
                end
            else
                for ii = 1:block.num_trace
                    [~, tmp4]   = unique(ind_y_flat(:, ii));
                    if (length(tmp4) > 1)                 
                        ind_surf_flat(ii) = interp1(ind_y_flat(tmp4, ii), tmp4, ind_surf(ii), 'nearest', 'extrap');
                    end
                end
            end
            ind_surf_flat   = round(ind_surf_flat);
            ind_surf_flat((ind_surf_flat < 1) | (ind_surf_flat > num_sample_trim)) ...
                            = NaN;
        end
        
        if bed_avail % flatten bed pick
            ind_bed_flat    = NaN(1, block.num_trace);
            if parallel_check
                parfor ii = 1:block.num_trace
                    [~, tmp4] = unique(ind_y_flat(:, ii));
                    if (length(tmp4) > 1)
                        ind_bed_flat(ii) = interp1(ind_y_flat(tmp4, ii), tmp4, ind_bed(ii), 'nearest', 'extrap'); %#ok<PFBNS>
                    end
                end
            else
                for ii = 1:block.num_trace
                    [~, tmp4] = unique(ind_y_flat(:, ii));
                    if (length(tmp4) > 1)
                        ind_bed_flat(ii) = interp1(ind_y_flat(tmp4, ii), tmp4, ind_bed(ii), 'nearest', 'extrap');
                    end
                end
            end
            ind_bed_flat    = round(ind_bed_flat);
            ind_bed_flat((ind_bed_flat < 1) | (ind_bed_flat > num_sample_trim)) ...
                            = NaN;
        end
        
        switch pk.keep_or_flat
            
            case 'keep'
                
                p_flat      = zeros(1, (pk.num_keep_phase + pk.num_keep_aresp + pk.num_man));
                
                if pk.num_keep_phase
                    for ii = 1:pk.num_keep_phase
                        p_flat(ii) = plot(block.dist_lin([1 block.num_trace]), (1e6 .* block.twtt(repmat(round(ind_y_flat(round(ind_y_phase(pk.ind_keep_phase(ii), ind_x_pk)), ind_x_pk)), 1, 2))), 'w:', ...
                                          'linewidth', 2);
                    end
                end
                if pk.num_keep_aresp
                    for ii = (pk.num_keep_phase + 1):(pk.num_keep_phase + pk.num_keep_aresp)
                        p_flat(ii) = plot(block.dist_lin([1 block.num_trace]), (1e6 .* block.twtt(repmat(round(ind_y_flat(round(ind_y_aresp(pk.ind_keep_aresp(ii), ind_x_pk)), ind_x_pk)), 1, 2))), 'w:', ...
                                         'linewidth', 2);
                    end
                end
                if pk.num_man
                    for ii = (pk.num_keep_phase + pk.num_keep_aresp + 1):(pk.num_keep_phase + pk.num_keep_aresp + pk.num_man)
                        p_flat(ii) = plot(block.dist_lin([1 block.num_trace]), ...
                                          (1e6 .* block.twtt(repmat(round(ind_y_flat(round(pk.ind_y_man((ii - pk.num_keep_phase - pk.num_keep_aresp), ind_x_pk)), ind_x_pk)), 1, 2))), 'w:', 'linewidth', 2);
                    end
                end
                
                if ~load_flat
                    
                    [p_pkflat, p_pksmoothflat] ...
                            = deal(zeros(1, pk.num_layer));
                    
                    % flatten ind_y and ind_y_smooth if they were recently loaded from an old format merge
                    for ii = 1:pk.num_layer
                        pk.layer(ii).ind_y_flat_mean    = NaN(1, num_mean);
                        pk.layer(ii).ind_y_flat_smooth  = NaN(1, block.num_trace);
                        for jj = 1:num_mean
                            [~, tmp3]                           = unique(ind_y_flat(:, pk.ind_x_mean(jj)));
                            pk.layer(ii).ind_y_flat_mean(jj)    = interp1(ind_y_flat(tmp3, pk.ind_x_mean(jj)), tmp3, pk.layer(ii).ind_y(pk.ind_x_mean(jj)), 'nearest', 'extrap');
                        end
                        for jj = 1:block.num_trace
                            [~, tmp3]                           = unique(ind_y_flat(:, jj));
                            pk.layer(ii).ind_y_flat_smooth(jj)  = interp1(ind_y_flat(tmp3, jj), tmp3, pk.layer(ii).ind_y_smooth(jj), 'nearest', 'extrap');
                        end
                        [pk.layer(ii).ind_y_flat_mean, pk.layer(ii).ind_y_flat_smooth] ...
                            = deal(round(pk.layer(ii).ind_y_flat_mean), round(pk.layer(ii).ind_y_flat_smooth));
                        pk.layer(ii).ind_y_flat_mean((pk.layer(ii).ind_y_flat_mean < 1) | (pk.layer(ii).ind_y_flat_mean > num_sample_trim)) ...
                            = NaN;
                        pk.layer(ii).ind_y_flat_smooth((pk.layer(ii).ind_y_flat_smooth < 1) | (pk.layer(ii).ind_y_flat_smooth > num_sample_trim)) ...
                            = NaN;
                        p_pkflat(ii)       = plot(block.dist_lin(pk.ind_x_mean(~isnan(pk.layer(ii).ind_y_flat_mean))), ...
                                                  (1e6 .* block.twtt(round(pk.layer(ii).ind_y_flat_mean(~isnan(pk.layer(ii).ind_y_flat_mean))))), 'r.', 'markersize', 12, 'visible', 'off');
                        p_pksmoothflat(ii) = plot(block.dist_lin(ind_decim(~isnan(pk.layer(ii).ind_y_flat_smooth(ind_decim)))), ...
                                                  (1e6 .* block.twtt(round(pk.layer(ii).ind_y_flat_smooth(ind_decim(~isnan(pk.layer(ii).ind_y_flat_smooth(ind_decim))))))), 'g.', 'markersize', 12, 'visible', 'off');
                    end
                    
                    pk_done = true;
                    smooth_done ...
                            = true(1, pk.num_layer);
                    
                elseif pk.num_layer % delete all layers
                    del_all
                    set(status_box, 'string', 'Deleted limited picked layers.')
                end
                
                set(disp_group, 'selectedobject', disp_check(5))
                
                mean_done   = false;
                mean_flat
                pause(time_pause)
                
            case 'flat'
                
                if (any(p_pkflat) && any(ishandle(p_pkflat)))
                    delete(p_pkflat(logical(p_pkflat) & ishandle(p_pkflat)))
                end
                if (any(p_pkflatmark) && any(ishandle(p_pkflatmark)))
                    delete(p_pkflatmark(logical(p_pkflatmark) & ishandle(p_pkflatmark)))
                end
                if (any(p_pk) && any(ishandle(p_pk)))
                    delete(p_pk(logical(p_pk) & ishandle(p_pk)))
                end
                if (any(p_pksmooth) && any(ishandle(p_pksmooth)))
                    delete(p_pksmooth(logical(p_pksmooth) & ishandle(p_pksmooth)))
                end
                if (any(p_pksmoothflat) && any(ishandle(p_pksmoothflat)))
                    delete(p_pksmoothflat(logical(p_pksmoothflat) & ishandle(p_pksmoothflat)))
                end
                
                [p_flat, p_pk, p_pkflat, p_pkflatmark, p_pksmooth, p_pksmoothflat] ...
                            = deal(zeros(1, pk.num_layer));
                
                mean_done   = false;
                mean_flat
                pause(time_pause)
                
                % flatten ind_y_flat_mean
                warning('off', 'MATLAB:interp1:NaNinY')
                tmp1        = NaN(pk.num_layer, num_mean);
                for ii = 1:num_mean
                    [~, tmp4] ...
                            = unique(ind_y_flat(:, pk.ind_x_mean(ii)));
                    tmp1(~isnan(tmp2(:, pk.ind_x_mean(ii))), ii) ...
                            = interp1(ind_y_flat(tmp4, pk.ind_x_mean(ii)), tmp4, tmp2(~isnan(tmp2(:, pk.ind_x_mean(ii))), pk.ind_x_mean(ii)), 'nearest', 'extrap');
                end
                for ii = 1:pk.num_layer
                    pk.layer(ii).ind_y_flat_mean ...
                            = tmp1(ii, :); % get back to the pk structure
                end
                
                tmp5        = [];
                % re-flatten ind_y_flat_mean and redo ind_y
                for ii = 1:pk.num_layer
                    
                    tmp1    = pk.layer(ii).ind_y_flat_mean;
                    % search again around min/max
                    for jj = find(~isnan(pk.layer(ii).ind_y_flat_mean))
                        try
                            [~, pk.layer(ii).ind_y_flat_mean(jj)] ...
                            = eval([pk.layer(ii).type '(amp_flat_mean((pk.layer(ii).ind_y_flat_mean(jj) - pk.num_win):(pk.layer(ii).ind_y_flat_mean(jj) + pk.num_win), jj));']);
                        catch
                            continue % have to do this because range may not be valid; not just a NaN issue
                        end
                    end
                    pk.layer(ii).ind_y_flat_mean ...
                            = pk.layer(ii).ind_y_flat_mean - (pk.num_win + 1) + tmp1; % adjust because of narrower window
                    if all(isnan(pk.layer(ii).ind_y_flat_mean))
                        tmp5= [tmp5 ii]; %#ok<AGROW>
                        continue
                    end
                    
                    p_flat(ii) = plot(block.dist_lin([1 block.num_trace]), (1e6 .* block.twtt(repmat(round(ind_y_flat(round(ind_y_pk(ii)), ind_x_pk)), 1, 2))), 'w:', 'linewidth', 2);
                    % show starter picks again (not really necessary but nice to be consistent)
                    if ~isnan(pk.layer(ii).ind_y_flat_mean(tmp3))
                        p_pkflatmark(ii) ...
                            = plot(block.dist_lin(pk.ind_x_mean(tmp3)), (1e6 * block.twtt(pk.layer(ii).ind_y_flat_mean(tmp3))), 'ko', 'markersize', 12, 'markerfacecolor', 'y', 'visisble', 'off');
                    else
                        p_pkflatmark(ii) = plot(block.dist_lin(pk.ind_x_mean(find(~isnan(pk.layer(ii).ind_y_flat_mean), 1, 'first'))), ...
                                               (1e6 * block.twtt(pk.layer(ii).ind_y_flat_mean(find(~isnan(pk.layer(ii).ind_y_flat_mean), 1, 'first')))), 'ko', 'markersize', 12, 'markerfacecolor', 'y', ...
                                               'visible', 'off');
                    end
                    p_pkflat(ii) = plot(block.dist_lin(pk.ind_x_mean(~isnan(pk.layer(ii).ind_y_flat_mean))), (1e6 .* block.twtt(pk.layer(ii).ind_y_flat_mean(~isnan(pk.layer(ii).ind_y_flat_mean)))), ...
                                       'r.', 'markersize', 12);
                    
                    % now redo ind_y in the same way as pk_flat
                    pk.layer(ii).ind_y  ...
                            = NaN(1, block.num_trace);
                    tmp1    = round(interp1(pk.ind_x_mean, pk.layer(ii).ind_y_flat_mean, 1:block.num_trace, 'linear', 'extrap'));
                    tmp1((tmp1 < 1) | (tmp1 > num_sample_trim)) ...
                            = NaN;
                    tmp2    = 1:block.num_trace;
                    tmp2    = tmp2(~isnan(tmp1));
                    pk.layer(ii).ind_y(~isnan(tmp1)) ...
                            = round(ind_y_flat(sub2ind([num_sample_trim block.num_trace], tmp1(~isnan(tmp1)), tmp2)));
                    pk.layer(ii).ind_y((pk.layer(ii).ind_y < 1) | (pk.layer(ii).ind_y > num_sample_trim)) ...
                            = NaN;
                    for jj = 1:(num_mean - 1)
                        if all(isnan(pk.layer(ii).ind_y_flat_mean(jj:(jj + 1))))
                            pk.layer(ii).ind_y(pk.ind_x_mean(jj):pk.ind_x_mean(jj + 1)) = NaN; % deal with breaks in layers
                        end
                    end
                    
                    % re-do best layer plots
                    p_pk(ii) = plot(block.dist_lin(ind_decim(~isnan(pk.layer(ii).ind_y(ind_decim)))), (1e6 .* block.twtt(round(pk.layer(ii).ind_y(ind_decim(~isnan(pk.layer(ii).ind_y(ind_decim))))))), ...
                                    'r.', 'markersize', 12, 'visible', 'off');
                    
                    pause(time_pause)
                    
                end
                warning('on', 'MATLAB:interp1:NaNinY')
                
                if ~isempty(tmp5) % get rid of screwy layers
                    [p_flat, p_pkflatmark, p_pk, p_pkflat, pk.layer, smooth_done] ...
                            = deal(p_flat(setdiff(1:pk.num_layer, tmp5)), p_pk(setdiff(1:pk.num_layer, tmp5)), p_pkflat(setdiff(1:pk.num_layer, tmp5)), p_pkflatmark(setdiff(1:pk.num_layer, tmp5)), ...
                                   pk.layer(setdiff(1:pk.num_layer, tmp5)), smooth_done(setdiff(1:pk.num_layer, tmp5)));
                    pk.num_layer ...
                            = pk.num_layer - length(tmp5);
                end
                
                pk_smooth
                
                set(layer_list, 'string', num2cell(1:pk.num_layer), 'value', 1)
                choose_layer
                set(pk_check, 'value', 1)
                set(smooth_check, 'value', 0)
                pk_done   = true;
                set(status_box, 'string', 'Layers re-flattened...')
                
        end
        
        % (re)flatten reference layers
        if ref_done
            if (any(p_refflat) && any(ishandle(p_refflat)))
                delete(p_refflat(logical(p_refflat) & ishandle(p_refflat)))
            end
            p_refflat       = zeros(1, pk_ref.num_layer);
            switch ref_start_end
                case 'start'
                    for ii = 1:pk_ref.num_layer
                        try
                            p_refflat(ii) = plot(block.dist_lin([1 block.ind_overlap(1)]), ...
                                            (1e6 .* block.twtt(repmat(round(ind_y_flat(round(pk_ref.layer(ii).ind_y_smooth(find(~isnan(pk_ref.layer(ii).ind_y_smooth), 1, 'last')) - pk.ind_trim_start + 1), ...
                                            (pk_ref.num_trace - pk_ref.ind_overlap(2) + 1))), 1, 2))), 'y:', 'linewidth', 2, 'visible', 'off');
                        catch
                            continue
                        end
                    end
                case 'end'
                    for ii = 1:pk_ref.num_layer
                        try
                            p_refflat(ii) = plot(block.dist_lin([block.ind_overlap(2) end]), ...
                                            (1e6 .* block.twtt(repmat(round(ind_y_flat(round(pk_ref.layer(ii).ind_y_smooth(find(~isnan(pk_ref.layer(ii).ind_y_smooth), 1, 'first')) - pk.ind_trim_start + 1), ...
                                            (pk_ref.num_trace - pk_ref.ind_overlap(2) + 1))), 1, 2))), 'y:', 'linewidth', 2, 'visible', 'off');
                        catch
                            continue
                        end
                    end
            end
        end
        
        set(flat_check, 'value', 1)
        disp_type           = 'flat';
        plot_flat
        set(status_box, 'string', 'Full flattened radargram and selected layers.')
        
    end

%% Horizontally average flattened amplitudes

    function mean_flat(source, eventdata)
        if ~flat_done
            set(status_box, 'string', 'No flattened data to average yet.')
            return
        end
        if mean_done
            ind_x_mean_old  = pk.ind_x_mean;
        end
        
        set(status_box, 'string', 'Horizontally averaging...')
        pause(time_pause)
        
        if (pk.num_ind_mean > 1)
            pk.ind_x_mean   = ceil((pk.num_ind_mean / 2) + 1):pk.num_ind_mean:(block.num_trace - ceil(pk.num_ind_mean / 2));
            num_mean        = length(pk.ind_x_mean);
            amp_flat_mean   = NaN(num_sample_trim, num_mean, 'single');
            tmp1            = floor(pk.num_ind_mean / 2);
            for ii = 1:num_mean
                amp_flat_mean(:, ii) ...
                            = nanmean(amp_flat(:, (pk.ind_x_mean(ii) - tmp1):(pk.ind_x_mean(ii) + tmp1)), 2);
            end
        elseif (pk.num_ind_mean == 1)
            pk.ind_x_mean   = 1:block.num_trace;
            num_mean        = block.num_trace;
            amp_flat_mean   = amp_flat;
        end
        
        flat_switch         = 'mean';
        set(disp_group, 'selectedobject', disp_check(5))
        disp_type           = 'flat';
        plot_flat
        set(mean_check, 'value', 1)
        
        if mean_done % redo averaging
            if (any(p_pkflat) && any(ishandle(p_pkflat)))
                delete(p_pkflat);
            end
            warning('off', 'MATLAB:interp1:NaNinY')
            for ii = 1:pk.num_layer
                pk.layer(ii).ind_y_flat_mean ...
                            = round(interp1(ind_x_mean_old, pk.layer(ii).ind_y_flat_mean, pk.ind_x_mean, 'linear', 'extrap'));
                p_pkflat(ii)= plot(block.dist_lin(pk.ind_x_mean(~isnan(pk.layer(ii).ind_y_flat_mean))), (1e6 .* block.twtt(pk.layer(ii).ind_y_flat_mean(~isnan(pk.layer(ii).ind_y_flat_mean)))), ...
                                   'r.', 'markersize', 12, 'visible', 'off');
            end
            warning('on', 'MATLAB:interp1:NaNinY')
        end
        
        mean_done           = true;
        show_pk
        set(status_box, 'string', 'Horizontally averaged, flattened radargram.')
    end

%% Pick best layers in flattened radargram

    function pk_flat(source, eventdata)
        
        % don't proceed if no horizontal averaging done
        if ~mean_done
            set(status_box, 'string', 'No attempt at horizontally averaging yet.')
            return
        elseif ~strcmp(flat_switch, 'mean') % check if mean has been done, if so display
            set(mean_check, 'value', 1)
            plot_flat
        end
        
        axes(ax_radar)
        if ~strcmp(disp_type, 'flat')
            set(disp_group, 'selectedobject', disp_check(5))
            disp_type       = 'flat';
            if pk_done
                set(mean_check, 'value', 1)
                set(pk_check, 'value', 1)
                show_mean
            end
        end
        
        pk_done             = false;
        if ~pk.num_layer % initialize for no layers
            smooth_done     = logical([]);
            [p_pksmooth, p_pksmoothflat] ...
                            = deal([]);
        end
        
        % layer picking callbacks
        tmp3                = pk.num_layer + 1; % used to shorten real/flatten conversion loop after this while picking loop
        set(pkgui, 'keypressfcn', [])
        while true
            
            set(status_box, 'string', 'H: Pick high; D: delete; U: undo; L: cut left; R: cut right; C: cut chunk; M: merge; Q: done.')
            
            % get pick and convert to indices
            [ind_x_pk, ind_y_pk, button] ...
                            = ginput(1);
            [ind_x_pk, ind_y_pk] ...
                            = deal(interp1(block.dist_lin(pk.ind_x_mean), 1:num_mean, ind_x_pk, 'nearest', 'extrap'), interp1(block.twtt, 1:num_sample_trim, (1e-6 * ind_y_pk), 'nearest', 'extrap'));
            
            if strcmpi(char(button), 'H') % trace peak, same approach as trough except using max instead of min
                
                curr_type   = 'max';
                pk.num_layer= pk.num_layer + 1;
                curr_layer  = pk.num_layer;
                [smooth_done, p_pksmooth, p_pksmoothflat] ...
                            = deal([smooth_done false], [p_pksmooth 0], [p_pksmoothflat 0]);
                prop_flat
                set(status_box, 'string', ['Layer #' num2str(curr_layer) ' picked.'])
                
            elseif strcmpi(char(button), 'D') % delete layer
                
                if (pk.num_layer > 0) % check if there are any layers
                    tmp1    = NaN(pk.num_layer, 1);
                    for ii = 1:pk.num_layer
                        tmp1(ii) ...
                            = pk.layer(ii).ind_y_flat_mean(ind_x_pk); % y index at x index pick for each layer
                    end
                    try
                        tmp1 = interp1(tmp1(~isnan(tmp1)), find(~isnan(tmp1)), ind_y_pk, 'nearest', 'extrap');
                    catch
                        set(status_box, 'string', 'Cannot determine which layer to delete. Pick a more distinct x index.')
                        pause(time_pause)
                        continue
                    end
                    delete([p_pkflatmark(tmp1) p_pkflat(tmp1)])
                    if smooth_done(tmp1)
                        delete([p_pksmooth(tmp1) p_pksmoothflat(tmp1)])
                    end
                    tmp2    = setdiff(1:pk.num_layer, tmp1);
                    [p_pkflatmark, p_pkflat, pk.layer, smooth_done, p_pksmooth, p_pksmoothflat, pk.num_layer] ...
                            = deal(p_pkflatmark(tmp2), p_pkflat(tmp2), pk.layer(tmp2), smooth_done(tmp2), p_pksmooth(tmp2), p_pksmoothflat(tmp2), (pk.num_layer - 1));
                    set(status_box, 'string', ['Deleted layer #' num2str(tmp1) '.'])
                    if (tmp1 < tmp3)
                        tmp3= tmp3 - 1; % remove 1 from the set of old layers
                    end
                end
                
            elseif strcmpi(char(button), 'U') % undo
                
                if (pk.num_layer >= tmp3)
                    delete([p_pkflatmark(end) p_pkflat(end)]);
                    [p_pkflatmark, p_pkflat, pk.layer, smooth_done, p_pksmooth, p_pksmoothflat, pk.num_layer] ...
                            = deal(p_pkflatmark(1:(end - 1)), p_pkflat(1:(end - 1)), pk.layer(1:(end - 1)), smooth_done(1:(end - 1)), p_pksmooth(1:(end - 1)), p_pksmoothflat(1:(end - 1)), (pk.num_layer - 1));
                    set(status_box, 'string', 'Undid last layer.')
                end
                
            elseif strcmpi(char(button), 'L') % delete layer to left of pick
                
                if (pk.num_layer >= tmp3)
                    tmp1    = NaN((pk.num_layer - tmp3 + 1), 1);
                    try
                        for ii = 1:(pk.num_layer - tmp3 + 1);
                            tmp1(ii)...
                            = pk.layer(tmp3 + ii - 1).ind_y_flat_mean(ind_x_pk); % y index at x index pick for each layer
                        end
                    catch
                        set(status_box, 'string', 'Something went wrong compiling new layer y-indices. Try again.')
                        pause(time_pause)
                        continue
                    end
                    tmp2    = find(~isnan(tmp1));
                    if ((length(tmp2) > 1) && (length(unique(tmp1(tmp2))) == length(tmp2)))
                        tmp1= interp1(tmp1(tmp2), tmp2, ind_y_pk, 'nearest', 'extrap') + tmp3 - 1;
                    elseif (length(tmp2) == 1)
                        tmp1= tmp3 - 1 + tmp2;
                    else
                        set(status_box, 'string', 'Cannot determine which new layer to edit. Pick a more distinct x index.')
                        continue
                    end
                    pk.layer(tmp1).ind_y_flat_mean(1:ind_x_pk) ...
                            = NaN;
                    delete(p_pkflat(tmp1))
                    if all(isnan(pk.layer(tmp1).ind_y_flat_mean))
                        delete(p_pkflatmark(tmp1));
                        tmp4= setdiff(1:pk.num_layer, tmp1);
                        [p_pkflatmark, p_pkflat, pk.layer, smooth_done, p_pksmooth, p_pksmoothflat, pk.num_layer] ...
                            = deal(p_pkflatmark(tmp4), p_pkflat(tmp4), pk.layer(tmp4), smooth_done(tmp4), p_pksmooth(tmp4), p_pksmoothflat(tmp4), (pk.num_layer - 1));
                        set(status_box, 'string', 'Deleted edited layer because it is now empty.')
                    else
                        p_pkflat(tmp1) ...
                            = plot(block.dist_lin(pk.ind_x_mean(~isnan(pk.layer(tmp1).ind_y_flat_mean))), (1e6 .* block.twtt(pk.layer(tmp1).ind_y_flat_mean(~isnan(pk.layer(tmp1).ind_y_flat_mean)))), ...
                                   '.', 'color', [1 0.7 0.7], 'markersize', 12);
                    end
                    set(status_box, 'string', ['Cut layer #' num2str(tmp1) ' left at ' num2str(ind_x_pk, '%3.1f') ' km.'])
                end
                
            elseif strcmpi(char(button), 'R') % delete layer to right of pick
                
                if (pk.num_layer >= tmp3)
                    tmp1    = NaN((pk.num_layer - tmp3 + 1), 1);
                    try
                        for ii = 1:(pk.num_layer - tmp3 + 1);
                            tmp1(ii)...
                            = pk.layer(tmp3 + ii - 1).ind_y_flat_mean(ind_x_pk); % y index at x index pick for each layer
                        end
                    catch
                        set(status_box, 'string', 'Something went wrong compiling new layer y-indices. Try again.')
                        pause(time_pause)
                        continue
                    end
                    tmp2    = find(~isnan(tmp1));
                    if ((length(tmp2) > 1) && (length(unique(tmp1(tmp2))) == length(tmp2)))
                        tmp1= interp1(tmp1(tmp2), tmp2, ind_y_pk, 'nearest', 'extrap') + tmp3 - 1;
                    elseif (length(tmp2) == 1)
                        tmp1= tmp3 - 1 + tmp2;
                    else
                        set(status_box, 'string', 'Cannot determine which new layer to edit. Pick a more distinct x index.')
                        pause(time_pause)
                        continue
                    end
                    pk.layer(tmp1).ind_y_flat_mean(ind_x_pk:end) ...
                            = NaN;
                    delete(p_pkflat(tmp1));
                    if all(isnan(pk.layer(tmp1).ind_y_flat_mean))
                        delete(p_pkflatmark(tmp1));
                        tmp4= setdiff(1:pk.num_layer, tmp1);
                        [p_pkflatmark, p_pkflat, pk.layer, smooth_done, p_pksmooth, p_pksmoothflat, pk.num_layer] ...
                            = deal(p_pkflatmark(tmp4), p_pkflat(tmp4), pk.layer(tmp4), smooth_done(tmp4), p_pksmooth(tmp4), p_pksmoothflat(tmp4), (pk.num_layer - 1));
                        set(status_box, 'string', 'Deleted edited layer because it is now empty.')
                    else
                        p_pkflat(tmp1) ...
                            = plot(block.dist_lin(pk.ind_x_mean(~isnan(pk.layer(tmp1).ind_y_flat_mean))), (1e6 .* block.twtt(pk.layer(tmp1).ind_y_flat_mean(~isnan(pk.layer(tmp1).ind_y_flat_mean)))), ...
                                   '.', 'color', [1 0.7 0.7], 'markersize', 12);
                        set(status_box, 'string', ['Cut layer #' num2str(tmp1) ' right at ' num2str(ind_x_pk, '%3.1f') ' km.'])
                    end
                end
                
            elseif strcmpi(char(button), 'C') % cut out a portion of layer
                
                if (pk.num_layer >= tmp3)
                    tmp1    = NaN((pk.num_layer - tmp3 + 1), 1);
                    try
                        for ii = 1:(pk.num_layer - tmp3 + 1)
                            tmp1(ii)...
                            = pk.layer(tmp3 + ii - 1).ind_y_flat_mean(ind_x_pk); % y index at x index pick for each layer
                        end
                    catch
                        set(status_box, 'string', 'Something went wrong compiling new layer y-indices. Try again.')
                        pause(time_pause)
                        continue
                    end
                    tmp2    = find(~isnan(tmp1));
                    if ((length(tmp2) > 1) && (length(unique(tmp1(tmp2))) == length(tmp2)))
                        tmp1= interp1(tmp1(tmp2), tmp2, ind_y_pk, 'nearest', 'extrap') + tmp3 - 1;
                    elseif (length(tmp2) == 1)
                        tmp1= tmp3 - 1 + tmp2;
                    else
                        set(status_box, 'string', 'Cannot determine which new layer to edit. Pick a more distinct x index.')
                        pause(time_pause)
                        continue
                    end
                    set(status_box, 'string', 'Now choose right end of cut...');
                    [tmp2, ~]   = ginput(1);
                    tmp2        = interp1(block.dist_lin(pk.ind_x_mean), 1:num_mean, tmp2, 'nearest', 'extrap');
                    pk.layer(tmp1).ind_y_flat_mean(ind_x_pk:tmp2) ...
                            = NaN;
                    delete(p_pkflat(tmp1))
                    if all(isnan(pk.layer(tmp1).ind_y_flat_mean))
                        delete(p_pkflatmark(tmp1))
                        tmp5= setdiff(1:pk.num_layer, tmp1);
                        [p_pkflatmark, p_pkflat, pk.layer, smooth_done, p_pksmooth, p_pksmoothflat, pk.num_layer] ...
                            = deal(p_pkflatmark(tmp5), p_pkflat(tmp5), pk.layer(tmp5), smooth_done(tmp5), p_pksmooth(tmp5), p_pksmoothflat(tmp5), (pk.num_layer - 1));
                        set(status_box, 'string', 'Deleted edited layer because it is now empty.')
                    else
                        p_pkflat(tmp1) ...
                            = plot(block.dist_lin(pk.ind_x_mean(~isnan(pk.layer(tmp1).ind_y_flat_mean))), (1e6 .* block.twtt(pk.layer(tmp1).ind_y_flat_mean(~isnan(pk.layer(tmp1).ind_y_flat_mean)))), ...
                                   '.', 'color', [1 0.7 0.7], 'markersize', 12);
                        set(status_box, 'string', ['Cut layer #' num2str(tmp1) ' between ' num2str(ind_x_pk, '%3.1f') ' and ' num2str(tmp2, '%3.1f') ' km.'])
                    end
                end
                
            elseif strcmpi(char(button), 'M') % merge two new layers
                
                if (pk.num_layer >= (tmp3 + 1))
                    
                    tmp1    = NaN((pk.num_layer - tmp3 + 1), 1);
                    try
                        for ii = 1:(pk.num_layer - tmp3 + 1);
                            tmp1(ii)...
                            = pk.layer(tmp3 + ii - 1).ind_y_flat_mean(ind_x_pk); % y index at x index pick for each layer
                        end
                    catch
                        set(status_box, 'string', 'Something went wrong compiling new layer y-indices. Try again.')
                        pause(time_pause)
                        continue
                    end
                    tmp2    = find(~isnan(tmp1));
                    if ((length(tmp2) > 1) && (length(unique(tmp1(tmp2))) == length(tmp2)))
                        tmp1= interp1(tmp1(tmp2), tmp2, ind_y_pk, 'nearest', 'extrap') + tmp3 - 1;
                    elseif (length(tmp2) == 1)
                        tmp1= tmp3 - 1 + tmp2;
                    else
                        set(status_box, 'string', 'Cannot determine which new layer to edit. Pick a more distinct x index.')
                        pause(time_pause)
                        continue
                    end
                    
                    set(status_box, 'string', 'Now choose layer to merge with...')
                    [ind_x_pk, ind_y_pk] ...
                            = ginput(1);
                    [ind_x_pk, ind_y_pk] ...
                            = deal(interp1(block.dist_lin(pk.ind_x_mean), 1:num_mean, ind_x_pk, 'nearest', 'extrap'), interp1(block.twtt, 1:num_sample_trim, (1e-6 * ind_y_pk), 'nearest', 'extrap'));
                    tmp2    = NaN((pk.num_layer - tmp3 + 1), 1);
                    try
                        for ii = 1:(pk.num_layer - tmp3 + 1);
                            tmp2(ii)...
                            = pk.layer(tmp3 + ii - 1).ind_y_flat_mean(ind_x_pk); % y index at x index pick for each layer
                        end
                    catch
                        set(status_box, 'string', 'Something went wrong compiling new layer y-indices. Try again.')
                        pause(time_pause)
                        continue
                    end
                    tmp4    = find(~isnan(tmp2));
                    if ((length(tmp4) > 1) && (length(unique(tmp2(tmp4))) == length(tmp4)))
                        tmp2= interp1(tmp2(tmp4), tmp4, ind_y_pk, 'nearest', 'extrap') + tmp3 - 1;
                    elseif (length(tmp4) == 1)
                        tmp2= tmp3 - 1 + tmp4;
                    else
                        set(status_box, 'string', 'Cannot determine which layer to merge with. Pick a more distinct x index.')
                        pause(time_pause)
                        continue
                    end
                    
                    pk.layer(tmp1).ind_y_flat_mean(isnan(pk.layer(tmp1).ind_y_flat_mean)) ...
                            = pk.layer(tmp2).ind_y_flat_mean(isnan(pk.layer(tmp1).ind_y_flat_mean));
                    delete(p_pkflat([tmp1 tmp2]))
                    delete(p_pkflatmark(tmp2))
                    tmp4    = setdiff(1:pk.num_layer, tmp2);
                    p_pkflat(tmp1) ...
                            = plot(block.dist_lin(pk.ind_x_mean(~isnan(pk.layer(tmp1).ind_y_flat_mean))), (1e6 .* block.twtt(pk.layer(tmp1).ind_y_flat_mean(~isnan(pk.layer(tmp1).ind_y_flat_mean)))), ...
                                   'c.', 'markersize', 12);
                    pause(time_pause)
                    set(p_pkflat(tmp1), 'color', [1 0.7 0.7])
                    [p_pkflatmark, p_pkflat, pk.layer, smooth_done, p_pksmooth, p_pksmoothflat, pk.num_layer] ...
                            = deal(p_pkflatmark(tmp4), p_pkflat(tmp4), pk.layer(tmp4), smooth_done(tmp4), p_pksmooth(tmp4), p_pksmoothflat(tmp4), (pk.num_layer - 1));
                    set(status_box, 'string', ['Layers #' num2str(tmp1) ' and ' num2str(tmp2) 'merged.'])
                    
                else
                    set(status_box, 'string', 'Not enough new layers to merge.')
                    continue
                end
                
            elseif strcmpi(char(button), 'Q') % done picking lines
                
                set(status_box, 'string', 'Done picking flattened layers...')
                break
                
            end
        end
        set(pkgui, 'keypressfcn', @keypress)
        
        if ~pk.num_layer
            set(status_box, 'string', 'No layers picked/left.')
            return
        end
        set(p_pkflat, 'color', 'r')
        
        % convert flat indices back into real twtt indices
        warning('off', 'MATLAB:interp1:NaNinY')
        for ii = 1:pk.num_layer
            pk.layer(ii).ind_y  = NaN(1, block.num_trace);
            tmp1                = round(interp1(pk.ind_x_mean, pk.layer(ii).ind_y_flat_mean, 1:block.num_trace, 'linear', 'extrap'));
            tmp1((tmp1 < 1) | (tmp1 > num_sample_trim)) ...
                                = NaN;
            tmp2                = 1:block.num_trace;
            tmp2                = tmp2(~isnan(tmp1));
            pk.layer(ii).ind_y(~isnan(tmp1)) ...
                                = ind_y_flat(sub2ind([num_sample_trim block.num_trace], tmp1(~isnan(tmp1)), tmp2));
            pk.layer(ii).ind_y((pk.layer(ii).ind_y < 1) | (pk.layer(ii).ind_y > num_sample_trim)) ...
                                = NaN;
            for jj = 1:(num_mean - 1)
                if all(isnan(pk.layer(ii).ind_y_flat_mean(jj:(jj + 1))))
                    pk.layer(ii).ind_y(pk.ind_x_mean(jj):pk.ind_x_mean(jj + 1)) ...
                                = NaN; % deal with breaks in layers
                end
            end
        end
        warning('on', 'MATLAB:interp1:NaNinY')
        
        % plot best layers in real twtt space
        if (any(p_pk) && any(ishandle(p_pk)))
            delete(p_pk(logical(p_pk) & ishandle(p_pk)))
        end
        set(pk_check, 'value', 1)
        p_pk                = zeros(1, pk.num_layer);
        for ii = 1:pk.num_layer
            tmp1            = pk.layer(ii).ind_y(ind_decim);
            tmp2            = block.dist_lin(ind_decim);
            p_pk(ii)        = plot(tmp2(~isnan(tmp1)), (1e6 .* block.twtt(round(tmp1(~isnan(tmp1))))), 'r.', 'markersize', 12, 'visible', 'off');
        end
        
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
        
        set(layer_list, 'string', num2cell(1:pk.num_layer), 'value', curr_layer);
        
        % fix unsolved bug, which may not exist anymore
        [p_pk, p_pkflat, p_pkflatmark, p_pksmooth, p_pksmoothflat] ...
                            = deal(p_pk(1:pk.num_layer), p_pkflat(1:pk.num_layer), p_pkflatmark(1:pk.num_layer), p_pksmooth(1:pk.num_layer), p_pksmoothflat(1:pk.num_layer));
        
        pk.keep_or_flat     = 'flat';
        pk_done             = true;
        match_done          = false;
        choose_layer
        set(status_box, 'string', 'Best layers mapped and re-ordered.')
    end

%% Propagate flattened layer from pick

    function prop_flat(source, eventdata)
        
        pk.layer(curr_layer).ind_y_flat_mean ...
                            = NaN(1, num_mean);
        [~, pk.layer(curr_layer).ind_y_flat_mean(ind_x_pk)] ...
                            = eval([curr_type '(amp_flat_mean((ind_y_pk - pk.num_win):(ind_y_pk + pk.num_win), ind_x_pk));']); % y index of nearest min/max
        pk.layer(curr_layer).ind_y_flat_mean(ind_x_pk) ...
                            = ind_y_pk - ((pk.num_win + 1) - pk.layer(curr_layer).ind_y_flat_mean(ind_x_pk)); % correct y index of min/max because search was done in a narrow window
        p_pkflatmark(curr_layer) ...
                            = plot(block.dist_lin(pk.ind_x_mean(ind_x_pk)), (1e6 * block.twtt(pk.layer(curr_layer).ind_y_flat_mean(ind_x_pk))), 'ko', 'markersize', 12, 'markerfacecolor', 'y'); ...
                                   % nearest min/max at picked x index
        
        for ii = (ind_x_pk - 1):-1:1 % loop for left of ind_x_pk
            try
                [~, pk.layer(curr_layer).ind_y_flat_mean(ii)] ...
                            = eval([curr_type '(amp_flat_mean((pk.layer(curr_layer).ind_y_flat_mean(ii + 1) - pk.num_win):(pk.layer(curr_layer).ind_y_flat_mean(ii + 1) + pk.num_win), ii));']);
                pk.layer(curr_layer).ind_y_flat_mean(ii) ...
                            = pk.layer(curr_layer).ind_y_flat_mean(ii + 1) - ((pk.num_win + 1) - pk.layer(curr_layer).ind_y_flat_mean(ii));
            catch %#ok<*CTCH>
                break
            end
        end
        for ii = (ind_x_pk + 1):num_mean % loop for right of ind_x_pk
            try
                [~, pk.layer(curr_layer).ind_y_flat_mean(ii)] ...
                            = eval([curr_type '(amp_flat_mean((pk.layer(curr_layer).ind_y_flat_mean(ii - 1) - pk.num_win):(pk.layer(curr_layer).ind_y_flat_mean(ii - 1) + pk.num_win), ii));']);
                pk.layer(curr_layer).ind_y_flat_mean(ii) ...
                            = pk.layer(curr_layer).ind_y_flat_mean(ii - 1) - ((pk.num_win + 1) - pk.layer(curr_layer).ind_y_flat_mean(ii));
            catch %#ok<CTCH>
                break
            end
        end
        
        p_pkflat(curr_layer) = plot(block.dist_lin(pk.ind_x_mean(~isnan(pk.layer(curr_layer).ind_y_flat_mean))), ...
                                   (1e6 .* block.twtt(pk.layer(curr_layer).ind_y_flat_mean(~isnan(pk.layer(curr_layer).ind_y_flat_mean)))), 'k.', 'color', [1 0.7 0.7], 'markersize', 12);
        pk.layer(curr_layer).type ...
                            = curr_type;
    end

%% Manually pick a layer to keep

    function pk_man(source, eventdata)
        
        if ~trim_done
            set(status_box, 'string', 'Trim excess data first.')
            return
        end
        if ~flat_done
            set(status_box, 'string', 'Cannot trace a layer manually if flattening not done.')
            return
        end
        if (~strcmp(disp_type, 'amp.') || (get(disp_group, 'selectedobject') ~= disp_check(1)))
            set(disp_group, 'selectedobject', disp_check(1))
            disp_type       = 'amp.';
            plot_db
        end
        
        axes(ax_radar)
        ii                  = 0;
        [ind_x_pk, ind_y_pk, button] ...
                            = deal([]);
        
        set(pkgui, 'keypressfcn', [])
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
                ind_y_pk(ii)= interp1(block.twtt, 1:num_sample_trim, (1e-6 .* ind_y_pk(ii)), 'nearest', 'extrap'); % interpolate traveltime pick onto traveltime vector
                if (ii > 1)
                    delete(tmp3) % get rid of old plot handle
                end
                tmp3        = plot(block.dist_lin(ind_x_pk), (1e6 .* block.twtt(ind_y_pk)), 'rx', 'markersize', 12); % original picks
                ii          = ii + 1;
                
            elseif (strcmpi(char(button), 'U') && (ii > 1))
                
                ind_x_pk    = ind_x_pk(1:(end - 1));
                ind_y_pk    = ind_y_pk(1:(end - 1));
                if (logical(tmp3) && ishandle(tmp3))
                    delete(tmp3)
                end
                if (ii > 2)
                    tmp3    = plot(block.dist_lin(ind_x_pk), (1e6 .* block.twtt(ind_y_pk)), 'rx', 'markersize', 12); % original picks
                end
                ii          = ii - 1;
                
            elseif (double(get(pkgui, 'currentcharacter')) == 13)
                
                if (ii < 4)
                    if (logical(tmp3) && ishandle(tmp3))
                        delete(tmp3)
                    end
                    set(status_box, 'string', 'Not enough picked points to make a manual layer. Start over.')
                    set(pkgui, 'keypressfcn', @keypress)
                    break
                end
                
                if ~issorted(ind_x_pk) % resort picks
                    [ind_x_pk, tmp1] = sort(ind_x_pk);
                    ind_y_pk= ind_y_pk(tmp1);
                end
                if (length(unique(ind_x_pk)) < length(ind_x_pk)) % don't keep picks that are accidentally at the same horizontal index, otherwise spline will fail
                    [ind_x_pk, tmp1] = unique(ind_x_pk);
                    ind_y_pk= ind_y_pk(tmp1);
                end
                
                pk.num_layer= pk.num_layer + 1;
                
                pk.layer(pk.num_layer).ind_y ...
                            = NaN(1, block.num_trace);
                pk.layer(pk.num_layer).ind_y(ind_x_pk(1):decim:ind_x_pk(end)) ...
                            = round(spline(ind_x_pk, ind_y_pk, ind_x_pk(1):decim:ind_x_pk(end))); % interpolate spline through picks
                for ii = interp1(ind_decim, 1:length(ind_decim), ind_x_pk, 'nearest', 'extrap') % check for local maxima
                    [~, tmp2] ...
                            = max(amp_mean((pk.layer(pk.num_layer).ind_y(ind_decim(ii)) - pk.num_win):(pk.layer(pk.num_layer).ind_y(ind_decim(ii)) + pk.num_win), ii));
                    pk.layer(pk.num_layer).ind_y(ind_decim(ii)) ...
                            = pk.layer(pk.num_layer).ind_y(ind_decim(ii)) - pk.num_win - 1 + tmp2; % adjust local maximum index
                end
                pk.layer(pk.num_layer).ind_y(ind_x_pk(1):ind_x_pk(end)) ...
                            = round(spline(ind_x_pk, pk.layer(pk.num_layer).ind_y(ind_x_pk), ind_x_pk(1):ind_x_pk(end))); % interpolate spline through improved picks
                
                tmp1        = find(~isnan(pk.layer(pk.num_layer).ind_y(ind_decim))); % 
                delete(tmp3)
                p_pk(pk.num_layer) ...
                            = plot(block.dist_lin(ind_decim(tmp1)), (1e6 .* block.twtt(round(pk.layer(pk.num_layer).ind_y(ind_decim(tmp1))))), 'r.', 'markersize', 12, 'visible', 'off');
                
                warning('off', 'MATLAB:interp1:NaNinY')
                pk.layer(pk.num_layer).ind_y_flat_mean ...
                            = NaN(1, num_mean);
                for ii = interp1(pk.ind_x_mean, 1:num_mean, ind_x_pk(1), 'nearest', 'extrap'):interp1(pk.ind_x_mean, 1:num_mean, ind_x_pk(end), 'nearest', 'extrap')
                    [~, tmp1] ...
                            = unique(ind_y_flat(:, pk.ind_x_mean(ii)));
                    pk.layer(pk.num_layer).ind_y_flat_mean(ii) ...
                            = interp1(ind_y_flat(tmp1, pk.ind_x_mean(ii)), tmp1, pk.layer(pk.num_layer).ind_y(pk.ind_x_mean(ii)), 'nearest', 'extrap');
                end
                warning('on', 'MATLAB:interp1:NaNinY')
                
                p_pkflatmark(pk.num_layer) ...
                            = plot(block.dist_lin(pk.ind_x_mean(find(~isnan(pk.layer(pk.num_layer).ind_y_flat_mean), 1))), ...
                                   (1e6 * block.twtt(pk.layer(pk.num_layer).ind_y_flat_mean(find(~isnan(pk.layer(pk.num_layer).ind_y_flat_mean), 1)))), ...
                                   'ko', 'markersize', 12, 'markerfacecolor', 'y', 'visible', 'off');
                p_pkflat(pk.num_layer) ...
                            = plot(block.dist_lin(pk.ind_x_mean(~isnan(pk.layer(pk.num_layer).ind_y_flat_mean))), ...
                                   (1e6 .* block.twtt(pk.layer(pk.num_layer).ind_y_flat_mean(~isnan(pk.layer(pk.num_layer).ind_y_flat_mean)))), 'r.', 'markersize', 12, 'visible', 'off');
                
                [pk.layer(pk.num_layer).type, smooth_done, p_pksmooth, p_pksmoothflat] ...
                            = deal('max', [smooth_done false], [p_pksmooth 0], [p_pksmoothflat 0]);
                
                pk_smooth
                pk_sort
                curr_layer  = find(tmp1 == pk.num_layer); % tmp1 from pk_sort
                set(layer_list, 'string', num2cell(1:pk.num_layer), 'value', curr_layer);
                choose_layer
                set(pk_check, 'value', 1)
                show_pk
                set(status_box, 'string', 'Manual layer successfully picked.')
                set(pkgui, 'keypressfcn', @keypress)
                return
            end
        end
    end

%% Sort layers from top to bottom based on their average vertical index

    function pk_sort(source, eventdata)
        tmp1                = zeros(pk.num_layer, 1);
        for ii = 1:pk.num_layer
            tmp1(ii)        = nanmean(pk.layer(ii).ind_y);
        end
        [~, tmp1]           = sort(tmp1);
        [pk.layer, p_pk, p_pkflatmark, p_pkflat, smooth_done, p_pksmooth, p_pksmoothflat] ...
                            = deal(pk.layer(tmp1), p_pk(tmp1), p_pkflatmark(tmp1), p_pkflat(tmp1), smooth_done(tmp1), p_pksmooth(tmp1), p_pksmoothflat(tmp1));
    end

%% Choose/highlight the current layer

    function choose_layer(source, eventdata)
        curr_layer          = get(layer_list, 'value');
        if pk.num_layer
            if (any(p_pkflat) && any(ishandle(p_pkflat)))
                set(p_pkflat(logical(p_pkflat) & ishandle(p_pkflat)), 'markersize', 12)
            end
            if (any(p_pk) && any(ishandle(p_pk)))
                set(p_pk(logical(p_pk) & ishandle(p_pk)), 'markersize', 12)
            end
            if (any(p_pksmooth) && any(ishandle(p_pksmooth)))
                set(p_pksmooth(logical(p_pksmooth) & ishandle(p_pksmooth)), 'markersize', 12)
            end
            if (any(p_pksmoothflat) && any(ishandle(p_pksmoothflat)))
                set(p_pksmoothflat(logical(p_pksmoothflat) & ishandle(p_pksmoothflat)), 'markersize', 12)
            end
            if (logical(p_pkflat(curr_layer)) && ishandle(p_pkflat(curr_layer)))
                set(p_pkflat(curr_layer), 'markersize', 24)
            end
            if (logical(p_pk(curr_layer)) && ishandle(p_pk(curr_layer)))
                set(p_pk(curr_layer), 'markersize', 24)
            end
            if smooth_done(curr_layer)
                if (logical(p_pksmooth(curr_layer)) && ishandle(p_pksmooth(curr_layer)))
                    set(p_pksmooth(curr_layer), 'markersize', 24)
                end
                if (logical(p_pksmoothflat(curr_layer)) && ishandle(p_pksmoothflat(curr_layer)))
                    set(p_pksmoothflat(curr_layer), 'markersize', 24)
                end
            end
        end
    end

%% Choose/select a layer interactively

    function choose_pk(source, eventdata)
        if ~pk_done
            set(status_box, 'string', 'No best layers to focus on.')
            return
        end
        if (strcmp(disp_type, 'flat') && ~strcmp(flat_switch, 'mean'))
            flat_switch     = 'mean';
            plot_flat
        end
        set(status_box, 'string', 'Choose a layer to highlight...')
        [ind_x_pk, ind_y_pk]= ginput(1);
        switch disp_type
            case 'amp.'
                ind_x_pk    = interp1(block.dist_lin(ind_decim), ind_decim, ind_x_pk, 'nearest', 'extrap');
            case 'flat'
                ind_x_pk    = pk.ind_x_mean(interp1(block.dist_lin(pk.ind_x_mean), 1:num_mean, ind_x_pk, 'nearest', 'extrap'));
        end
        ind_y_pk            = interp1(block.twtt, 1:num_sample_trim, (1e-6 * ind_y_pk), 'nearest', 'extrap');
        tmp1                = NaN(pk.num_layer, 1);
        for ii = 1:pk.num_layer
            switch disp_type
                case 'amp.'
                    tmp1(ii)= pk.layer(ii).ind_y(ind_x_pk);
                case 'flat'
                    tmp1(ii)= pk.layer(ii).ind_y_flat_mean(interp1(pk.ind_x_mean, 1:num_mean, ind_x_pk, 'nearest', 'extrap'));
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
        if (isempty(curr_layer) || any(curr_layer < 1) || any(curr_layer > pk.num_layer));
            curr_layer      = 1;
        end
        set(layer_list, 'value', curr_layer)
        choose_layer
        set(status_box, 'string', ['Layer #' num2str(curr_layer) ' chosen.'])
    end

%% Focus on a layer vertically

    function focus_layer(source, eventdata)
        if ~pk_done
            set(status_box, 'string', 'No best layers to focus on.')
            return
        end
        axes(ax_radar)
        switch disp_type
            case 'amp.'
                ylim(1e6 .* block.twtt(round([min(pk.layer(curr_layer).ind_y) max(pk.layer(curr_layer).ind_y)])))
            case 'flat'
                ylim(1e6 .* block.twtt(round([min(pk.layer(curr_layer).ind_y_flat_mean) max(pk.layer(curr_layer).ind_y_flat_mean)])))
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
        set(status_box, 'string', ['Focused on Layer #' num2str(curr_layer) '.'])
        narrow_cb
    end

%% Switch to previous layer in list

    function pk_last(source, eventdata)
        if (curr_layer > 1)
            curr_layer      = curr_layer - 1;
            set(layer_list, 'value', curr_layer)
            choose_layer
        end
    end

%% Switch to next layer in the list

    function pk_next(source, eventdata)
        if (curr_layer < pk.num_layer)
            curr_layer      = curr_layer + 1;
            set(layer_list, 'value', curr_layer)
            choose_layer
        end
    end

%% Delete layer

    function del_layer(source, eventdata)
        if (pk.num_layer && curr_layer)
            set(status_box, 'string', 'Delete current layer? Y: yes; otherwise: no.')
            waitforbuttonpress
            if ~strcmpi(get(pkgui, 'currentcharacter'), 'Y')
                set(status_box, 'string', 'Layer deletion cancelled by user.')
                return
            end
            tmp1             = setdiff(1:pk.num_layer, curr_layer);
            if (any(p_pkflatmark(curr_layer)) && any(ishandle(p_pkflatmark(curr_layer))))
                delete(p_pkflatmark(curr_layer))
            end
            if (logical(p_pkflat(curr_layer)) && ishandle(p_pkflat(curr_layer)))
                delete(p_pkflat(curr_layer))
            end
            if (logical(p_pk(curr_layer)) && ishandle(p_pk(curr_layer)))
                delete(p_pk(curr_layer))
            end
            if (logical(p_pksmooth(curr_layer)) && ishandle(p_pksmooth(curr_layer)))
                delete(p_pksmooth(curr_layer))
            end
            if (logical(p_pksmoothflat(curr_layer)) && ishandle(p_pksmoothflat(curr_layer)))
                delete(p_pksmoothflat(curr_layer))
            end
            [pk.num_layer, pk.layer, p_pkflatmark, p_pkflat, p_pk, smooth_done, p_pksmooth, p_pksmoothflat] ...
                            = deal((pk.num_layer - 1), pk.layer(tmp1), p_pkflatmark(tmp1), p_pkflat(tmp1), p_pk(tmp1), smooth_done(tmp1), p_pksmooth(tmp1), p_pksmoothflat(tmp1));
            set(status_box, 'string', ['Layer #' num2str(curr_layer) ' deleted.'])
            curr_layer      = curr_layer - 1;
            if (isempty(curr_layer) || any(curr_layer < 1) || any(curr_layer > pk.num_layer))
                curr_layer  = 1;
            end
            if pk.num_layer
                set(layer_list, 'string', num2cell(1:pk.num_layer), 'value', curr_layer)
            else
                set(layer_list, 'string', 'N/A', 'value', 1)
            end
            match_done      = false;
            choose_layer
        else
            set(status_box, 'string', 'No best layers to delete yet.')
        end
    end

%% Delete all layers

    function del_all(source, eventdata)
        set(status_box, 'string', 'Delete all layers? Y: yes; otherwise: no.')
        waitforbuttonpress
        if ~strcmpi(get(pkgui, 'currentcharacter'), 'Y')
            set(status_box, 'string', 'Layer deletion cancelled by user.')
            return
        end
        if pk.num_layer
            if (any(p_pk) && any(ishandle(p_pk)))
                delete(p_pk((logical(p_pk) & ishandle(p_pk))))
            end
            if (any(p_pkflat) && any(ishandle(p_pkflat)))
                delete(p_pkflat((logical(p_pkflat) & ishandle(p_pkflat))))
            end
            if (any(p_pkflatmark) && any(ishandle(p_pkflatmark)))
                delete(p_pkflatmark(logical(p_pkflatmark) & ishandle(p_pkflatmark)))
            end
            if (any(p_pksmooth) && any(ishandle(p_pksmooth)))
                delete(p_pksmooth(logical(p_pksmooth) & ishandle(p_pksmooth)))
            end
            if (any(p_pksmoothflat) && any(ishandle(p_pksmoothflat)))
                delete(p_pksmoothflat(logical(p_pksmoothflat) & ishandle(p_pksmoothflat)))
            end
            [pk.num_layer, p_pkflatmark, p_pkflat, p_pk, p_pksmooth, p_pksmoothflat] ...
                            = deal(0);
            pk.layer        = struct;
            smooth_done     = logical([]);
            set(layer_list, 'string', 'N/A', 'value', 1)
            set(status_box, 'string', 'All layers deleted.')
            match_done      = false;
        else
            set(status_box, 'string', 'No best layers to delete yet.')
        end
    end

%% Adjust/edit current layer

    function adj_layer(source, eventdata)
        
        if ~pk_done
            set(status_box, 'string', 'No picked layers to adjust yet.')
            return
        end
        
        axes(ax_radar)
        tmp5                = 0;
        
        set(pkgui, 'keypressfcn', [])
        while true
            
            set(status_box, 'string', 'L: delete left; R: delete right; C: cut (left start); U: undo; Q: quit.')
            
            % get pick and convert to indices
            [ind_x_pk, ~, button] ...
                            = ginput(1);
            if (strcmp(flat_switch, 'mean') && strcmp(disp_type, 'flat'));
                ind_x_pk    = pk.ind_x_mean(interp1(block.dist_lin(pk.ind_x_mean), 1:num_mean, ind_x_pk, 'nearest', 'extrap'));
            else
                ind_x_pk    = interp1(block.dist_lin(ind_decim), ind_decim, ind_x_pk, 'nearest', 'extrap');
            end
            
            if strcmpi(char(button), 'L')
                
                [tmp2, tmp3]= deal(button, cell(1, 4));
                [tmp3{1}, tmp3{2}] ...
                            = deal(pk.layer(curr_layer).ind_y(1:ind_x_pk), pk.layer(curr_layer).ind_y_flat_mean(1:interp1(pk.ind_x_mean, 1:num_mean, ind_x_pk, 'nearest', 'extrap')));
                [pk.layer(curr_layer).ind_y(1:ind_x_pk), pk.layer(curr_layer).ind_y_flat_mean(1:interp1(pk.ind_x_mean, 1:num_mean, ind_x_pk, 'nearest', 'extrap'))] ...
                            = deal(NaN);
                if smooth_done(curr_layer)
                    [tmp3{3}, tmp3{4}] ...
                            = deal(pk.layer(curr_layer).ind_y_flat_smooth(1:ind_x_pk), pk.layer(curr_layer).ind_y_smooth(1:ind_x_pk));
                    [pk.layer(curr_layer).ind_y_flat_smooth(1:ind_x_pk), pk.layer(curr_layer).ind_y_smooth(1:ind_x_pk)] ...
                            = deal(NaN);
                end
                fix_layer
                set(status_box, 'string', ['Layer #' num2str(curr_layer) ' cut left at ' num2str(block.dist_lin(ind_x_pk), '%3.1f') ' km.'])
                tmp5        = ind_x_pk;
                
            elseif strcmpi(char(button), 'R')
                
                [tmp2, tmp3]= deal(button, cell(1, 4));
                [tmp3{1}, tmp3{2}] ...
                            = deal(pk.layer(curr_layer).ind_y(ind_x_pk:end), pk.layer(curr_layer).ind_y_flat_mean(interp1(pk.ind_x_mean, 1:num_mean, ind_x_pk, 'nearest', 'extrap'):end));
                [pk.layer(curr_layer).ind_y(ind_x_pk:end), pk.layer(curr_layer).ind_y_flat_mean(interp1(pk.ind_x_mean, 1:num_mean, ind_x_pk, 'nearest', 'extrap'):end)] ...
                            = deal(NaN);
                if smooth_done(curr_layer)
                    [tmp3{3}, tmp3{4}] ...
                            = deal(pk.layer(curr_layer).ind_y_flat_smooth(ind_x_pk:end), pk.layer(curr_layer).ind_y_smooth(ind_x_pk:end));
                    [pk.layer(curr_layer).ind_y_flat_smooth(ind_x_pk:end), pk.layer(curr_layer).ind_y_smooth(ind_x_pk:end)] ...
                            = deal(NaN);
                end
                fix_layer
                set(status_box, 'string', ['Layer #' num2str(curr_layer) ' cut right at ' num2str(block.dist_lin(ind_x_pk), '%3.1f') ' km.'])
                tmp5        = ind_x_pk;
                
            elseif strcmpi(char(button), 'C')
                
                [tmp2, tmp3]= deal(button, cell(1, 4));
                set(status_box, 'string', 'Now choose right end of cut...')
                [tmp1, ~]   = ginput(1);
                if (strcmp(flat_switch, 'mean') && strcmp(disp_type, 'flat'))
                    tmp1    = pk.ind_x_mean(interp1(block.dist_lin(pk.ind_x_mean), 1:num_mean, tmp1, 'nearest', 'extrap'));
                else
                    tmp1    = interp1(block.dist_lin(ind_decim), ind_decim, tmp1, 'nearest', 'extrap');
                end
                
                [tmp3{1}, tmp3{2}] ...
                            = deal(pk.layer(curr_layer).ind_y(ind_x_pk:tmp1), ...
                                   pk.layer(curr_layer).ind_y_flat_mean(interp1(pk.ind_x_mean, 1:num_mean, ind_x_pk, 'nearest', 'extrap'):interp1(pk.ind_x_mean, 1:num_mean, tmp1, 'nearest', 'extrap')));
                [pk.layer(curr_layer).ind_y(ind_x_pk:tmp1), ...
                    pk.layer(curr_layer).ind_y_flat_mean(interp1(pk.ind_x_mean, 1:num_mean, ind_x_pk, 'nearest', 'extrap'):interp1(pk.ind_x_mean, 1:num_mean, tmp1, 'nearest', 'extrap'))] ...
                            = deal(NaN);
                if smooth_done(curr_layer)
                    [tmp3{3}, tmp3{4}] ...
                            = deal(pk.layer(curr_layer).ind_y_flat_smooth(ind_x_pk:tmp1), pk.layer(curr_layer).ind_y_smooth(ind_x_pk:tmp1));
                    [pk.layer(curr_layer).ind_y_flat_smooth(ind_x_pk:tmp1), pk.layer(curr_layer).ind_y_smooth(ind_x_pk:tmp1)] ...
                            = deal(NaN);
                end
                [tmp4, tmp5]= deal(tmp1, ind_x_pk); % in case of undo
                
                fix_layer
                set(status_box, 'string', ['Layer #' num2str(curr_layer) ' cut between ' num2str(block.dist_lin(ind_x_pk), '%3.1f') ' and ' ...
                    num2str(block.dist_lin(tmp1(find(~isnan(tmp1), 1, 'last'))), '%3.1f') ' km.'])
                
            elseif strcmpi(char(button), 'U') % undo adjustment done in current set
                
                if ~tmp5
                    continue
                end
                
                if strcmpi(char(tmp2), 'L')
                    
                    [pk.layer(curr_layer).ind_y(1:tmp5), pk.layer(curr_layer).ind_y_flat_mean(1:interp1(pk.ind_x_mean, 1:num_mean, tmp5, 'nearest', 'extrap'))] ...
                            = deal(tmp3{1}, tmp3{2});
                    if smooth_done(curr_layer)
                        [pk.layer(curr_layer).ind_y_flat_smooth(1:tmp5), pk.layer(curr_layer).ind_y_smooth(1:tmp5)] ...
                            = deal(tmp3{3}, tmp3{4});
                    end
                    
                elseif strcmpi(char(tmp2), 'R')
                    
                    [pk.layer(curr_layer).ind_y(tmp5:end), pk.layer(curr_layer).ind_y_flat_mean(interp1(pk.ind_x_mean, 1:num_mean, tmp5, 'nearest', 'extrap'):end)] ...
                            = deal(tmp3{1}, tmp3{2});
                    if smooth_done(curr_layer)
                        [pk.layer(curr_layer).ind_y_flat_smooth(tmp5:end), pk.layer(curr_layer).ind_y_smooth(tmp5:end)] ...
                            = deal(tmp3{3}, tmp3{4});
                    end
                    
                elseif strcmpi(char(tmp2), 'C')
                    
                    [pk.layer(curr_layer).ind_y(tmp5:tmp4), ...
                        pk.layer(curr_layer).ind_y_flat_mean(interp1(pk.ind_x_mean, 1:num_mean, tmp5, 'nearest', 'extrap'):interp1(pk.ind_x_mean, 1:num_mean, tmp4, 'nearest', 'extrap'))] ...
                            = deal(tmp3{1}, tmp3{2});
                    if smooth_done(curr_layer)
                        [pk.layer(curr_layer).ind_y_flat_smooth(tmp5:tmp4), pk.layer(curr_layer).ind_y_smooth(tmp5:tmp4)] ...
                            = deal(tmp3{3}, tmp3{4});
                    end
                    
                end
                
                tmp3        = 0;
                fix_layer
                set(status_box, 'string', 'Undid previous adjustment.')
                
            elseif strcmpi(char(button), 'Q')
                
                pk_sort
                tmp1        = find((tmp1 == curr_layer), 1); % tmp1 now from pk_sort
                if ((tmp1 ~= curr_layer) && ~isempty(tmp1))
                    set(status_box, 'string', ['Done adjusting Layer #' num2str(curr_layer) ' (now Layer #' num2str(tmp1) ').'])
                    curr_layer = tmp1;
                else
                    set(status_box, 'string', ['Done adjusting Layer #' num2str(curr_layer) '.'])
                end
                if (isempty(curr_layer) || any(curr_layer < 1) || any(curr_layer > pk.num_layer))
                    curr_layer = 1;
                end
                set(layer_list, 'value', curr_layer)
                choose_layer
                match_done  = false;
                set(pkgui, 'keypressfcn', @keypress)
                return
                
            end
        end
    end

%% Fix layer based on interactive adjustments

    function fix_layer(source, eventdata)
        tmp1                = ind_decim(~isnan(pk.layer(curr_layer).ind_y(ind_decim)));
        if isempty(tmp1)
            del_layer
            set(status_box, 'string', 'Layer now empty so it was deleted.')
            return
        end
        if (logical(p_pkflat(curr_layer)) && ishandle(p_pkflat(curr_layer)))
            delete(p_pkflat(curr_layer))
            p_pkflat(curr_layer) ...
                            = plot(block.dist_lin(pk.ind_x_mean(~isnan(pk.layer(curr_layer).ind_y_flat_mean))), ...
                                   (1e6 .* block.twtt(pk.layer(curr_layer).ind_y_flat_mean(~isnan(pk.layer(curr_layer).ind_y_flat_mean)))), 'r.', 'markersize', 24, 'visible', 'off');
        end
        show_flat
        if (logical(p_pk(curr_layer)) && ishandle(p_pk(curr_layer)))
            delete(p_pk(curr_layer))
            p_pk(curr_layer)= plot(block.dist_lin(tmp1), (1e6 .* block.twtt(round(pk.layer(curr_layer).ind_y(tmp1)))), 'r.', 'markersize', 24, 'visible', 'off');
        end
        show_pk
        if smooth_done(curr_layer)
            delete([p_pksmooth(curr_layer) p_pksmoothflat(curr_layer)])
            tmp1            = ind_decim(~isnan(pk.layer(curr_layer).ind_y_smooth(ind_decim)));
            p_pksmooth(curr_layer) ...
                            = plot(block.dist_lin(tmp1), (1e6 .* block.twtt(round(pk.layer(curr_layer).ind_y_smooth(tmp1)))), 'g.', 'markersize', 24, 'visible', 'off');
            tmp1            = ind_decim(~isnan(pk.layer(curr_layer).ind_y_flat_smooth(ind_decim)));
            p_pksmoothflat(curr_layer) ...
                            = plot(block.dist_lin(tmp1), (1e6 .* block.twtt(round(pk.layer(curr_layer).ind_y_flat_smooth(tmp1)))), 'g.', 'markersize', 24, 'visible', 'off');
            show_smooth
        end
        choose_layer
    end

%% Merge two layers

    function merge_layer(source, eventdata)
        
        if ~pk_done
            set(status_box, 'string', 'No best layers to merge yet.')
            return
        end
        if (pk.num_layer < 2)
            set(status_box, 'string', 'Not enough layers to merge.')
            return
        end
        if (strcmp(disp_type, 'flat') && strcmp(flat_switch, 'full'))
            flat_switch     = 'mean';
            plot_flat
        end
        
        set(status_box, 'string', ['Pick layer to merge with layer #' num2str(curr_layer) ' (Q: cancel)...'])
        
        axes(ax_radar)
        
        % get pick and convert to indices
        [ind_x_pk, ind_y_pk, button] ...
                            = ginput(1);
        if strcmpi(char(button), 'Q')
            set(status_box, 'string', 'Layer merging cancelled by user.')
            return
        end
        
        if strcmp(disp_type, 'flat')
            ind_x_pk        = pk.ind_x_mean(interp1(block.dist_lin(pk.ind_x_mean), 1:num_mean, ind_x_pk, 'nearest', 'extrap'));
        else
            ind_x_pk        = interp1(block.dist_lin(ind_decim), ind_decim, ind_x_pk, 'nearest', 'extrap');
        end
        ind_y_pk            = interp1(block.twtt, 1:num_sample_trim, (1e-6 * ind_y_pk), 'nearest', 'extrap');
        
        % get current layer positions at ind_x_pk, depending on what we're working with
        tmp1                = NaN(pk.num_layer, 1);
        for ii = 1:pk.num_layer
            switch disp_type
                case 'amp.'
                    tmp1(ii)= pk.layer(ii).ind_y(ind_x_pk);
                case 'flat'
                    tmp1(ii)= pk.layer(ii).ind_y_flat_mean(interp1(pk.ind_x_mean, 1:num_mean, ind_x_pk, 'nearest', 'extrap'));
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
        if ((smooth_done(curr_layer) && ~smooth_done(tmp1)) || (~smooth_done(curr_layer) && smooth_done(tmp1)))
            set(status_box, 'string', 'Only one of layers to be merged has been smoothed. Smoothing first...')
            tmp3            = tmp1;
            pk_smooth
            tmp1            = tmp3;
            set(status_box, 'string', 'Continuing merging...')
        end
        
        % replace NaN values in first layer with those in the second layer
        [pk.layer(curr_layer).ind_y_flat_mean(isnan(pk.layer(curr_layer).ind_y_flat_mean)), pk.layer(curr_layer).ind_y(isnan(pk.layer(curr_layer).ind_y))] ...
                            = deal(pk.layer(tmp1).ind_y_flat_mean(isnan(pk.layer(curr_layer).ind_y_flat_mean)), pk.layer(tmp1).ind_y(isnan(pk.layer(curr_layer).ind_y)));
        if all(smooth_done([curr_layer tmp1]))
            [pk.layer(curr_layer).ind_y_smooth(isnan(pk.layer(curr_layer).ind_y_smooth)), pk.layer(curr_layer).ind_y_flat_smooth(isnan(pk.layer(curr_layer).ind_y_flat_smooth))] ...
                            = deal(pk.layer(tmp1).ind_y_smooth(isnan(pk.layer(curr_layer).ind_y_smooth)), pk.layer(tmp1).ind_y_flat_smooth(isnan(pk.layer(curr_layer).ind_y_flat_smooth)));
        end
        
        % fix plots
        if (any(p_pkflatmark(tmp1)) && any(ishandle(p_pkflatmark(tmp1)))) % keep original p_pkflat, delete merged one
            delete(p_pkflatmark(tmp1))
        end
        if (any(p_pkflat([curr_layer tmp1])) && any(ishandle(p_pkflat([curr_layer tmp1]))))
            delete(p_pkflat([curr_layer tmp1]))
            p_pkflat(curr_layer) ...
                            = plot(block.dist_lin(pk.ind_x_mean(~isnan(pk.layer(curr_layer).ind_y_flat_mean))), ...
                                   (1e6 .* block.twtt(pk.layer(curr_layer).ind_y_flat_mean(~isnan(pk.layer(curr_layer).ind_y_flat_mean)))), 'r.', 'markersize', 24, 'visible', 'off');
        end
        show_flat
        tmp2                = ind_decim(~isnan(pk.layer(curr_layer).ind_y(ind_decim)));
        if (any(p_pk([curr_layer tmp1])) && any(ishandle(p_pk([curr_layer tmp1]))))
            delete(p_pk([curr_layer tmp1]))
            p_pk(curr_layer)= plot(block.dist_lin(tmp2), (1e6 .* block.twtt(round(pk.layer(curr_layer).ind_y(tmp2)))), 'r.', 'markersize', 24, 'visible', 'off');
        end
        if all(smooth_done([curr_layer tmp1]))
            delete([p_pksmooth([curr_layer tmp1]) p_pksmoothflat([curr_layer tmp1])])
            tmp2            = ind_decim(~isnan(pk.layer(curr_layer).ind_y_smooth(ind_decim)));
            p_pksmooth(curr_layer) ...
                            = plot(block.dist_lin(tmp2), (1e6 .* block.twtt(round(pk.layer(curr_layer).ind_y_smooth(tmp2)))), 'g.', 'markersize', 24, 'visible', 'off');
            tmp2            = ind_decim(~isnan(pk.layer(curr_layer).ind_y_flat_smooth(ind_decim)));
            p_pksmoothflat(curr_layer) ...
                            = plot(block.dist_lin(tmp2), (1e6 .* block.twtt(round(pk.layer(curr_layer).ind_y_flat_smooth(tmp2)))), 'g.', 'markersize', 24, 'visible', 'off');
        end
        
        tmp3                = setdiff(1:pk.num_layer, tmp1);
        [p_pkflatmark, p_pkflat, p_pk, p_pksmooth, p_pksmoothflat] ...
                            = deal(p_pkflatmark(tmp3), p_pkflat(tmp3), p_pk(tmp3), p_pksmooth(tmp3), p_pksmoothflat(tmp3));
        show_pk
        smooth_done(curr_layer) ...
                            = false;
        [pk.layer, smooth_done, pk.num_layer] ...
                            = deal(pk.layer(tmp3), smooth_done(tmp3), (pk.num_layer - 1));
        if any(smooth_done)
            show_smooth
        end
        
        set(layer_list, 'string', num2cell(1:pk.num_layer), 'value', 1)
        set(status_box, 'string', ['Layers #' num2str(curr_layer) ' and #' num2str(tmp1) ' merged.'])
        if (curr_layer > tmp1)
            curr_layer      = curr_layer - 1;
        end
        pk_sort
        curr_layer          = find((tmp1 == curr_layer), 1); % using tmp1 from pk_sort
        if (isempty(curr_layer) || any(curr_layer < 1) || any(curr_layer > pk.num_layer))
            curr_layer      = 1;
        end
        set(layer_list, 'value', curr_layer)
        match_done          = false;
        choose_layer
        
    end

%% Smooth best flattened layers

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
        
        tic
        tmp1                = [];
        warning('off', 'MATLAB:interp1:NaNinY')
        for ii = find(~smooth_done)
            pk.layer(ii).ind_y_flat_smooth ...
                            = round(smooth_lowess(interp1(pk.ind_x_mean, pk.layer(ii).ind_y_flat_mean, 1:block.num_trace, 'linear', 'extrap'), round(pk.length_smooth / nanmean(diff(block.dist))))');
            pk.layer(ii).ind_y_smooth ...
                            = round(smooth_lowess(pk.layer(ii).ind_y, round(pk.length_smooth / nanmean(diff(block.dist))))');
            pk.layer(ii).ind_y_flat_smooth((pk.layer(ii).ind_y_flat_smooth < 1) | (pk.layer(ii).ind_y_flat_smooth > num_sample_trim)) ...
                            = NaN;
            pk.layer(ii).ind_y_smooth((pk.layer(ii).ind_y_smooth < 1) | (pk.layer(ii).ind_y_smooth > num_sample_trim)) ...
                            = NaN;
            if (all(isnan(pk.layer(ii).ind_y_smooth)) || all(isnan(pk.layer(ii).ind_y_flat_smooth)))
                tmp1        = [tmp1 ii]; %#ok<AGROW>
            end
        end
        warning('on', 'MATLAB:interp1:NaNinY')
        
        % remove layers that are empty when smoothed
        axes(ax_radar)
        if ~isempty(tmp1)
            tmp2            = setdiff(1:pk.num_layer, tmp1);
            delete([p_pk(tmp1) p_pkflat(tmp1) p_pkflatmark(tmp1)])
            for ii = tmp1
                if (logical(p_pksmooth(ii)) && ishandle(p_pksmooth(ii)))
                    delete(p_pksmooth(ii))
                end
                if (logical(p_pksmoothflat(ii)) && ishandle(p_pksmoothflat(ii)))
                    delete(p_pksmoothflat(ii))
                end
            end
            curr_layer      = interp1(1:length(tmp2), tmp2, curr_layer, 'nearest', 'extrap');
            [pk.num_layer, pk.layer, p_pkflatmark, p_pkflat, p_pk, smooth_done, p_pksmooth, p_pksmoothflat] ...
                            = deal(length(tmp2), pk.layer(tmp2), p_pkflatmark(tmp2), p_pkflat(tmp2), p_pk(tmp2), smooth_done(tmp2), p_pksmooth(tmp2), p_pksmoothflat(tmp2));
            set(layer_list, 'string', num2cell(1:pk.num_layer), 'value', curr_layer)
            set(status_box, 'string', ['Layer(s) #' num2str(tmp1) ' deleted because they were empty after smoothing.'])
        end
        
        for ii = find(~smooth_done)
            if (logical(p_pksmooth(ii)) && ishandle(p_pksmooth(ii)))
                delete(p_pksmooth(ii))
            end
            if (logical(p_pksmoothflat(ii)) && ishandle(p_pksmoothflat(ii)))
                delete(p_pksmoothflat(ii))
            end
            tmp1            = pk.layer(ii).ind_y_smooth(ind_decim);
            tmp2            = block.dist_lin(ind_decim);
            p_pksmooth(ii)  = plot(tmp2(~isnan(tmp1)), (1e6 .* block.twtt(round(tmp1(~isnan(tmp1))))), 'g.', 'markersize', 12, 'visible', 'off');
            tmp1            = pk.layer(ii).ind_y_flat_smooth(ind_decim);
            p_pksmoothflat(ii) ...
                            = plot(tmp2(~isnan(tmp1)), (1e6 .* block.twtt(round(tmp1(~isnan(tmp1))))), 'g.', 'markersize', 12, 'visible', 'off');
        end
        set(smooth_check, 'value', 1)
        set(status_box, 'string', ['Smoothed ' num2str(length(find(~smooth_done))) ' layer(s) in ' num2str(toc, '%.0f') ' s.'])
        smooth_done(~smooth_done) ...
                            = true;
        match_done          = false;
        show_smooth
    end

%% Match layers in current block to reference layers

    function pk_match(source, eventdata)
        
        if ~pk.num_layer
            set(status_box, 'string', 'No picked layers to match yet.')
            return
        end
        if ~all(smooth_done)
            set(status_box, 'string', 'Not all layers have been smoothed yet.')
            return
        end
        if ~ref_done
            set(status_box, 'string', 'First block of the transect? Y: yes; otherwise: no.')
            waitforbuttonpress
            if strcmpi(get(pkgui, 'currentcharacter'), 'Y')
                pk.ind_match ...
                            = (1:pk.num_layer)';
                match_done  = true;
                set(status_box, 'string', 'Layer numbers set.')
                return
            else
                set(status_box, 'string', 'No reference layers loaded to match to.')
                return
            end
        end
        if strcmp(ref_start_end, 'end')
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
        set(p_pkflat(logical(p_pkflat) & ishandle(p_pkflat)), 'color', [1 0.7 0.7])
        set(p_pk(logical(p_pk) & ishandle(p_pk)), 'color', [1 0.7 0.7])
        set(p_pksmooth(logical(p_pksmooth) & ishandle(p_pksmooth)), 'color', [0.7 1 0.7])
        set(p_pksmoothflat(logical(p_pksmoothflat) & ishandle(p_pksmoothflat)), 'color', [0.7 1 0.7])
        set(p_ref(logical(p_ref) & ishandle(p_ref)), 'color', [1 1 0.7])
        
        % search for matching layers near reference layers
        for ii = 1:pk.num_layer
            
            tmp2            = NaN(pk_ref.num_layer, block.ind_overlap(1));
            tmp2(:, ~isnan(pk.layer(ii).ind_y_smooth(1:block.ind_overlap(1)))) ...
                            = tmp1(:, ~isnan(pk.layer(ii).ind_y_smooth(1:block.ind_overlap(1)))) - ...
                                   repmat(block.twtt(round(pk.layer(ii).ind_y_smooth(~isnan(pk.layer(ii).ind_y_smooth(1:block.ind_overlap(1))))))', pk_ref.num_layer, 1);
            tmp3            = NaN(1, pk_ref.num_layer);
            for jj = 1:pk_ref.num_layer
                tmp3(jj)    = abs(nanmean(tmp2(jj, :)));
            end
            
            if (length(find(tmp3 < pk.twtt_match)) == 1) % found a single matching layer within the threshold
                pk.ind_match(ii) ...
                            = pk_ref.ind_match(logical(tmp3 < pk.twtt_match));
                % brighten layer colors once successful
                set(p_pkflat(ii), 'color', 'r')
                set(p_pk(ii), 'color', 'r')
                set(p_pksmooth(ii), 'color', 'g')
                set(p_pksmoothflat(ii), 'color', 'g')
                set(p_ref(logical(tmp3 < pk.twtt_match)), 'color', 'y')
            end
        end
        
        set(status_box, 'string', ['Matched ' num2str(length(find(~isnan(pk.ind_match)))) '/' num2str(pk.num_layer) ' new layers from ' num2str(length(find(sum(~isnan(tmp1), 2)))) ' overlapping reference layers.'])
        
        % add new layer numbers for those that didn't match any reference layers
        if any(isnan(pk.ind_match))
            pk.ind_match(isnan(pk.ind_match)) ...
                            = (max(pk_ref.ind_match) + 1):(max(pk_ref.ind_match) + length(find(isnan(pk.ind_match))));
        end
        match_done          = true;
    end

%% Save layer picks

    function pk_save(source, eventdata)
        
        % want everything done before saving
        if ~match_done
            set(status_box, 'string', 'Layers not matched yet, or need re-matching.')
            return
        end
        
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
        if ~isempty(path_pk)
            [file_pk, path_pk] = uiputfile('*.mat', 'Save picks:', [path_pk file_data(1:(end - 4)) '_pk.mat']);
        elseif ~isempty(path_data)
            [file_pk, path_pk] = uiputfile('*.mat', 'Save picks:', [path_data file_data(1:(end - 4)) '_pk.mat']);
        else
            [file_pk, path_pk] = uiputfile('*.mat', 'Save picks:', [file_data(1:(end - 4)) '_pk.mat']);
        end
        
        if ~ischar(file_pk)
            [file_pk, path_pk] = deal('', tmp1);
        else
            
            set(status_box, 'string', 'Saving picks...')
            
            tmp1            = pk;
            
            % save many variables in pk structure for easy reference independent of data later on
            [pk.lat, pk.lon, pk.x, pk.y, pk.num_sample, pk.num_trace, pk.file_in, pk.file_block, pk.twtt_min_ref, pk.twtt_max_ref, pk.dist, pk.dist_lin, pk.ind_overlap, pk.elev_air] ...
                            = deal(block.lat, block.lon, block.x, block.y, block.num_sample, block.num_trace, block.file_in, file_data(1:(end - 4)), twtt_min_ref, twtt_max_ref, block.dist, block.dist_lin, ...
                                   block.ind_overlap, block.elev_air);
            if surf_avail
                pk.twtt_surf= block.twtt_surf;
            end
            if bed_avail
                pk.twtt_bed = block.twtt_bed;
            end
            
            % get traveltimes and echo intensities from indices, and adjust indices as appropriate assuming trimming has occurred
            if surf_avail
                pk.elev_surf= block.elev_air - (block.twtt_surf .* (speed_vacuum / 2)); % ice-sheet surface elevation (used to calculate layer elevations)
                if bed_avail
                    pk.elev_bed ...
                            = pk.elev_surf - ((block.twtt_bed - block.twtt_surf) .* (speed_ice / 2)); % bed elevation
                end
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
                [pk.layer(ii).int(~isnan(pk.layer(ii).ind_y)), pk.layer(ii).int_smooth(~isnan(pk.layer(ii).ind_y_smooth))] ...
                            = deal(block.amp(sub2ind([num_sample_trim block.num_trace], round(pk.layer(ii).ind_y(~isnan(pk.layer(ii).ind_y))), find(~isnan(pk.layer(ii).ind_y)))), ...
                                   block.amp(sub2ind([num_sample_trim block.num_trace], round(pk.layer(ii).ind_y_smooth(~isnan(pk.layer(ii).ind_y_smooth))), find(~isnan(pk.layer(ii).ind_y_smooth)))));
                [pk.layer(ii).ind_y, pk.layer(ii).ind_y_smooth, pk.layer(ii).ind_y_flat_mean, pk.layer(ii).ind_y_flat_smooth] ...
                            = deal((pk.layer(ii).ind_y + pk.ind_trim_start), (pk.layer(ii).ind_y_smooth + pk.ind_trim_start), (pk.layer(ii).ind_y_flat_mean + pk.ind_trim_start), ...
                                   (pk.layer(ii).ind_y_flat_smooth + pk.ind_trim_start)); % adjust due to trimming
            end
            
            pk              = orderfields(pk);
            pk.layer        = orderfields(pk.layer);
            
            save([path_pk file_pk], '-v7.3', 'pk')
            
            % make a simple figure that also gets saved
            set(0, 'DefaultFigureWindowStyle', 'default')
            pkfig           = figure('position', [10 10 1600 1000]);
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
            print(pkfig, '-dpng', [path_pk file_pk(1:(end - 4)) '.png'])
            
            pk              = tmp1;
            set(0, 'DefaultFigureWindowStyle', 'docked')
            
            set(status_box, 'string', ['Picks saved as ' file_pk(1:(end - 4)) '.'])
            
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

%% Update minimum dB/phase

    function slide_db_min(source, eventdata)
        switch disp_type
            case {'amp.' 'flat'}
                if (get(cb_min_slide, 'value') < db_max)
                    if get(cbfix_check1, 'value')
                        tmp1 = db_max - db_min;
                    end
                    db_min  = get(cb_min_slide, 'value');
                    if get(cbfix_check1, 'value')
                        db_max = db_min + tmp1;
                        if (db_max > db_max_ref)
                            db_max = db_max_ref;
                            db_min = db_max - tmp1;
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
            case 'phase'
                if (get(cb_min_slide, 'value') < phase_diff_max)
                    if get(cbfix_check1, 'value')
                        tmp1 = phase_diff_max - phase_diff_min;
                    end
                    phase_diff_min ...
                            = get(cb_min_slide, 'value');
                    if get(cbfix_check1, 'value')
                        phase_diff_max = phase_diff_min + tmp1;
                        if (phase_diff_max > phase_diff_max_ref)
                            phase_diff_max = phase_diff_max_ref;
                            phase_diff_min = phase_diff_max - tmp1;
                            if (phase_diff_min < phase_diff_min_ref)
                                phase_diff_min = phase_diff_min_ref;
                            end
                            if (phase_diff_min < get(cb_min_slide, 'min'))
                                set(cb_min_slide, 'value', get(cb_min_slide, 'min'))
                            else
                                set(cb_min_slide, 'value', phase_diff_min)
                            end
                        end
                        set(cb_max_edit, 'string', sprintf('%3.2f', phase_diff_max))
                        if (phase_diff_max > get(cb_max_slide, 'max'))
                            set(cb_max_slide, 'value', get(cb_max_slide, 'max'))
                        else
                            set(cb_max_slide, 'value', phase_diff_max)
                        end
                    end
                    set(cb_min_edit, 'string', sprintf('%3.2f', phase_diff_min))
                    update_db_range
                else
                    if (phase_diff_min < get(cb_min_slide, 'min'))
                        set(cb_min_slide, 'value', get(cb_min_slide, 'min'))
                    else
                        set(cb_min_slide, 'value', phase_diff_min)
                    end
                end
            case 'ARESP'
                if (get(cb_min_slide, 'value') < aresp_max)
                    if get(cbfix_check1, 'value')
                        tmp1 = aresp_max - aresp_min;
                    end
                    aresp_min ...
                            = get(cb_min_slide, 'value');
                    if get(cbfix_check1, 'value')
                        aresp_max = aresp_min + tmp1;
                        if (aresp_max > aresp_max_ref)
                            aresp_max = aresp_max_ref;
                            aresp_min = aresp_max - tmp1;
                            if (aresp_min < aresp_min_ref)
                                aresp_min = aresp_min_ref;
                            end
                            if (aresp_min < get(cb_min_slide, 'min'))
                                set(cb_min_slide, 'value', get(cb_min_slide, 'min'))
                            else
                                set(cb_min_slide, 'value', aresp_min)
                            end
                        end
                        set(cb_max_edit, 'string', sprintf('%3.2f', aresp_max))
                        if (aresp_max > get(cb_max_slide, 'max'))
                            set(cb_max_slide, 'value', get(cb_max_slide, 'max'))
                        else
                            set(cb_max_slide, 'value', aresp_max)
                        end
                    end
                    set(cb_min_edit, 'string', sprintf('%3.2f', aresp_min))
                    update_db_range
                else
                    if (aresp_min < get(cb_min_slide, 'min'))
                        set(cb_min_slide, 'value', get(cb_min_slide, 'min'))
                    else
                        set(cb_min_slide, 'value', aresp_min)
                    end
                end
            case 'spec.'
                if (get(cb_min_slide, 'value') < spec_max)
                    if get(cbfix_check1, 'value')
                        tmp1 = spec_max - spec_min;
                    end
                    spec_min ...
                            = get(cb_min_slide, 'value');
                    if get(cbfix_check1, 'value')
                        spec_max = spec_min + tmp1;
                        if (spec_max > spec_max_ref)
                            spec_max = spec_max_ref;
                            spec_min = spec_max - tmp1;
                            if (spec_min < spec_min_ref)
                                spec_min = spec_min_ref;
                            end
                            if (spec_min < get(cb_min_slide, 'min'))
                                set(cb_min_slide, 'value', get(cb_min_slide, 'min'))
                            else
                                set(cb_min_slide, 'value', spec_min)
                            end
                        end
                        set(cb_max_edit, 'string', sprintf('%3.2f', spec_max))
                        if (spec_max > get(cb_max_slide, 'max'))
                            set(cb_max_slide, 'value', get(cb_max_slide, 'max'))
                        else
                            set(cb_max_slide, 'value', spec_max)
                        end
                    end
                    set(cb_min_edit, 'string', sprintf('%3.2f', spec_min))
                    update_db_range
                else
                    if (spec_min < get(cb_min_slide, 'min'))
                        set(cb_min_slide, 'value', get(cb_min_slide, 'min'))
                    else
                        set(cb_min_slide, 'value', spec_min)
                    end
                end
        end
        set(cb_min_slide, 'enable', 'off')
        drawnow
        set(cb_min_slide, 'enable', 'on')
    end

%% Update maximum dB/phase

    function slide_db_max(source, eventdata)
        switch disp_type
            case {'amp.' 'flat'}
                if (get(cb_max_slide, 'value') > db_min)
                    if get(cbfix_check1, 'value')
                        tmp1 = db_max - db_min;
                    end
                    db_max = get(cb_max_slide, 'value');
                    if get(cbfix_check1, 'value')
                        db_min = db_max - tmp1;
                        if (db_min < db_min_ref)
                            db_min = db_min_ref;
                            db_max = db_min + tmp1;
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
            case 'phase'
                if (get(cb_max_slide, 'value') > phase_diff_min)
                    if get(cbfix_check1, 'value')
                        tmp1= phase_diff_max - phase_diff_min;
                    end
                    phase_diff_max ...
                            = get(cb_max_slide, 'value');
                    if get(cbfix_check1, 'value')
                        phase_diff_min = phase_diff_max - tmp1;
                        if (phase_diff_min < phase_diff_min_ref)
                            phase_diff_min = phase_diff_min_ref;
                            phase_diff_max = phase_diff_min + tmp1;
                            if (phase_diff_max > phase_diff_max_ref)
                                phase_diff_max = phase_diff_max_ref;
                            end
                            if (phase_diff_max > get(cb_max_slide, 'max'))
                                set(cb_max_slide, 'value', get(cb_max_slide, 'max'))
                            else
                                set(cb_max_slide, 'value', phase_diff_max)
                            end
                        end
                        set(cb_min_edit, 'string', sprintf('%3.2f', phase_diff_min))
                        if (phase_diff_min < get(cb_min_slide, 'min'))
                            set(cb_min_slide, 'value', get(cb_min_slide, 'min'))
                        else
                            set(cb_min_slide, 'value', phase_diff_min)
                        end
                    end
                    set(cb_max_edit, 'string', sprintf('%3.2f', phase_diff_max))
                    update_db_range
                else
                    if (phase_diff_max > get(cb_max_slide, 'max'))
                        set(cb_max_slide, 'value', get(cb_max_slide, 'max'))
                    else
                        set(cb_max_slide, 'value', phase_diff_max)
                    end
                end
            case 'ARESP'
                if (get(cb_max_slide, 'value') > aresp_min)
                    if get(cbfix_check1, 'value')
                        tmp1= aresp_max - aresp_min;
                    end
                    aresp_max ...
                            = get(cb_max_slide, 'value');
                    if get(cbfix_check1, 'value')
                        aresp_min = aresp_max - tmp1;
                        if (aresp_min < aresp_min_ref)
                            aresp_min = aresp_min_ref;
                            aresp_max = aresp_min + tmp1;
                            if (aresp_max > aresp_max_ref)
                                aresp_max = aresp_max_ref;
                            end
                            if (aresp_max > get(cb_max_slide, 'max'))
                                set(cb_max_slide, 'value', get(cb_max_slide, 'max'))
                            else
                                set(cb_max_slide, 'value', aresp_max)
                            end
                        end
                        set(cb_min_edit, 'string', sprintf('%3.2f', aresp_min))
                        if (aresp_min < get(cb_min_slide, 'min'))
                            set(cb_min_slide, 'value', get(cb_min_slide, 'min'))
                        else
                            set(cb_min_slide, 'value', aresp_min)
                        end
                    end
                    set(cb_max_edit, 'string', sprintf('%3.2f', aresp_max))
                    update_db_range
                else
                    if (aresp_max > get(cb_max_slide, 'max'))
                        set(cb_max_slide, 'value', get(cb_max_slide, 'max'))
                    else
                        set(cb_max_slide, 'value', aresp_max)
                    end
                end
            case 'spec.'
                if (get(cb_max_slide, 'value') > spec_min)
                    if get(cbfix_check1, 'value')
                        tmp1 = spec_max - spec_min;
                    end
                    spec_max ...
                        = get(cb_max_slide, 'value');
                    if get(cbfix_check1, 'value')
                        spec_min = spec_max - tmp1;
                        if (spec_min < spec_min_ref)
                            spec_min = spec_min_ref;
                            spec_max = spec_min + tmp1;
                            if (spec_max > spec_max_ref)
                                spec_max = spec_max_ref;
                            end
                            if (spec_max > get(cb_max_slide, 'max'))
                                set(cb_max_slide, 'value', get(cb_max_slide, 'max'))
                            else
                                set(cb_max_slide, 'value', spec_max)
                            end
                        end
                        set(cb_min_edit, 'string', sprintf('%3.2f', spec_min))
                        if (spec_min < get(cb_min_slide, 'min'))
                            set(cb_min_slide, 'value', get(cb_min_slide, 'min'))
                        else
                            set(cb_min_slide, 'value', spec_min)
                        end
                    end
                    set(cb_max_edit, 'string', sprintf('%3.2f', spec_max))
                    update_db_range
                else
                    if (spec_max > get(cb_max_slide, 'max'))
                        set(cb_max_slide, 'value', get(cb_max_slide, 'max'))
                    else
                        set(cb_max_slide, 'value', spec_max)
                    end
                end
        end
        set(cb_max_slide, 'enable', 'off')
        drawnow
        set(cb_max_slide, 'enable', 'on')
    end

%% Reset minimum dB/phase

    function reset_db_min(source, eventdata)
        switch disp_type
            case {'amp.' 'flat'}
                if (db_min_ref < get(cb_min_slide, 'min'))
                    set(cb_min_slide, 'value', get(cb_min_slide, 'min'))
                else
                    set(cb_min_slide, 'value', db_min_ref)
                end
                set(cb_min_edit, 'string', num2str(db_min_ref))
                db_min      = db_min_ref;
            case 'phase'
                if (phase_diff_min_ref < get(cb_min_slide, 'min'))
                    set(cb_min_slide, 'value', get(cb_min_slide, 'min'))
                else
                    set(cb_min_slide, 'value', phase_diff_min_ref)
                end
                set(cb_min_edit, 'string', num2str(phase_diff_min_ref))
                phase_diff_min ...
                            = phase_diff_min_ref;
            case 'ARESP'
                if (aresp_min_ref < get(cb_min_slide, 'min'))
                    set(cb_min_slide, 'value', get(cb_min_slide, 'min'))
                else
                    set(cb_min_slide, 'value', aresp_min_ref)
                end
                set(cb_min_edit, 'string', num2str(aresp_min_ref))
                aresp_min   = aresp_min_ref;
            case 'spec.'
                if (spec_min_ref < get(cb_min_slide, 'min'))
                    set(cb_min_slide, 'value', get(cb_min_slide, 'min'))
                else
                    set(cb_min_slide, 'value', spec_min_ref)
                end
                set(cb_min_edit, 'string', num2str(spec_min_ref))
                spec_min    = phase_diff_min_ref;
        end
        update_db_range
    end

%% Reset maximum dB/phase

    function reset_db_max(source, eventdata)
        switch disp_type
            case {'amp.' 'flat'}
                if (db_max_ref > get(cb_max_slide, 'max'))
                    set(cb_max_slide, 'value', get(cb_max_slide, 'max'))
                else
                    set(cb_max_slide, 'value', db_max_ref)
                end
                set(cb_max_edit, 'string', num2str(db_max_ref))
                db_max      = db_max_ref;
            case 'phase'
                if (phase_diff_max_ref > get(cb_max_slide, 'max'))
                    set(cb_max_slide, 'value', get(cb_max_slide, 'max'))
                else
                    set(cb_max_slide, 'value', phase_diff_max_ref)
                end
                set(cb_max_edit, 'string', num2str(phase_diff_max_ref))
                phase_diff_max ...
                            = phase_diff_max_ref;
            case 'ARESP'
                if (aresp_max_ref > get(cb_max_slide, 'max'))
                    set(cb_max_slide, 'value', get(cb_max_slide, 'max'))
                else
                    set(cb_max_slide, 'value', aresp_max_ref)
                end
                set(cb_max_edit, 'string', num2str(aresp_max_ref))
                aresp_max   = aresp_max_ref;
            case 'spec.'
                if (spec_max_ref > get(cb_max_slide, 'max'))
                    set(cb_max_slide, 'value', get(cb_max_slide, 'max'))
                else
                    set(cb_max_slide, 'value', spec_max_ref)
                end
                set(cb_max_edit, 'string', num2str(spec_max_ref))
                spec_max    = spec_max_ref;
        end
        update_db_range
    end

%% Update dB/phase range

    function update_db_range(source, eventdata)
        axes(ax_radar)
        switch disp_type
            case {'amp.' 'flat'}
                caxis([db_min db_max])
            case 'phase'
                caxis([phase_diff_min phase_diff_max])
            case 'ARESP'
                caxis([aresp_min aresp_max])
            case 'spec.'
                caxis([spec_min spec_max])
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
            reset_dist_min;
        else
            if (tmp1(1) < get(dist_min_slide, 'min'))
                set(dist_min_slide, 'value', get(dist_min_slide, 'min'))
            else
                set(dist_min_slide, 'value', tmp1(1))
            end
            set(dist_min_edit, 'string', sprintf('%3.1f', tmp1(1)))
            dist_min        = tmp1(1);
        end
        if (tmp1(2) > dist_max_ref);
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

%% Switch display type

    function disp_radio(~, eventdata)
        disp_type           = get(eventdata.NewValue, 'string');
        switch disp_type
            case 'amp.'
                plot_db
            case 'phase'
                plot_phase_diff
            case 'ARESP'
                plot_aresp
            case 'spec.'
                plot_spec
            case 'flat'
                plot_flat
        end
    end

%% Plot data in dB (amplitude)

    function plot_db(source, eventdata)
        if ~load_done
            set(status_box, 'string', 'Data not yet loaded.')
            return
        end
        if (logical(p_data) && ishandle(p_data)) % get rid of old plotted data
            delete(p_data)
        end
        if (any(p_pkflat) && any(ishandle(p_pkflat)))
            set(p_pkflat(logical(p_pkflat) & ishandle(p_pkflat)), 'visible', 'off')
        end
        if (any(p_pksmoothflat) && any(ishandle(p_pksmoothflat)))
            set(p_pksmoothflat(logical(p_pksmoothflat) & ishandle(p_pksmoothflat)), 'visible', 'off')
        end
        if (any(p_refflat) && any(ishandle(p_refflat)))
            set(p_refflat(logical(p_refflat) & ishandle(p_refflat)), 'visible', 'off')
        end
        if (any(p_pkflatmark) && any(ishandle(p_pkflatmark)))
            set(p_pkflatmark(logical(p_pkflatmark) & ishandle(p_pkflatmark)), 'visible', 'off')
        end
        if (any(p_flat) && any(ishandle(p_flat)))
            set(p_flat(logical(p_flat) & ishandle(p_flat)), 'visible', 'off')
        end
        axes(ax_radar)
        if (get(cmap_list, 'value') ~= 1)
            set(cmap_list, 'value', 1)
            change_cmap
        end
        p_data              = imagesc(block.dist_lin(ind_decim), (1e6 .* block.twtt), amp_mean, [db_min db_max]);
        if (logical(p_surf) && ishandle(p_surf))
            set(p_surf, 'visible', 'on')
            uistack(p_surf, 'top')
        elseif surf_avail
            p_surf          = plot(block.dist_lin(ind_decim), (1e6 .* block.twtt_surf(ind_decim)), 'g--', 'linewidth', 2);
            if any(isnan(ind_surf))
                set(p_surf, 'marker', '.', 'linestyle', 'none', 'markersize', 12)
            end
        end
        if (logical(p_bed) && ishandle(p_bed))
            set(p_bed, 'visible', 'on')
            uistack(p_bed, 'top')
        elseif bed_avail
            p_bed           = plot(block.dist_lin(ind_decim), (1e6 .* block.twtt_bed(ind_decim)), 'g--', 'linewidth', 2);
            if any(isnan(ind_bed))
                set(p_bed, 'marker', '.', 'linestyle', 'none', 'markersize', 12)
            end
        end
        disp_type           = 'amp.';
        narrow_cb
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
            disp_type       = 'amp.';
            plot_db
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
        if (logical(p_surf) && ishandle(p_surf))
            set(p_surf, 'visible', 'on')
            uistack(p_surf, 'top')
        end
        if (logical(p_bed) && ishandle(p_bed))
            set(p_bed, 'visible', 'on')
            uistack(p_bed, 'top')
        end
        set(cbl, 'string', '(rad)')
        set(cb_min_slide, 'min', phase_diff_min_ref, 'max', phase_diff_max_ref, 'value', phase_diff_min)
        set(cb_max_slide, 'min', phase_diff_min_ref, 'max', phase_diff_max_ref, 'value', phase_diff_max)
        set(cb_min_edit, 'string', sprintf('%2.1f', phase_diff_min))
        set(cb_max_edit, 'string', sprintf('%2.1f', phase_diff_max))
        disp_type           = 'phase';
        show_phase
        narrow_cb
    end

%% Plot ARESP-derived layer slope

    function plot_aresp(source, eventdata)
        if ~aresp_avail
            set(disp_group, 'selectedobject', disp_check(1))
            disp_type       = 'amp.';
            plot_db
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
        p_data              = imagesc(block.dist_lin(1:decim:size(block.slope_aresp, 2)), (1e6 .* block.twtt), block.slope_aresp(:, 1:decim:end), [aresp_min aresp_max]);
        if (logical(p_surf) && ishandle(p_surf))
            set(p_surf, 'visible', 'on')
            uistack(p_surf, 'top')
        end
        if (logical(p_bed) && ishandle(p_bed))
            set(p_bed, 'visible', 'on')
            uistack(p_bed, 'top')
        end
        set(cbl, 'string', '(tan \theta)')
        set(cb_min_slide, 'min', aresp_min_ref, 'max', aresp_max_ref, 'value', aresp_min)
        set(cb_max_slide, 'min', aresp_min_ref, 'max', aresp_max_ref, 'value', aresp_max)
        set(cb_min_edit, 'string', sprintf('%2.1f', aresp_min))
        set(cb_max_edit, 'string', sprintf('%2.1f', aresp_max))
        disp_type           = 'ARESP';
        show_aresp
        narrow_cb
    end

%% Plot specularity

    function plot_spec(source, eventdata)
        if ~spec_avail
            set(disp_group, 'selectedobject', disp_check(1))
            disp_type       = 'amp.';
            plot_db
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
        p_data              = imagesc(block.dist_lin(1:decim:size(block.phase_diff_filt, 2)), (1e6 .* block.twtt), block.specularity(:, 1:decim:end), [spec_min spec_max]);
        if (logical(p_surf) && ishandle(p_surf))
            set(p_surf, 'visible', 'on')
            uistack(p_surf, 'top')
        elseif surf_avail
            p_surf          = plot(block.dist_lin(ind_decim), (1e6 .* block.twtt_surf(ind_decim)), 'g--', 'linewidth', 2);
            if any(isnan(ind_surf))
                set(p_surf, 'marker', '.', 'linestyle', 'none', 'markersize', 12)
            end
        end
        if (logical(p_bed) && ishandle(p_bed))
            set(p_bed, 'visible', 'on')
            uistack(p_bed, 'top')
        elseif bed_avail
            p_bed           = plot(block.dist_lin(ind_decim), (1e6 .* block.twtt_bed(ind_decim)), 'g--', 'linewidth', 2);
            if any(isnan(ind_bed))
                set(p_bed, 'marker', '.', 'linestyle', 'none', 'markersize', 12)
            end
        end
        set(cbl, 'string', '()')
        set(cb_min_slide, 'min', spec_min_ref, 'max', spec_max_ref, 'value', spec_min)
        set(cb_max_slide, 'min', spec_min_ref, 'max', spec_max_ref, 'value', spec_max)
        set(cb_min_edit, 'string', sprintf('%3.2f', spec_min))
        set(cb_max_edit, 'string', sprintf('%3.2f', spec_max))
        disp_type       = 'spec.';
        narrow_cb
    end

%% Plot layer-flattened radargram

    function plot_flat(source, eventdata)
        if ~flat_done
            set(disp_group, 'selectedobject', disp_check(1))
            disp_type       = 'amp.';
            plot_db
            return
        end
        if (logical(p_data) && ishandle(p_data)) % get rid of old plotted data
            delete(p_data)
        end
        if (logical(p_surf) && ishandle(p_surf))
            set(p_surf, 'visible', 'off')
        end
        if (logical(p_bed) && ishandle(p_bed))
            set(p_bed, 'visible', 'off')
        end
        if (logical(p_startphase) && ishandle(p_startphase))
            set(p_startphase, 'visible', 'off')
        end
        if (logical(p_startaresp) && ishandle(p_startaresp))
            set(p_startaresp, 'visible', 'off')
        end
        if (any(p_phase) && any(ishandle(p_phase)))
            set(p_phase(logical(p_phase) & ishandle(p_phase)), 'visible', 'off')
        end
        if (any(p_aresp) && any(ishandle(p_aresp)))
            set(p_aresp(logical(p_aresp) & ishandle(p_aresp)), 'visible', 'off')
        end
        if (any(p_ref) && any(ishandle(p_ref)))
            set(p_ref(logical(p_ref) & ishandle(p_ref)), 'visible', 'off')
        end
        if (any(p_man(:)) && any(ishandle(p_man(:))))
            set(p_man(logical(p_man) & ishandle(p_man)), 'visible', 'off')
        end
        if (any(p_pk) && any(ishandle(p_pk)))
            set(p_pk(logical(p_pk) & ishandle(p_pk)), 'visible', 'off')
        end
        if (any(p_pksmooth) && any(ishandle(p_pksmooth)))
            set(p_pksmooth(logical(p_pksmooth) & ishandle(p_pksmooth)), 'visible', 'off')
        end
        if (logical(p_surfflat) && ishandle(p_surfflat))
            delete(p_surfflat)
        end
        if (logical(p_bedflat) && ishandle(p_bedflat))
            delete(p_bedflat)
        end
        axes(ax_radar)
        if (get(cmap_list, 'value') ~= 1)
            set(cmap_list, 'value', 1)
            change_cmap
        end
        switch flat_switch
            case 'full'
                p_data      = imagesc(block.dist_lin(ind_decim), (1e6 .* block.twtt), amp_flat(:, ind_decim), [db_min db_max]);
                if surf_avail
                    if any(isnan(ind_surf_flat))
                        tmp1        = ind_surf_flat(ind_decim);
                        p_surfflat  = plot(block.dist_lin(ind_decim(~isnan(tmp1))), (1e6 .* block.twtt(tmp1(~isnan(tmp1)))), 'g.', 'markersize', 12);
                    else
                        p_surfflat  = plot(block.dist_lin(ind_decim), (1e6 .* block.twtt(ind_surf_flat(ind_decim))), 'g--', 'linewidth', 2);
                    end
                end
                if bed_avail
                    if any(isnan(ind_bed_flat))
                        tmp1        = ind_bed_flat(ind_decim);
                        p_bedflat   = plot(block.dist_lin(ind_decim(~isnan(tmp1))), (1e6 .* block.twtt(tmp1(~isnan(tmp1)))), 'g.', 'markersize', 12);
                    else
                        p_bedflat   = plot(block.dist_lin(ind_decim), (1e6 .* block.twtt(ind_bed_flat(ind_decim))), 'g--', 'linewidth', 2);
                    end
                end
            case 'mean'
                p_data  = imagesc(block.dist_lin(pk.ind_x_mean), (1e6 .* block.twtt), amp_flat_mean, [db_min db_max]);
                if surf_avail
                    if any(isnan(ind_surf_flat))
                        tmp1        = ind_surf_flat(pk.ind_x_mean);
                        p_surfflat  = plot(block.dist_lin(pk.ind_x_mean(~isnan(tmp1))), (1e6 .* block.twtt(tmp1(~isnan(tmp1)))), 'g.', 'markersize', 12);
                    else
                        p_surfflat  = plot(block.dist_lin(pk.ind_x_mean), (1e6 .* block.twtt(ind_surf_flat(pk.ind_x_mean))), 'g--', 'linewidth', 2);
                    end
                end
                if bed_avail
                    if any(isnan(ind_bed_flat))
                        tmp1        = ind_bed_flat(pk.ind_x_mean);
                        p_bedflat   = plot(block.dist_lin(pk.ind_x_mean(~isnan(tmp1))), (1e6 .* block.twtt(tmp1(~isnan(tmp1)))), 'g.', 'markersize', 12);
                    else
                        p_bedflat   = plot(block.dist_lin(pk.ind_x_mean), (1e6 .* block.twtt(ind_bed_flat(pk.ind_x_mean))), 'g--', 'linewidth', 2);
                    end
                end
        end
        disp_type       = 'flat';
        narrow_cb
        show_flat
        show_ref
        show_pk
        show_smooth
        set(cbl, 'string', '(dB)')
        set(cb_min_slide, 'min', db_min_ref, 'max', db_max_ref, 'value', db_min)
        set(cb_max_slide, 'min', db_min_ref, 'max', db_max_ref, 'value', db_max)
        set(cb_min_edit, 'string', sprintf('%3.0f', db_min))
        set(cb_max_edit, 'string', sprintf('%3.0f', db_max))
    end

%% Show phase-tracked layers

    function show_phase(source, eventdata)
        if phase_done
            if (get(phase_check, 'value') && strcmp(disp_type, 'amp.'))
                if (any(p_startphase) && any(ishandle(p_startphase)))
                    set(p_startphase(logical(p_startphase) & ishandle(p_startphase)), 'visible', 'on')
                    uistack(p_startphase(logical(p_startphase) & ishandle(p_startphase)), 'top')
                end
                if (any(p_phase) && any(ishandle(p_phase)))
                    set(p_phase(logical(p_phase) & ishandle(p_phase)), 'visible', 'on')
                    uistack(p_phase(logical(p_phase) & ishandle(p_phase)), 'top')
                end
            else
                if (any(p_startphase) && any(ishandle(p_startphase)))
                    set(p_startphase(logical(p_startphase) & ishandle(p_startphase)), 'visible', 'off')
                end
                if (any(p_phase) && any(ishandle(p_phase)))
                    set(p_phase(logical(p_phase) & ishandle(p_phase)), 'visible', 'off')
                end
            end
        elseif get(phase_check, 'value')
            set(phase_check, 'value', 0)
        end
    end

%% Show ARESP-tracked layers

    function show_aresp(source, eventdata)
        if aresp_done
            if (get(aresp_check, 'value') && strcmp(disp_type, 'amp.'))
                if (any(p_startaresp) && any(ishandle(p_startaresp)))
                    set(p_startaresp(logical(p_startaresp) & ishandle(p_startaresp)), 'visible', 'on')
                    uistack(p_startaresp(logical(p_startaresp) & ishandle(p_startaresp)), 'top')
                end
                if (any(p_aresp) && any(ishandle(p_aresp)))
                    set(p_aresp(logical(p_aresp) & ishandle(p_aresp)), 'visible', 'on')
                    uistack(p_aresp(logical(p_aresp) & ishandle(p_aresp)), 'top')
                end
            else
                if (any(p_startaresp) && any(ishandle(p_startaresp)))
                    set(p_startaresp(logical(p_startaresp) & ishandle(p_startaresp)), 'visible', 'off')
                end
                if (any(p_aresp) && any(ishandle(p_aresp)))
                    set(p_aresp(logical(p_aresp) & ishandle(p_aresp)), 'visible', 'off')
                end
            end
        elseif get(aresp_check, 'value')
            set(aresp_check, 'value', 0)
        end
    end

%% Show flattened predictions

    function show_flat(source, eventdata)
        if flat_done
            if (get(flat_check, 'value') && strcmp(disp_type, 'flat'))
                if (any(p_flat) && any(ishandle(p_flat)))
                    set(p_flat(logical(p_flat) & ishandle(p_flat)), 'visible', 'on')
                    uistack(p_flat(logical(p_flat) & ishandle(p_flat)), 'top')
                end
            else
                if (any(p_flat) && any(ishandle(p_flat)))
                    set(p_flat(logical(p_flat) & ishandle(p_flat)), 'visible', 'off')
                end
            end
        elseif get(flat_check, 'value')
            set(flat_check, 'value', 0)
        end
    end

%% Show horizontally averaged, flattened data

    function show_mean(source, eventdata)
        if mean_done
            if get(mean_check, 'value')
                flat_switch = 'mean';
            else
                flat_switch = 'full';
            end
            set(disp_group, 'selectedobject', disp_check(5))
            disp_type       = 'flat';
            plot_flat
        elseif get(mean_check, 'value')
            set(mean_check, 'value', 0)
            flat_switch     = 'full';
        end
    end

%% Show picked layers

    function show_pk(source, eventdata)
        if pk_done
            if get(pk_check, 'value')
                switch disp_type
                    case 'amp.'
                        if (any(p_pk) && any(ishandle(p_pk)))
                            set(p_pk(logical(p_pk) & ishandle(p_pk)), 'visible', 'on')
                            uistack(p_pk(logical(p_pk) & ishandle(p_pk)), 'top')
                        end
                        if (any(p_pkflat) && any(ishandle(p_pkflat)))
                            set(p_pkflat(logical(p_pkflat) & ishandle(p_pkflat)), 'visible', 'off')
                        end
                        if (any(p_pkflatmark) && any(ishandle(p_pkflatmark)))
                            set(p_pkflatmark(logical(p_pkflatmark) & ishandle(p_pkflatmark)), 'visible', 'off')
                        end
                    case 'flat'
                        if (any(p_pkflat) && any(ishandle(p_pkflat)))
                            set(p_pkflat(logical(p_pkflat) & ishandle(p_pkflat)), 'visible', 'on')
                            uistack(p_pkflat(logical(p_pkflat) & ishandle(p_pkflat)), 'top')
                        end
                        if (any(p_pkflatmark) && any(ishandle(p_pkflatmark)))
                            set(p_pkflatmark(logical(p_pkflatmark) & ishandle(p_pkflatmark)), 'visible', 'on')
                            uistack(p_pkflatmark(logical(p_pkflatmark) & ishandle(p_pkflatmark)), 'top')
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
                if (any(p_pkflatmark) && any(ishandle(p_pkflatmark)))
                    set(p_pkflatmark(logical(p_pkflatmark) & ishandle(p_pkflatmark)), 'visible', 'off')
                end
            end
        elseif get(pk_check, 'value')
            set(pk_check, 'value', 0)
        end
    end

%% Show smoothed layers

    function show_smooth(source, eventdata)
        if any(smooth_done)
            if get(smooth_check, 'value')
                switch disp_type
                    case 'amp.'
                        if (any(p_pksmooth) && any(ishandle(p_pksmooth)))
                            set(p_pksmooth(smooth_done), 'visible', 'on')
                            uistack(p_pksmooth(smooth_done), 'top')
                            set(p_pksmoothflat(smooth_done), 'visible', 'off')
                        end
                    case 'flat'
                        if (any(p_pksmoothflat) && any(ishandle(p_pksmoothflat)))
                            set(p_pksmoothflat(smooth_done), 'visible', 'on')
                            uistack(p_pksmoothflat(smooth_done), 'top')
                            set(p_pksmooth(smooth_done), 'visible', 'off')
                        end
                end
            else
                if (any(p_pksmooth) && any(ishandle(p_pksmooth)))
                    set(p_pksmooth(smooth_done), 'visible', 'off')
                end
                if (any(p_pksmoothflat) && any(ishandle(p_pksmoothflat)))
                    set(p_pksmoothflat(smooth_done), 'visible', 'off')
                end
            end
        elseif get(smooth_check, 'value')
            set(smooth_check, 'value', 0)
        end
    end

%% Show reference layers

    function show_ref(source, eventdata)
        if ref_done
            if get(ref_check, 'value')
                switch disp_type
                    case 'amp.'
                        if (any(p_ref) && any(ishandle(p_ref)))
                            set(p_ref(logical(p_ref) & ishandle(p_ref)), 'visible', 'on')
                            uistack(p_ref(logical(p_ref) & ishandle(p_ref)), 'top')
                        end
                        if (any(p_refflat) && any(ishandle(p_refflat)))
                            set(p_refflat(logical(p_refflat) & ishandle(p_refflat)), 'visible', 'off')
                        end
                    case 'flat'
                        if (any(p_refflat) && any(ishandle(p_refflat)))
                            set(p_refflat(logical(p_refflat) & ishandle(p_refflat)), 'visible', 'on')
                            uistack(p_refflat(logical(p_refflat) & ishandle(p_refflat)), 'top')
                        end
                        if (any(p_ref) && any(ishandle(p_ref)))
                            set(p_ref(logical(p_ref) & ishandle(p_ref)), 'visible', 'off')
                        end
                end
            else
                set(p_ref(logical(p_ref) & ishandle(p_ref)), 'visible', 'off')
                set(p_refflat(logical(p_refflat) & ishandle(p_refflat)), 'visible', 'off')
            end
        elseif get(ref_check, 'value')
            set(ref_check, 'value', 0)
        end
    end

%% Show manual layers

    function show_man(source, eventdata)
        if pk.num_man
            if get(man_check, 'value')
                if (strcmp(disp_type, 'amp.') && any(p_man(:)) && any(ishandle(p_man(:))))
                    set(p_man(logical(p_man(:)) & ishandle(p_man(:))), 'visible', 'on')
                    uistack(p_man(logical(p_man(:)) & ishandle(p_man(:))), 'top')
                end
            else
                if (any(p_man(:)) && any(ishandle(p_man(:))))
                    set(p_man(logical(p_man(:)) & ishandle(p_man(:))), 'visible', 'off')
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
            if (decim > 1)
                amp_mean    = NaN(num_sample_trim, length(ind_decim), 'single');
                tmp1        = floor(decim / 2);
                for ii = 1:length(ind_decim)
                    amp_mean(:, ii) ...
                            = nanmean(block.amp(:, (ind_decim(ii) - tmp1):(ind_decim(ii) + tmp1)), 2);
                end
            else
                amp_mean    = single(block.amp);
            end
            switch disp_type
                case 'amp.'
                    plot_db
                case 'phase'
                    plot_phase_diff
                case 'ARESP'
                    plot_aresp
                case 'spec.'
                    plot_spec
                case 'flat'
                    plot_flat
            end
            set(status_box, 'string', ['Displayed data decimated to 1/' num2str(decim) ' samples.'])
        else
            set(status_box, 'string', ['Decimation set to 1/' num2str(decim) ' samples.'])
        end
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
            axes(ax_radar) %#ok<*LAXES>
            p_phase         = zeros(1, pk.num_phase);
            for ii = 1:pk.num_phase
                p_phase(ii)     = plot(block.dist_lin(ind_decim), (1e6 .* block.twtt(round(ind_y_phase(ii, ind_decim)))), 'b', 'linewidth', 1); % plot the phase-tracked layers
            end
            if keep_phase_done
                set(p_phase(pk.ind_keep_phase), 'color', 'w', 'linewidth', 2)
            end
            set(status_box, 'string', ['Adjusted phase-tracked layers to ' num2str(1e-6 * pk.freq) ' MHz'])
        end
    end

%% Adjust number of indices over which to horizontally average

    function adj_num_ind_mean(source, eventdata)
        pk.num_ind_mean     = abs(round(str2double(get(num_ind_mean_edit, 'string'))));
        set(status_box, 'string', ['Horizontal averaging length adjusted to ' num2str(pk.num_ind_mean) ' samples.'])
        set(num_ind_mean_edit, 'string', num2str(pk.num_ind_mean))
        if mean_done
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
        set(status_box, 'string', ['Vertical search window adjusted to +/- ' num2str(pk.num_win) ' samples.'])
        set(num_win_edit, 'string', num2str(pk.num_win))
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
    end

%% Adjust matching range

    function adj_twtt_match(source, eventdata)
        pk.twtt_match       = 1e-6 * abs(str2double(get(twtt_match_edit, 'string')));
        set(status_box, 'string', ['Matching traveltime threshold set to ' num2str(1e6 * pk.twtt_match) ' us.'])
        set(twtt_match_edit, 'string', num2str(1e6 * pk.twtt_match))
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
        figure('position', [10 10 1200 800]);
        axis ij tight
        axis([dist_min dist_max (1e6 .* [twtt_min twtt_max])])
        hold on
        colormap(cmaps{get(cmap_list, 'value')})
        switch disp_type
            case 'amp.'
                imagesc(block.dist_lin(ind_decim), (1e6 .* block.twtt), amp_mean, [db_min db_max])
                caxis([db_min db_max])
                if surf_avail
                    tmp1    = plot(block.dist_lin(ind_decim), (1e6 .* block.twtt_surf(ind_decim)), 'g--', 'linewidth', 2);
                    if any(isnan(ind_surf))
                        set(tmp1, 'marker', '.', 'linestyle', 'none', 'markersize', 12);
                    end
                end
                if bed_avail
                    tmp1    = plot(block.dist_lin(ind_decim), (1e6 .* block.twtt_bed(ind_decim)), 'g--', 'linewidth', 2);
                    if any(isnan(ind_bed))
                        set(tmp1, 'marker', '.', 'linestyle', 'none', 'markersize', 12)
                    end
                end
                if get(ref_check, 'value')
                    switch ref_start_end
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
                        tmp1= plot(block.dist_lin(ind_decim), (1e6 .* block.twtt(round(ind_y_aresp(ii, ind_decim)))), 'b', 'linewidth', 1);
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
                        if ~isempty(ind_decim(~isnan(pk.layer(ii).ind_y(ind_decim))))
                            tmp1 = plot(block.dist_lin(ind_decim(~isnan(pk.layer(ii).ind_y(ind_decim)))), (1e6 .* block.twtt(round(pk.layer(ii).ind_y(ind_decim(~isnan(pk.layer(ii).ind_y(ind_decim))))))), ...
                                        'r.', 'markersize', 12);
                            if (ii == curr_layer)
                                set(tmp1, 'markersize', 24)
                            end
                        end
                    end
                end
                if get(smooth_check, 'value')
                    for ii = 1:pk.num_layer
                        if ~isempty(ind_decim(~isnan(pk.layer(ii).ind_y_smooth(ind_decim))))
                            tmp1    = plot(block.dist_lin(ind_decim(~isnan(pk.layer(ii).ind_y_smooth(ind_decim)))), ...
                                           (1e6 .* block.twtt(round(pk.layer(ii).ind_y_smooth(ind_decim(~isnan(pk.layer(ii).ind_y_smooth(ind_decim))))))), 'g.', 'markersize', 12);
                        end
                        if (ii == curr_layer)
                            set(tmp1, 'markersize', 24)
                        end
                    end
                end
                text((dist_max + (0.015 * (dist_max - dist_min))), (1e6 * (twtt_min - (0.04 * (twtt_max - twtt_min)))), '(dB)', 'color', 'k', 'fontsize', 20)
            case 'phase'
                imagesc(block.dist_lin(1:decim:size(block.phase_diff_filt, 2)), (1e6 .* block.twtt), block.phase_diff_filt(:, 1:decim:end), [phase_diff_min phase_diff_max])
                if surf_avail
                    tmp1    = plot(block.dist_lin(ind_decim), (1e6 .* block.twtt_surf(ind_decim)), 'g--', 'linewidth', 2);
                    if any(isnan(ind_surf))
                        set(tmp1, 'marker', '.', 'linestyle', 'none', 'markersize', 12);
                    end
                end
                if bed_avail
                    tmp1    = plot(block.dist_lin(ind_decim), (1e6 .* block.twtt_bed(ind_decim)), 'g--', 'linewidth', 2);
                    if any(isnan(ind_bed))
                        set(tmp1, 'marker', '.', 'linestyle', 'none', 'markersize', 12)
                    end
                end
                text((dist_max + (0.015 * (dist_max - dist_min))), (1e6 * (twtt_min - (0.04 * (twtt_max - twtt_min)))), '(rad)', 'color', 'k', 'fontsize', 20)
            case 'ARESP'
                imagesc(block.dist_lin(1:decim:size(block.slope_aresp, 2)), (1e6 .* block.twtt), block.slope_aresp(:, 1:decim:end), [aresp_min aresp_max])
                if surf_avail
                    tmp1    = plot(block.dist_lin(ind_decim), (1e6 .* block.twtt_surf(ind_decim)), 'g--', 'linewidth', 2);
                    if any(isnan(ind_surf))
                        set(tmp1, 'marker', '.', 'linestyle', 'none', 'markersize', 12);
                    end
                end
                if bed_avail
                    tmp1    = plot(block.dist_lin(ind_decim), (1e6 .* block.twtt_bed(ind_decim)), 'g--', 'linewidth', 2);
                    if any(isnan(ind_bed))
                        set(tmp1, 'marker', '.', 'linestyle', 'none', 'markersize', 12)
                    end
                end
                text((dist_max + (0.015 * (dist_max - dist_min))), (1e6 * (twtt_min - (0.04 * (twtt_max - twtt_min)))), '(tan \theta)', 'color', 'k', 'fontsize', 20)
            case 'spec.'
                imagesc(block.dist_lin(1:decim:size(block.phase_diff_filt, 2)), (1e6 .* block.twtt), block.specularity(:, 1:decim:end), [spec_min spec_max])
                if surf_avail
                    tmp1    = plot(block.dist_lin(ind_decim), (1e6 .* block.twtt_surf(ind_decim)), 'g--', 'linewidth', 2);
                    if any(isnan(ind_surf))
                        set(tmp1, 'marker', '.', 'linestyle', 'none', 'markersize', 12);
                    end
                end
                if bed_avail
                    tmp1    = plot(block.dist_lin(ind_decim), (1e6 .* block.twtt_bed(ind_decim)), 'g--', 'linewidth', 2);
                    if any(isnan(ind_bed))
                        set(tmp1, 'marker', '.', 'linestyle', 'none', 'markersize', 12)
                    end
                end
            case 'flat'
                switch flat_switch
                    case 'full'
                        imagesc(block.dist_lin(ind_decim), (1e6 .* block.twtt), amp_flat(:, ind_decim), [db_min db_max]);
                        if surf_avail
                            if any(isnan(ind_surf_flat))
                                tmp1        = ind_surf_flat(ind_decim);
                                plot(block.dist_lin(ind_decim(~isnan(tmp1))), (1e6 .* block.twtt(tmp1(~isnan(tmp1)))), 'g.', 'markersize', 12);
                            else
                                plot(block.dist_lin(ind_decim), (1e6 .* block.twtt(ind_surf_flat(ind_decim))), 'g--', 'linewidth', 2);
                            end
                        end
                        if bed_avail
                            if any(isnan(ind_bed_flat))
                                tmp1        = ind_bed_flat(ind_decim);
                                plot(block.dist_lin(ind_decim(~isnan(tmp1))), (1e6 .* block.twtt(tmp1(~isnan(tmp1)))), 'g.', 'markersize', 12);
                            else
                                plot(block.dist_lin(ind_decim), (1e6 .* block.twtt(ind_bed_flat(ind_decim))), 'g--', 'linewidth', 2);
                            end
                        end
                    case 'mean'
                        imagesc(block.dist_lin(pk.ind_x_mean), (1e6 .* block.twtt), amp_flat_mean, [db_min db_max]);
                        if surf_avail
                            if any(isnan(ind_surf_flat))
                                tmp1        = ind_surf_flat(pk.ind_x_mean);
                                plot(block.dist_lin(pk.ind_x_mean(~isnan(tmp1))), (1e6 .* block.twtt(tmp1(~isnan(tmp1)))), 'g.', 'markersize', 12);
                            else
                                plot(block.dist_lin(pk.ind_x_mean), (1e6 .* block.twtt(ind_surf_flat(pk.ind_x_mean))), 'g--', 'linewidth', 2);
                            end
                        end
                        if bed_avail
                            if any(isnan(ind_bed_flat))
                                tmp1        = ind_bed_flat(pk.ind_x_mean);
                                plot(block.dist_lin(pk.ind_x_mean(~isnan(tmp1))), (1e6 .* block.twtt(tmp1(~isnan(tmp1)))), 'g.', 'markersize', 12);
                            else
                                plot(block.dist_lin(pk.ind_x_mean), (1e6 .* block.twtt(ind_bed_flat(pk.ind_x_mean))), 'g--', 'linewidth', 2);
                            end
                        end
                end
                if get(ref_check, 'value')
                    switch ref_start_end
                        case 'start'
                            for ii = 1:pk_ref.num_layer
                                try
                                    plot(block.dist_lin([1 block.ind_overlap(1)]), ...
                                        (1e6 .* block.twtt(repmat(round(ind_y_flat(round(pk_ref.layer(ii).ind_y_smooth(find(~isnan(pk_ref.layer(ii).ind_y_smooth), 1, 'last')) - pk.ind_trim_start + 1), ...
                                        (pk_ref.num_trace - pk_ref.ind_overlap(2) + 1))), 1, 2))), 'y:', 'linewidth', 2);
                                catch
                                    continue
                                end
                            end
                        case 'end'
                            for ii = 1:pk_ref.num_layer
                                try
                                    plot(block.dist_lin([block.ind_overlap(2) end]), ...
                                        (1e6 .* block.twtt(repmat(round(ind_y_flat(round(pk_ref.layer(ii).ind_y_smooth(find(~isnan(pk_ref.layer(ii).ind_y_smooth), 1, 'first')) - pk.ind_trim_start + 1), ...
                                        (pk_ref.num_trace - pk_ref.ind_overlap(2) + 1))), 1, 2))), 'y:', 'linewidth', 2);
                                catch
                                    continue
                                end
                            end
                    end
                end
                if get(flat_check, 'value')
                    for ii = 1:pk.num_layer
                        if ~isnan(pk.layer(ii).ind_y(ind_x_pk))
                            plot(block.dist_lin([1 block.num_trace]), (1e6 .* block.twtt(repmat(round(ind_y_flat(round(pk.layer(ii).ind_y(ind_x_pk)), ind_x_pk)), 1, 2))), 'w:', 'linewidth', 2)
                        else
                            plot(block.dist_lin([1 block.num_trace]), (1e6 .* block.twtt(repmat(round(ind_y_flat(round(pk.layer(ii).ind_y(find(~isnan(pk.layer(ii).ind_y), 1, 'first'))), ...
                                 find(~isnan(pk.layer(ii).ind_y), 1, 'first'))), 1, 2))), 'w:', 'linewidth', 2)
                        end
                    end
                end
                if get(pk_check, 'value')
                    for ii = 1:pk.num_layer
                        tmp1 = plot(block.dist_lin(pk.ind_x_mean(~isnan(pk.layer(ii).ind_y_flat_mean))), (1e6 .* block.twtt(round(pk.layer(ii).ind_y_flat_mean(~isnan(pk.layer(ii).ind_y_flat_mean))))), ...
                                    'r.', 'markersize', 12);
                        if (ii == curr_layer)
                            set(tmp1, 'markersize', 24)
                        end
                    end
                end
                if get(smooth_check, 'value')
                    for ii = 1:pk.num_layer
                        tmp1 = plot(block.dist_lin(ind_decim(~isnan(pk.layer(ii).ind_y_flat_smooth(ind_decim)))), ...
                                    (1e6 .* block.twtt(round(pk.layer(ii).ind_y_flat_smooth(ind_decim(~isnan(pk.layer(ii).ind_y_flat_smooth(ind_decim))))))), 'g.', 'markersize', 12);
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
            tmp1(1, :)      = interp1(block.dist_lin(ind_decim), 1:length(ind_decim), [dist_min dist_max], 'nearest', 'extrap');
            tmp1(2, :)      = interp1(block.twtt, 1:num_sample_trim, [twtt_min twtt_max], 'nearest', 'extrap');
            switch disp_type
                case 'amp.'
                    tmp1    = amp_mean(tmp1(2, 1):tmp1(2, 2), tmp1(1, 1):tmp1(1, 2));
                case 'phase'
                    tmp1    = block.phase_diff_filt(tmp1(2, 1):tmp1(2, 2), tmp1(1, 1):decim:tmp1(1, 2));
                case 'ARESP'
                    tmp1    = block.slope_aresp(tmp1(2, 1):tmp1(2, 2), tmp1(1, 1):decim:tmp1(1, 2));
                case 'spec.'
                    tmp1    = block.specularity(tmp1(2, 1):tmp1(2, 2), tmp1(1, 1):decim:tmp1(1, 2));
                case 'flat'
                    switch flat_switch
                        case 'full'
                            tmp1 = amp_flat(tmp1(2, 1):tmp1(2, 2), tmp1(1, 1):tmp1(1, 2));
                        case 'mean'
                            tmp1(1, :) = interp1(block.dist_lin(pk.ind_x_mean), 1:num_mean, [dist_min dist_max], 'nearest', 'extrap');
                            tmp1 = amp_flat_mean(tmp1(2, 1):tmp1(2, 2), tmp1(1, 1):tmp1(1, 2));
                    end
            end
            tmp2            = NaN(1, 2);
            [tmp2(1), tmp2(2)] ...
                            = deal(nanmean(tmp1(~isinf(tmp1))), nanstd(tmp1(~isinf(tmp1))));
            tmp1            = zeros(1, 2);
            switch disp_type
                case {'amp.' 'flat'}
                    if ((tmp2(1) - (2 * tmp2(2))) < db_min_ref)
                        tmp1(1) = db_min_ref;
                    else
                        tmp1(1) = tmp2(1) - (2 * tmp2(2));
                    end
                    if ((tmp2(1) + (2 * tmp2(2))) > db_max_ref)
                        tmp1(2) = db_max_ref;
                    else
                        tmp1(2) = tmp2(1) + (2 * tmp2(2));
                    end
                    [db_min, db_max] ...
                            = deal(tmp1(1), tmp1(2));
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
                case 'phase'
                    if ((tmp2(1) - (2 * tmp2(2))) < phase_diff_min_ref)
                        tmp1(1) = phase_diff_min_ref;
                    else
                        tmp1(1) = tmp2(1) - (2 * tmp2(2));
                    end
                    if ((tmp2(1) + (2 * tmp2(2))) > phase_diff_max_ref)
                        tmp1(2) = phase_diff_max_ref;
                    else
                        tmp1(2) = tmp2(1) + (2 * tmp2(2));
                    end
                    [phase_diff_min, phase_diff_max] ...
                            = deal(tmp1(1), tmp1(2));
                    if (phase_diff_min < get(cb_min_slide, 'min'))
                        set(cb_min_slide, 'value', get(cb_min_slide, 'min'))
                    else
                        set(cb_min_slide, 'value', phase_diff_min)
                    end
                    if (phase_diff_max > get(cb_max_slide, 'max'))
                        set(cb_max_slide, 'value', get(cb_max_slide, 'max'))
                    else
                        set(cb_max_slide, 'value', phase_diff_max)
                    end
                    set(cb_min_edit, 'string', sprintf('%3.2f', phase_diff_min))
                    set(cb_max_edit, 'string', sprintf('%3.2f', phase_diff_max))
                    caxis([phase_diff_min phase_diff_max])
                case 'ARESP'
                    if ((tmp2(1) - (2 * tmp2(2))) < aresp_min_ref)
                        tmp1(1) = aresp_min_ref;
                    else
                        tmp1(1) = tmp2(1) - (2 * tmp2(2));
                    end
                    if ((tmp2(1) + (2 * tmp2(2))) > aresp_max_ref)
                        tmp1(2) = aresp_max_ref;
                    else
                        tmp1(2) = tmp2(1) + (2 * tmp2(2));
                    end
                    [aresp_min, aresp_max] ...
                            = deal(tmp1(1), tmp1(2));
                    if (aresp_min < get(cb_min_slide, 'min'))
                        set(cb_min_slide, 'value', get(cb_min_slide, 'min'))
                    else
                        set(cb_min_slide, 'value', aresp_min)
                    end
                    if (aresp_max > get(cb_max_slide, 'max'))
                        set(cb_max_slide, 'value', get(cb_max_slide, 'max'))
                    else
                        set(cb_max_slide, 'value', aresp_max)
                    end
                    set(cb_min_edit, 'string', sprintf('%3.2f', aresp_min))
                    set(cb_max_edit, 'string', sprintf('%3.2f', aresp_max))
                    caxis([aresp_min aresp_max])
                case 'spec.'
                    if ((tmp2(1) - (2 * tmp2(2))) < spec_min_ref)
                        tmp1(1) = spec_min_ref;
                    else
                        tmp1(1) = tmp2(1) - (2 * tmp2(2));
                    end
                    if ((tmp2(1) + (2 * tmp2(2))) > spec_max_ref)
                        tmp1(2) = spec_max_ref;
                    else
                        tmp1(2) = tmp2(1) + (2 * tmp2(2));
                    end
                    [spec_min, spec_max] ...
                            = deal(tmp1(1), tmp1(2));
                    if (spec_min < get(cb_min_slide, 'min'))
                        set(cb_min_slide, 'value', get(cb_min_slide, 'min'))
                    else
                        set(cb_min_slide, 'value', spec_min)
                    end
                    if (spec_max > get(cb_max_slide, 'max'))
                        set(cb_max_slide, 'value', get(cb_max_slide, 'max'))
                    else
                        set(cb_max_slide, 'value', spec_max)
                    end
                    set(cb_min_edit, 'string', sprintf('%3.2f', spec_min))
                    set(cb_max_edit, 'string', sprintf('%3.2f', spec_max))
                    caxis([spec_min spec_max])
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
                if get(man_check, 'value')
                    set(man_check, 'value', 0)
                else
                    set(man_check, 'value', 1)
                end
                show_man
            case '4'
                if get(flat_check, 'value')
                    set(flat_check, 'value', 0)
                else
                    set(flat_check, 'value', 1)
                end
                show_flat
            case '5'
                if get(mean_check, 'value')
                    set(mean_check, 'value', 0)
                else
                    set(mean_check, 'value', 1)
                end
                show_mean
            case '6'
                if get(pk_check, 'value')
                    set(pk_check, 'value', 0)
                else
                    set(pk_check, 'value', 1)
                end
                show_pk
            case '7'
                if get(smooth_check, 'value')
                    set(smooth_check, 'value', 0)
                else
                    set(smooth_check, 'value', 1)
                end
                show_smooth
            case '8'
                if get(aresp_check, 'value')
                    set(aresp_check, 'value', 0)
                else
                    set(aresp_check, 'value', 1)
                end
                show_aresp
            case '9'
                pk_man
            case 'a'
                adj_layer
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
                del_layer
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
                pk_keep_phase
            case 'l'
                load_data
            case 'm'
                merge_layer
            case 'n'
                pk_next
            case 'o'
                focus_layer
            case 'p'
                pk_flat
            case 'q'
                pop_fig
            case 'r'
                load_ref
            case 's'
                pk_smooth
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
                choose_pk
            case 'downarrow'
                tmp1        = twtt_max - twtt_min;
                tmp2        = twtt_max + (0.25 * tmp1);
                if (tmp2 > twtt_max_ref)
                    twtt_max= twtt_max_ref;
                else
                    twtt_max= tmp2;
                end
                twtt_min    = twtt_max - tmp1;
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
                if (tmp2 > dist_max_ref)
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
                tmp1        = twtt_max - twtt_min;
                tmp2        = twtt_min - (0.25 * tmp1);
                if (tmp2 < twtt_min_ref)
                    twtt_min= twtt_min_ref;
                else
                    twtt_min= tmp2;
                end
                twtt_max= twtt_min + tmp1;
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
            case 'space'
                if (get(disp_group, 'selectedobject') == disp_check(1))
                    if flat_done
                        set(disp_group, 'selectedobject', disp_check(5))
                        disp_type = 'flat';
                        plot_flat
                    elseif phase_avail
                        set(disp_group, 'selectedobject', disp_check(2))
                        disp_type = 'phase';
                        plot_phase_diff
                    elseif aresp_avail
                        set(disp_group, 'selectedobject', disp_check(3))
                        disp_type = 'ARESP';
                        plot_aresp
                    elseif spec_avail
                        set(disp_group, 'selectedobject', disp_check(4))
                        disp_type = 'spec.';
                        plot_spec                        
                    end
                else
                    set(disp_group, 'selectedobject', disp_check(1))
                    disp_type = 'amp.';
                    plot_db
                end
        end
    end

%% Test something

    function misctest(source, eventdata)
        pk.ind_trim_start
    end

%%
end