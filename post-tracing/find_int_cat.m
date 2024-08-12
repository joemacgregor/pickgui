% FIND_INT_CAT Determine intersections of radar segments.
% 
%   FIND_INT_CAT determines the positions of intersections between each
%   segment and all other segments.
% 
% Joe MacGregor (NASA)
% Last updated: 16 July 2024

clear

dir_mat						= '/Users/jamacgre/OneDrive - NASA/research/matlab/pickgui_v2/mat/';
do_int                      = true;
do_save						= true;
dist_density                = 10e3; % maximum intersection density, km
angle_threshold             = 5; % minimum intersection angle, deg
decim                       = 25;
dist_diff_max				= 5; % multiple of dist_diff median above which to ignore intersection as a big gap within radargram
plotting                    = false;

% get x/y for traced segments
load([dir_mat 'xyz_all.mat'], 'num_campaign')
pk_cat						= load([dir_mat 'pk_cat.mat'], 'dist', 'file_pk', 'ind_campaign', 'ind_pk_merge', 'num_file_pk', 'x', 'y');
[dist, file_pk, ind_campaign, ind_pk_merge, num_file_pk, x, y] ...
							= deal(pk_cat.dist, pk_cat.file_pk, pk_cat.ind_campaign, pk_cat.ind_pk_merge(:, 1:2), pk_cat.num_file_pk, pk_cat.x, pk_cat.y);
dist_diff					= cell(1, pk_cat.num_file_pk);
for ii = 1:pk_cat.num_file_pk
	file_pk{ii}				= file_pk{ii}(6:(end - 7));
	dist_diff{ii}			= diff(dist{ii});
end

%%

% BedMachine v5
BM5                         = struct;
BM5.x                       = double(ncread('/Users/jamacgre/OneDrive - NASA/research/data/greenland/BedMachine/BedMachineGreenland-v5.nc', 'x'))'; % projected x, km
BM5.y                       = double(flipud(ncread('/Users/jamacgre/OneDrive - NASA/research/data/greenland/BedMachine/BedMachineGreenland-v5.nc', 'y'))); % projected y, km
BM5.mask_gris               = double(rot90(ncread('/Users/jamacgre/OneDrive - NASA/research/data/greenland/BedMachine/BedMachineGreenland-v5.nc', 'mask'))); % ice mask

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

% Mouginot 2019 ice-sheet drainage polygon basins, with peripheral ice masses manually removed in QGIS, then rasterized to BM5 grid, took me forever to figure out but got it with SAGA in QGIS
M19_ice_mask				= readgeoraster('/Users/jamacgre/OneDrive - NASA/data/GreenValley/Greenland_Basins_PS_v1.4.2_simple_raster.tif');
M19_ice_mask				= flipud(logical(M19_ice_mask));

% now BM5 mask 2=peripheral ice massess and 3=ice sheet (previously 2)
BM5.mask_combo_plot			= BM5.mask_combo;
BM5.mask_combo((BM5.mask_combo == 2) & M19_ice_mask) ...
							= 3;
BM5.mask_gris(BM5.mask_combo == 2) ...
							= false;

clear M19_ice_mask

% X/Y limits for standard map display, m
[x_min, x_max, y_min, y_max]= deal(-632e3, 846e3, -3344e3, -670e3);

% reference 1-km X/Y grid, m
[x_grd, y_grd]              = meshgrid(x_min:1e3:x_max, (y_min:1e3:y_max)');
[num_y_grd, num_x_grd]		= size(x_grd);

mask_ice_grd				= false(size(x_grd));
mask_ice_grd(interp2(BM5.x, BM5.y, double(BM5.mask_combo), x_grd, y_grd) == 3) ...
							= true;
% clear BM5

%% simplify identification of potentially matching segments

ind_empty_segment			= [];
ind_grd						= cell(1, num_file_pk);

for ii = 1:num_file_pk
	ind_grd{ii}				= unique([interp1(y_grd(:, 1)', 1:num_y_grd, y{ii}, 'nearest', 'extrap')' interp1(x_grd(1, :), 1:num_x_grd, x{ii}, 'nearest', 'extrap')'], 'rows'); % subscript good x/ys
	ind_grd{ii}				= sub2ind([num_y_grd num_x_grd], ind_grd{ii}(:, 1), ind_grd{ii}(:, 2)); % subscript to indices
	ind_grd{ii}				= ind_grd{ii}(mask_ice_grd(ind_grd{ii})); % get rid of non-ice values	
end

%% Calculate intersections between each segment all others (all others in that campaign and all other campaigns)

if do_int

    disp('Finding all other intersections...')
    
    num_int                 = zeros(num_file_pk); % number of intersections
    int_all                 = [];
    
	for ii = 1:(num_file_pk - 1)
		
        disp([file_pk{ii} ' (' num2str(ii) ' / ' num2str(num_file_pk - 1) ')...'])
        pause(0.25)

		ind_tmp1            = 1:decim:length(x{ii});
		ind_tmp1(end)       = length(x{ii});
        
		for jj = (ii + 1):num_file_pk
            
			if isempty(intersect(ind_grd{ii}, ind_grd{jj}))
				continue
			end
			
            disp(['...compared with ' file_pk{jj} '...'])
            
            int_all_sz_old  = size(int_all, 1);
        	ind_tmp2		= 1:decim:length(x{jj});
        	ind_tmp2(end)	= length(x{jj});
			
			[~, ~, ind_int1, ind_int2]  = intersecti(x{ii}(ind_tmp1), y{ii}(ind_tmp1), x{jj}(ind_tmp2), y{jj}(ind_tmp2));
			
			if ~isempty(ind_int1)
				
				for kk = 1:length(ind_int1)
					
					if (ind_int1(kk) == 1)
						ind_tmp3                = 1:ind_tmp1(2);
					elseif (ind_int1(kk) == length(ind_tmp1))
						ind_tmp3                = ind_tmp1(end - 1):ind_tmp1(end);
					else
						ind_tmp3                = ind_tmp1(ind_int1(kk) - 1):ind_tmp1(ind_int1(kk) + 1);
					end
					
					if (ind_int2(kk) == 1)
						ind_tmp4                = 1:ind_tmp2(2);
					elseif (ind_int2(kk) == length(ind_tmp2))
						ind_tmp4                = ind_tmp2(end - 1):ind_tmp2(end);
					else
						ind_tmp4                = ind_tmp2(ind_int2(kk) - 1):ind_tmp2(ind_int2(kk) + 1);
					end
					
					% intersection of two vectors using linear interpolation
					[~, ~, ind_int3, ind_int4]  = intersecti(x{ii}(ind_tmp3), y{ii}(ind_tmp3), x{jj}(ind_tmp4), y{jj}(ind_tmp4));
					
					if ~isempty(ind_int3) % only when there is an intersection add a row for ind_match segment
						[ind_int3, ind_int4]    = deal((ind_int3 + ind_tmp3(1) - 1), (ind_int4 + ind_tmp4(1) - 1));
						int_all					= [int_all; repmat(ii, length(ind_int3), 1) ind_int3 NaN(length(ind_int3), 1) repmat(jj, length(ind_int3), 1) ind_int4 NaN(length(ind_int3), 1)]; %#ok<AGROW>
					end
				end
			end
			
            ind2remove      = [];
			
            % remove intersections due to segment jumps or small intersection angles, or high intersection density
            for kk = (int_all_sz_old + 1):size(int_all, 1)
				
                ind_tmp5    = (int_all(kk, 2) - 1):(int_all(kk, 2) + 1);
                ind_tmp5    = ind_tmp5((ind_tmp5 > 0) & (ind_tmp5 <= length(dist_diff{ii})));
                if any(dist_diff{ii}(ind_tmp5) > (dist_diff_max * median(dist_diff{ii})))
                    ind2remove ...
                            = [ind2remove kk]; %#ok<AGROW>
					continue
                end
				
                ind_tmp6    = (int_all(kk, 5) - 1):(int_all(kk, 5) + 1);
                ind_tmp6    = ind_tmp6((ind_tmp6 > 0) & (ind_tmp6 <= length(dist_diff{jj})));
                if any(dist_diff{jj}(ind_tmp6) > (dist_diff_max * median(dist_diff{jj})))
                    ind2remove ...
                            = [ind2remove kk]; %#ok<AGROW>
					continue
                end
				
				if (kk < size(int_all, 1))
					if (abs(diff([dist{ii}(int_all(kk, 2)) dist{ii}(int_all((kk + 1), 2))])) < dist_density) % ensure density of consecutive intersections is not too great
            			ind2remove ...
							= [ind2remove kk]; %#ok<AGROW>
						continue
					end
					if (abs(diff([dist{jj}(int_all(kk, 5)) dist{jj}(int_all((kk + 1), 5))])) < dist_density)
            			ind2remove ...
							= [ind2remove kk]; %#ok<AGROW>
						continue
					end
				end
				
                int_all(kk, [3 6]) ...
                            = [atand(diff(y{ii}(ind_tmp5([1 end]))) ./ diff(x{ii}(ind_tmp5([1 end])))) atand(diff(y{jj}(ind_tmp6([1 end]))) ./ diff(x{jj}(ind_tmp6([1 end]))))]; %#ok<SAGROW>
                % segment intersection angle must be greater than angle_threshold
				if (abs(diff(int_all(kk, [3 6]))) < angle_threshold)
                    ind2remove ...
                            = [ind2remove kk]; %#ok<AGROW>
					continue
				end
            end
			
            if ~isempty(ind2remove)
                int_all     = int_all(setdiff(1:size(int_all, 1), ind2remove), :); % the segment intersection where data actually exist
            end
			
			num_int(ii, jj) = length((int_all_sz_old + 1):size(int_all, 1));
		end
	end

%%	
	if do_save
    	save([dir_mat 'int_all_cat'], 'int_all', 'num_int')
    	disp(['Done calculating all intersections and saved in ' dir_mat '.'])
	end
	
else
    load([dir_mat 'int_all_cat'], 'int_all', 'num_int') %#ok<UNRCH>
    disp(['Loaded all intersections from ' dir_mat '.'])
end

%%
if plotting

%% all segments and intersections

    figure
    hold on
	imagesc(BM5.x, BM5.y, BM5.mask_combo_plot)
	colormap([([135 206 235] ./ 255); ([159 89 39] ./ 255); 0.9 0.9 0.9])
    % colors                  = jet(num_campaign);
	[pt, put]				= deal(gobjects(0));
	for ii = 1:num_file_pk
		line(x{ii}(1:decim:end), y{ii}(1:decim:end), 'LineStyle', '-', 'LineWidth', 1, 'Color', 'k')%colors(ind_campaign(ii), :), 'Tag', file_pk{ii})
		try %#ok<TRYNC>
			% line(x{ii}(int_all((int_all(:, 1) == ii), 2)), y{ii}(int_all((int_all(:, 1) == ii), 2)), 'Color', 'k', 'Marker', 's', 'MarkerSize', 10, 'MarkerFaceColor', 'm', 'LineWidth', 1, 'LineStyle', 'none', 'Tag', num2str(ii))
			line(x{ii}(int_all((int_all(:, 4) == ii), 5)), y{ii}(int_all((int_all(:, 4) == ii), 5)), 'Color', 'k', 'Marker', '^', 'MarkerSize', 6, 'MarkerFaceColor', 'm', 'LineWidth', 1, 'LineStyle', 'none', 'Tag', num2str(ii))
		end
	end
    set(gca, 'FontSize', 20)
    xlabel('Polar stereographic X (m)')
    ylabel('Polar stereographic Y (m)')
    axis xy image
    grid on
    box on
    
%%
end

% %% NO LONGER NEEDED BECAUSE I FIXED THIS IN RADFRAMEPROC JUNE 2024 AND REDID IT ALL
% 	% adjust to pk reference frame (because non-unique times/positions were removed)
% 	int_all_pk				= int_all;
% 	for ii = 1:size(int_all_pk, 1)
% 		ind_tmp1			= find(ind_trace_ref{int_all_pk(ii, 1)} > 0);
% 		ind_tmp2			= find(ind_trace_ref{int_all_pk(ii, 4)} > 0);
% 		int_all_pk(ii, 2)	= ind_tmp1(int_all_pk(ii, 2));
% 		int_all_pk(ii, 5)	= ind_tmp2(int_all_pk(ii, 5));
% 	end