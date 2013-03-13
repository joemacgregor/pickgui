% ARESP_WRAPPER Perform ARESP on a set of transects.
% 
% Joe MacGregor (UTIG)
% Last updated: 03/13/13

clear

if license('test', 'distrib_computing_toolbox')
    if matlabpool('size')
        matlabpool close
    end
    parallel_check          = true;
    matlabpool
else
    parallel_check          = false;
end

% user defined parameters
dir_block                   = '2011_p3/block/';
decim_x                     = [10 10];
decim_z                     = 5;
z_smooth                    = [6 36];
num_overlap                 = [2 2];
area_max                    = 3; % maximum ratio of object area to (ind_sep * z_smooth)
ellipse_min                 = 2; % minimum object ellipticity (MajorAxisLength/MinorAxisLength)
angle_max                   = 20; % maximum object orientation
dist_max                    = 20; % maximum number of indices
num_obj_max                 = 80; % maximum number of objects

% check each directory is indeed a directory (only directories in dir_block please!)
dir_all                     = dir(dir_block);
ind_dir                     = false(1, length(dir_all));
for ii = 1:length(dir_all)
    if dir_all(ii).isdir
        ind_dir(ii)         = true;
    end
end

% extract names and ditch the ones we don't need
dir_all                     = dir_all(ind_dir);
dir_all                     = {dir_all.name}; % just need names now, ignore the rest
dir_all                     = dir_all(~strcmp(dir_all, '.'));
dir_all                     = dir_all(~strcmp(dir_all, '..'));
num_trans                   = length(dir_all);

disp(['Performing ARESP on all (' num2str(num_trans) ') transects in ' dir_block '...'])

% call aresp
for ii = 10%1:num_trans
    disp([dir_all{ii} ' (' num2str(ii) ' / ' num2str(num_trans) ')...'])
    aresp([dir_block dir_all{ii} '/'], '*', decim_x, decim_z, z_smooth, num_overlap, parallel_check, area_max, ellipse_min, angle_max, dist_max, num_obj_max)
end

if parallel_check
    matlabpool close
end

disp(['DONE with ARESP for all focused transects in ' dir_block '.'])