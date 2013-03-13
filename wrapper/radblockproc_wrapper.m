% ****** RADBLOCKPROC_WRAPPER.M ******
% 
% Run radblockproc on a set of directories.
% 
% Joe MacGregor (UTIG)
% Last updated: 02/01/13

clear

% user defined parameters
dir_data                    = '2012_p3/data_slice_new/';
dir_block                   = '2012_p3/block_new/';
% dir_pk                      = '2012_p3/pk/';
num_file_block              = 6;
num_overlap                 = 1;
lat_std                     = 70;

% check each directory is indeed a directory (only directories in dir_data please!)
dir_all                     = dir(dir_data);
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
num_dir                     = length(dir_all);

% make equivalent directories in dir_block and dir_pk
for ii = 1:num_dir
%     eval(['!mkdir ' dir_block dir_all{ii}]);
%     eval(['!mkdir ' dir_pk dir_all{ii}]);
end

disp(['Processing all (' num2str(num_dir) ') transects in ' dir_data '...'])

% call radblockproc
for ii = 8:num_dir
    disp([dir_all{ii} ' (' num2str(ii) ' / ' num2str(num_dir) ')...'])
    try
        radblockproc([dir_data dir_all{ii} '/'], '*', [dir_all{ii} '_block' num2str(num_file_block)], num_file_block, num_overlap, lat_std, [dir_block dir_all{ii} '/'])
    catch %#ok<CTCH>
        continue
    end
end

disp(['DONE processing all transects in ' dir_data '.'])