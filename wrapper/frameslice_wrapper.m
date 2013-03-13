% FRAMESLICE_WRAPPER Run frameslice on a set of transects.
% 
% Joe MacGregor (UTIG)
% Last updated: 02/01/13

clear

% user defined parameters
dir_in                      = '2012_p3/data_new/';
dir_slice                   = '2012_p3/data_slice_new/';
num_slice                   = 5;

% check each directory is indeed a directory (only directories in dir_in please!)
dir_all                     = dir(dir_in);
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

% make equivalent directories in dir_slice
% for ii = 1:num_dir
%     eval(['mkdir ' dir_slice dir_all{ii}])
% end

disp(['Slicing all (' num2str(num_dir) ') transects in ' dir_in '...'])

% call frameslice
for ii = 1:num_dir
    disp([dir_all{ii} ' (' num2str(ii) ' / ' num2str(num_dir) ')...'])
    try
        frameslice([dir_in dir_all{ii} '/'], '*', num_slice, [dir_slice dir_all{ii} '/'])
    catch %#ok<CTCH>
        continue
    end
end

disp(['DONE slicing all transects in ' dir_in '.'])