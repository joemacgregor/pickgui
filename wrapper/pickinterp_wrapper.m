% PICKINTERP_WRAPPER Run pickinterp on a set of transects.
% 
% Joe MacGregor (UTIG)
% Last updated: 11/19/12

clear

% user defined parameters
dir_merge                   = '2011_p3/merge/';
dir_data                    = '2010_p3/data/';
dir_pk                      = '2010_p3/pk_orig/';

% only merged pk files should be in dir_merge
file_all                    = dir([dir_merge, '*.mat']);
file_all                    = {file_all.name};
num_file                    = length(file_all);

% make equivalent directories in dir_pk
for ii = 1:num_dir
    eval(['mkdir ' dir_pk file_all{ii}(1:11)])
end

disp(['Interpolating all (' num2str(num_file) ') merged picks in ' dir_merge ' to original CReSIS data in ' dir_data '...'])

% call pickinterp
for ii = 1:num_file
    disp([file_all{ii} ' (' num2str(ii) ' / ' num2str(num_file) ')...'])
    pickinterp(dir_merge, file_all{ii}, [dir_data file_all{ii}(1:11) '/'], '*', [dir_pk file_all{ii}(1:11) '/'])
end

disp(['DONE interpolating all picks in ' dir_merge '.'])