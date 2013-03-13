% PHASEINTERP_WRAPPER Run phaseinterp on a set of transects.
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
dir_data_old                = '2012_p3/data_old/';
dir_block_new               = '2012_p3/block/';
sz_filt                     = 512;
wavenum_max                 = 0.1;
freq                        = 195e6;
lat_std                     = 70;
decim                       = 10;

% check each directory is indeed a directory
dir_new                     = dir(dir_block_new);
ind_dir                     = false(1, length(dir_new));
for ii = 1:length(dir_new)
    if dir_new(ii).isdir
        ind_dir(ii)         = true;
    end
end

% extract names and ditch the ones we don't need
dir_new                     = dir_new(ind_dir);
dir_new                     = {dir_new.name};
dir_new                     = dir_new(~strcmp(dir_new, '.'));
dir_new                     = dir_new(~strcmp(dir_new, '..'));
num_trans                   = length(dir_new);

disp(['Filtering and interpolating horizontal phase gradient in all (' num2str(num_trans) ' transects) single-channel blocks in ' dir_data_old ' to focused blocks in ' dir_block_new '...'])

% call phaseinterp
if parallel_check
    parfor ii = 1:num_trans
        disp([dir_new{ii} ' (' num2str(ii) ' / ' num2str(num_trans) ')...'])
        if ~exist([dir_data_old dir_new{ii}], 'dir')
            disp([dir_new{ii} ' does not exist in ' dir_data_old '. Moving on to next transect.'])
            continue
        end
        phaseinterp([dir_data_old dir_new{ii} '/'], sz_filt, wavenum_max, freq, lat_std, decim, [dir_block_new dir_new{ii} '/'], '*')
    end
else
    for ii = 1:num_trans
        disp([dir_new{ii} ' (' num2str(ii) ' / ' num2str(num_trans) ')...'])
        if ~exist([dir_data_old dir_new{ii}], 'dir')
            disp([dir_new{ii} ' does not exist in ' dir_data_old '. Moving on to next transect.'])
            continue
        end
        phaseinterp([dir_data_old dir_new{ii} '/'], sz_filt, wavenum_max, freq, lat_std, decim, [dir_block_new dir_new{ii} '/'], '*')
    end    
end

if parallel_check
    matlabpool close
end

disp(['DONE interpolating phase for all focused transects in ' dir_block_new '.'])