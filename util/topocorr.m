function [data_corr, ind_corr]  = topocorr(data, depth, elev_surf, do_strip)
% TOPOCORR Correct radargram for surface topography.
%   [DATA_CORR,IND_CORR] = TOPOCORR(DATA,DEPTH,ELEV_SURF,DO_STRIP) corrects the
%   radar data matrix DATA for the spatial variation in surface topography
%   described by vector ELEV_SURF, given the elevation vector ELEV. It
%   returns the corrected data matrix DATA_CORR and the vector of indices
%   IND_CORR by which each trace was shifted by.
%   
% Joe MacGregor (NASA)
% Last updated: 2 November 2024

[m, n]                          = size(data);
ind_corr                        = round((elev_surf - min(elev_surf)) ./ (depth(end) - depth(end - 1))); % index correction
max_corr                        = max(ind_corr);
data_corr                       = NaN((m + max_corr), n, 'single'); % initializes topo-corrected data to a size large enough for largest correction
for ii = 1:n
    data_corr((1 - ind_corr(ii) + max_corr):(m - ind_corr(ii) + max_corr), ii) ...
                                = data(:, ii); % topo corrected data
end
if do_strip
	data_corr                   = data_corr(1:m, :); % strips away bottom part of data
end