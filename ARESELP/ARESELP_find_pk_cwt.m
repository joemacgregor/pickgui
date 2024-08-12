function pk_cwt				= ARESELP_find_pk_cwt(amp, wavelet_scale, wavelet_shp, ind_surf, ind_bed, pk_prom_min, twtt_sep_min_ind, num_pk_max)
% ARESELP_FIND_PK_CWT Identify peaks in the wavelet transform of a single radar trace.
% 
%	Only use ARESELP_FIND_PK_CWT as called by ARESELP.
% 
%	ARESELP_FIND_PK_CWT is based on code archived at:
%	https://github.com/xiongsiting/ARESELP/
% 
% Joe MacGregor (NASA)
% Last updated: 12 August 2024

% continuous wavelet transform coefficients, column vector becomes a row
cwt_coef					= cwt(amp, wavelet_scale, wavelet_shp);

if isnan(ind_surf)
	ind_surf				= 0;
end
if isnan(ind_bed)
	ind_bed					= length(amp);
end

pk_cwt						= zeros(length(amp), 1);

% loop through each wavelet scale
for ii = 1:size(cwt_coef, 1)
	
	% find peaks in wavelet transform
	[~, ind_pk]				= findpeaks(cwt_coef(ii, (ind_surf + twtt_sep_min_ind):(ind_bed - twtt_sep_min_ind)), 'MinPeakProminence', pk_prom_min, 'MinPeakDistance', twtt_sep_min_ind, 'NPeaks', num_pk_max);
	
	% correct for surface start and buffer
	ind_pk					= ind_pk + ind_surf + twtt_sep_min_ind - 1;
	
	% remove peaks above surface and below bottom
	ind_pk((ind_pk < (ind_surf + twtt_sep_min_ind)) | (ind_pk > (ind_bed - twtt_sep_min_ind))) ...
							= [];
	
	% move along if no peaks
	if isempty(ind_pk)
		continue
	end
	
	% adjust peak locations based on wavelet scale
	half_scale				= floor(wavelet_scale(ii) / 2);
	for jj = 1:length(ind_pk)
		[~, ind_max]		= max(amp((ind_pk(jj) - half_scale):(ind_pk(jj) + half_scale)));
		if ((ind_pk(jj) - half_scale - 1  + ind_max) < ind_bed)
			ind_pk(jj)		= ind_pk(jj) - half_scale - 1  + ind_max;
		end
	end
	
	% preserve CWT coefficients if peaks good	
	pk_cwt(ind_pk)			= pk_cwt(ind_pk) + cwt_coef(ii, ind_pk)';
end