function [x, y]             = ll2ps(lat, lon, lat_std)
% LL2PS Convert geographic coordinates to polar stereographic coordinates.
%   
%   [X,Y] = LL2PS(LAT,LON) converts latitude LAT and longitude LON in
%   decimal degrees into polar stereographic coordinates in meters using
%   the WGS84 ellipsoid. The default standard parallel is set at 71S, which
%   is the standard SCAR projection for Antarctica.
%   
%   [X,Y] = LL2PS(LAT,LON,LAT_STD) uses LAT_STD in decimal degrees as the
%   standard parallel. North is positive, and south is negative. LAT_STD =
%   70 produces an acceptable projection for Greenland, but not one that
%   has been corrected to the commonly used central meridian (45W or -45).
%   
%   Modified from Ben Smith (UW-APL)
%   
% Joe MacGregor (UTIG)
% Last updated: 03/20/13

if ~any(nargin == [2 3])
    error('ll2ps:nargin', 'Incorrect number of inputs (must be 2 or 3).')
end
if ~isnumeric(lat)
    error('ll2ps:latnum', 'LAT is not numeric.')
end
if any(diff(sign(lat)))
    error('ll2ps:latsign', 'Values of LAT must all have same sign.')
end
if ~isnumeric(lon)
    error('ll2ps:lonnum', 'LON is not numeric.')
end
if (length(lat) ~= length(lon))
    error('ll2ps:latlonlength', 'Lengths of LAT and LON must be the same.')
end
if (nargin < 3)
    lat_std                 = -71;
elseif (~isnumeric(lat_std) || (length(lat_std) ~= 1))
    error('ll2ps:latstd', 'LAT_STD is not a numeric scalar.')
elseif (sign(lat_std) ~= sign(lat(find(~isnan(lat), 1))))
    error('ll2ps:latsigncomp', 'LAT and LAT_STD do not have the same sign.')
end
if (nargout ~= 2)
    error('ll2ps:nargout', 'Incorrect number of outputs (must be 2).')
end

sn                          = sign(lat_std);
lat_std                     = abs(lat_std);
dr                          = pi / 180;
cdr                         = 1 / dr;
rad_earth                   = 6378137; % WGS84
ecc                         = 0.08199188997903; % WGS84?
ecc2                        = 0.00672267002233; % WGS84?
t                           = tan((pi / 4) - (lat_std / (2 * cdr))) / (((1 - (ecc * sin(dr * lat_std))) / (1 + (ecc * sin(dr * lat_std)))) ^ (ecc / 2));
cm                          = cos(dr * lat_std) / sqrt(1 - (ecc2 * (sin(dr * lat_std) ^ 2)));
lat                         = lat .* sn;
t1                          = tan((pi / 4) - (lat / (2 * cdr))) ./ (((1 - (ecc .* sin(dr .* lat))) ./ (1 + (ecc .* sin(dr .* lat)))) .^ (ecc / 2));
rho                         = ((rad_earth * cm) / t) .* t1;
x                           = -real(rho.* sn .* sin(dr .* lon));
y                           = -sn .* real(rho .* cos(dr .* lon));
[x, y]                      = deal((1e-3 .* x), (1e-3 .* y));