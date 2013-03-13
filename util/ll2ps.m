function [x, y]     = ll2ps(lat, lon, lat_std)
% LL2PS Geographic to polar stereographic coordinate conversion.
%   [X,Y] = LL2PS(LAT,LON) converts latitude LAT and longitude LON in
%   decimal degrees into polar stereographic coordinates in kilometers,
%   using the WGS84 ellipsoid. The standard parallel is set at 70 degrees
%   south.
% 
%   [X,Y] = LL2PS(LAT,LON,LAT_STD) uses LAT_STD in decimal degrees as the
%   standard parallel. North is positive, and South is negative. Using
%   LAT_STD=-71 produces the standard SCAR projection for Antarctica. Using
%   LAT_STD=70 produces an acceptable projection for Greenland, but not one
%   that has been corrected to the commonly used central meridian (-45).
% 
%   Modified from Ben Smith (UW-APL)
% 
% Joe MacGregor (UTIG)
% Last updated: 02/17/12
 
if (nargin < 3);
    lat_std         = -70;
end;
sn                  = sign(lat_std);
lat_std             = abs(lat_std);

DR                  = pi / 180;
CDR                 = 1 / DR;
rad_earth           = 6378137; % WGS84
ecc                 = 0.08199188997903; % WGS84?
ecc2                = 0.00672267002233; % WGS84?
t                   = tan((pi / 4) - (lat_std / (2 * CDR))) / (((1 - (ecc * sin(DR * lat_std))) / (1 + (ecc * sin(DR * lat_std)))) ^ (ecc / 2));
cm                  = cos(DR * lat_std) / sqrt(1 - (ecc2 * (sin(DR * lat_std) ^ 2)));

lat                 = lat .* sn;
t1                  = tan((pi / 4) - (lat / (2 * CDR))) ./ (((1 - (ecc .* sin(DR .* lat))) ./ (1 + (ecc .* sin(DR .* lat)))) .^ (ecc / 2));
rho                 = ((rad_earth * cm) / t) .* t1;
x                   = -real(rho.* sn .* sin(DR .* lon));
y                   = -sn .* real(rho .* cos(DR .* lon));