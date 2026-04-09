function [Loc] = locX(Lat, Lon, refLat, refLon)
% LOCX Converts WGS84 Longitude and Latitude into Easting and Northing.
% 
% Author: Miroslav HALLO
% Charles University in Prague, 2017/06
% 
% Copyright (C) 2017 Miroslav Hallo
% This program is published under the GNU General Public License (GNU GPL).
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
% Lat - WGS84 Latitude [deg]
% Lon - WGS84 Longitude [deg]
% refLat - Reference Latitude [deg]
% refLon - Reference Longitude [deg]
% 
% OUTPUT:
% Loc(:,1) - Easting [km]
% Loc(:,2) - Northing [km]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Earth diameter [km]
EarthD = 12756.320;   

Loc = zeros(length(Lon),2);

for sta = 1:length(Lon)
    Loc(sta,1) = (Lon(sta)-refLon)*EarthD*pi*cosd(refLat)/360;
    Loc(sta,2) = (Lat(sta)-refLat)*EarthD*pi/360;
end

end