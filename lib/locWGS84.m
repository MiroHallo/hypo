function [Loc] = locWGS84(X, Y, refLat, refLon)
% LOCWGS84 Converts Easting and Northing into WGS84 Longitude and Latitude.
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
% X - Easting [km]
% Y - Northing [km]
% refLat - Reference Latitude [deg]
% refLon - Reference Longitude [deg]
% 
% OUTPUT:
% Loc(:,1) - WGS84 Latitude [deg]
% Loc(:,2) - WGS84 Longitude [deg]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Earth diameter [km]
EarthD = 12756.320;

Loc = zeros(length(X),2);

for sta = 1:length(X)
    Loc(sta,1) = refLat + (Y(sta) / (EarthD*pi/360));
    Loc(sta,2) = refLon + (X(sta) / (EarthD*pi*cosd(refLat)/360));
end

end