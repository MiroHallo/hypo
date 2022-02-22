function [Loc] = locX(Lat,Lon,refLat,refLon)

EarthD = 12756.320;   

Loc = zeros(length(Lon),2);

for sta = 1:length(Lon)
    Loc(sta,1) = (Lon(sta)-refLon)*EarthD*pi*cosd(refLat)/360;
    Loc(sta,2) = (Lat(sta)-refLat)*EarthD*pi/360;
end

return