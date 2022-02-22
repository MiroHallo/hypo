function [Loc] = locWGS84(X,Y,refLat,refLon)

EarthD = 12756.320;   

Loc = zeros(length(X),2);

for sta = 1:length(X)
    Loc(sta,1) = refLat + (Y(sta) / (EarthD*pi/360));
    Loc(sta,2) = refLon + (X(sta) / (EarthD*pi*cosd(refLat)/360));
end

return