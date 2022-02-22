%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Convert Latitude Longitude to xyz local coordinate system
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Code author: Miroslav Hallo
% Charles University in Prague, Faculty of Mathematics and Physics
% Web: http://geo.mff.cuni.cz/~hallo/
% E-mail: hallo@karel.troja.mff.cuni.cz
% Revision 6/2017: The first version of the function.
%
% Copyright (C) 2017  Miroslav Hallo
%
% This program is published under the GNU General Public License (GNU GPL).
%
% This program is free software: you can modify it and/or redistribute it
% or any derivative version under the terms of the GNU General Public
% License as published by the Free Software Foundation, either version 3
% of the License, or (at your option) any later version.
%
% This code is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY. We would like to kindly ask you to acknowledge the authors
% and don't remove their names from the code.
%
% You should have received copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
addpath([pwd,'/MATLAB']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUT

% location of stations in WGS84 Latitude Longitude (in degrees) Elevation (in kilometers)
LatLonEle = [
35.0094  135.7304  -0.674
34.7159  135.5199  -1.002
34.8983  135.7461  -0.791
34.6628  135.3896  -1.993
34.7630  135.7052  0.040
35.0469  135.4846  0.073
34.9286  135.4677  0.640
35.0713  135.5091  0.180
34.7785  135.6970  0.210
34.8621  135.5694  0.138
35.0625  135.7630  0.180
34.8521  135.9048  0.290
34.6568  135.6822  0.260
34.8635  135.5706  0.138
34.8209  135.3330  -0.211
34.6591  135.5118  -0.526
34.8900  135.3744  -0.736
];

% location of reference point in WGS84 Latitude Longitude (in degrees) Elevation (in kilometers)
Ref_LatLonEle = [34.844 135.622 0.140];

% END INPUT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Loc] = locX(LatLonEle(:,1),LatLonEle(:,2),Ref_LatLonEle(1),Ref_LatLonEle(2));

xyz = zeros(length(LatLonEle(:,1)),3);
xyz(:,1:2) = Loc;
xyz(:,3) = LatLonEle(:,3)-Ref_LatLonEle(3);

% station position in local coordinate system in kilometers
display('Position of stations in the local coordinate system in kilometers')
display('(Easting, Northing, Elevation)')

xyz





return