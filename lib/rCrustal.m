function [veloc_mod] = rCrustal(path)
% RCRUSTAL Read ASCII files with velocity model in DWN format.
% 
% Author: Miroslav HALLO
% Charles University in Prague, 2018/12
% 
% Copyright (C) 2018 Miroslav Hallo
% This program is published under the GNU General Public License (GNU GPL).
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
% path - filaname of the velocity model in DWN format
% 
% OUTPUT:
% veloc_mod(:,1) - Depth of layer top [km]
% veloc_mod(:,2) - Vp [km/s]
% veloc_mod(:,3) - Vs [km/s]
% veloc_mod(:,4) - Rho(g/c^3)
% veloc_mod(:,5) - Qp [-]
% veloc_mod(:,6) - Qs [-]
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileID = fopen(path,'r');

tline = fgets(fileID);
tline = fgets(fileID);

tline = fgets(fileID);
nLayers = cell2mat(textscan(tline,'%d'));

tline = fgets(fileID);
tline = fgets(fileID);

veloc_mod = zeros(nLayers,6);
for l =1:nLayers
    tline = fgets(fileID);
    veloc_mod(l,1:6) = cell2mat(textscan(tline,'%f %f %f %f %f %f'));
end

fclose(fileID);

end
