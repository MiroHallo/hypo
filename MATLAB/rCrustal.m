% Read ASCII files with velocity model in DWN form
% Author: Miroslav Hallo (hallo.miroslav@gmail.com)

function [veloc_mod] = rCrustal(path)

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




