% RUN2_HYPO_LOC Search in 3D space for the earthquake hypocenter.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Performs an exhaustive search in 3D space for the earthquake hypocenter, 
% supporting both homogeneous and layered 1D velocity models. The location 
% uncertainty is evaluated within a rigorous Bayesian framework, utilizing 
% arrival time uncertainties previously computed by run1_pert_model script.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Miroslav HALLO
% Charles University in Prague, Faculty of Mathematics and Physics
% E-mail: hallo@karel.troja.mff.cuni.cz
% Revision 2017/06: First version
% Revision 2018/12: Enhanced version
% Revision 2019/03: Extended documentation
% Revision 2026/04: New Matlab version
% Tested in Matlab R2025b
% Method:
% Tarantola, A. (2005, Chapter 7.1): Inverse Problem Theory and Methods for
%     Model Parameter Estimation, Society for Industrial and Applied 
%     Mathematics, Philadelphia, USA.
% Hallo,M., Oprsal,I., Asano,K., Gallovic,F. (2019): Seismotectonics of the
%     2018 Northern Osaka M6.1 earthquake and its aftershocks: joint
%     movements on strike-slip and reverse faults in inland Japan, Earth,
%     Planets and Space, 71:34. https://doi.org/10.1186/s40623-019-1016-8
%
% Copyright (C) 2017-2019 Miroslav Hallo
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
% INIT:
close all;
clearvars;
projRoot = fileparts(which(mfilename));
addpath(fullfile(projRoot, 'lib'));
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:

% Reference point in Latitude, Longitude [deg], Elevation [km] (positive up)
Ref_LatLonEle = [34.844 135.622 0.140];

% Input text file with wave arrival time data
% [Latitude[deg], Longitude[deg], Elevation[km], P-wave[s], S-wave[s]]
% Note: If any arrival time is missing, add a negative number (e.g. -1234)
loc_file = 'example_loc.txt';

% Input text file with wave arrival uncertainty (polynomial of 3rd degree)
% Sigma[s] = A * Dist[km]^3 * + B * Dist[km]^2 + C* Dist[km] + D
% First line for the P-wave arrival, the second line for the S-wave arrival
% Note: This text file is an output from the run1_pert_model.m script
% Note: The parser ignore commented lines
unc_file = 'example_pert_model_uncertainty.txt';

% Set the grid search in kilometers (Z direction is depth = positive down)
gridX = -2.0 : 0.1 : 1.0;
gridY = -2.5 : 0.1 : 0.5;
gridZ = 8.5 : 0.1 : 11.5;

% Layered velocity model in the DWN format
% 3rd line: number of layers; from 6th line: data
% Depth of layer top[km]  Vp[km/s]  Vs[km/s]  Rho[g/cm^3]  Qp[-]  Qs[-]
crustName = 'example_crustal.dat';

% Compute in homogeneous model (layered=0) or 1D layered model (layered=1)
layered = 1;

% In the case of homogeneous model, set seismic wave velocities [km/s]
vp = 6.4; % P-wave velocity (homogeneous model)
vs = vp/1.73; % S-wave velocity (homogeneous model)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create directory for results
resDir = fullfile(projRoot, 'results');
if ~exist(resDir, 'dir')
    mkdir(resDir);
end

% Prepare timestamp
timestamp = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
outfile = [char(timestamp),'_hypo_loc'];

% -------------------------------------------------------------------
% Read input file with wave arrival time data
try
    fid = fopen(loc_file,'r');
    data = textscan(fid,'%f %f %f %f %f %[^\n]', 'CommentStyle', '#');
    fclose(fid);
catch
    data = [];
end
total = length(data{1,1});

% Get location of stations in kilometers (x = easting; y = northing; z = elevation)
[Loc] = locX(data{1,1}(:,1), data{1,2}(:,1), Ref_LatLonEle(1), Ref_LatLonEle(2));

% Get relative elevation [km] (positive up)
xyz(1:total,1:2) = Loc;
xyz(1:total,3) = data{1,3}(:,1) - Ref_LatLonEle(3);

% P-wave arrival times [s], if missing use -1234
tp = data{1,4}(:,1)';

% S-wave arrival times [s], if missing use -1234
ts = data{1,5}(:,1)';

% -------------------------------------------------------------------
% Read file with wave arrival uncertainty (polynomial of 3rd degree)
try
    fid = fopen(unc_file,'r');
    uncdata = textscan(fid,'%f %f %f %f %[^\n]', 'CommentStyle', '#');
    fclose(fid);
catch
    uncdata = [];
end
if length(uncdata{1,1}) < 2
    disp('ERROR: File with arrival uncertainty must have 2 data lines (P and S waves)')
    return
end
% Coefficients for a polynomial of 3rd degree for arrival times uncertainties
pCoeff = [uncdata{1,1}(1,1) uncdata{1,2}(1,1) uncdata{1,3}(1,1) uncdata{1,4}(1,1)];
sCoeff = [uncdata{1,1}(2,1) uncdata{1,2}(2,1) uncdata{1,3}(2,1) uncdata{1,4}(2,1)];

% -------------------------------------------------------------------
% Read velocity models (only for location in the 1D layered model)
if layered == 1
    [veloc_mod] = rCrustal(crustName);
    VelMod.nlay = length(veloc_mod(:,1));
    VelMod.depth(1:VelMod.nlay) = [veloc_mod(2:VelMod.nlay,1);veloc_mod(VelMod.nlay,1)+10];
    VelMod.Vp = veloc_mod(1:VelMod.nlay,2); % P-waves velocities [km/s]
    VelMod.Vs = veloc_mod(1:VelMod.nlay,3); % S-waves velocities [km/s]
end


%% -------------------------------------------------------------------
% Computing PDF in the 3D space

% Number of stations
N = length(xyz(:,1));

PDF = zeros(length(gridX), length(gridY), length(gridZ));
tp_synth = zeros(1,N);
ts_synth = zeros(1,N);

% Used arrivals (arrival time >= 0)
st_p_used = tp >= 0;
st_s_used = ts >= 0;

for x = 1 : length(gridX)
    disp(['Done: ',num2str(100* (x/length(gridX)),'%4.1f' ),'%'])
    for y = 1 : length(gridY)
        for z = 1 : length(gridZ)
            
            diff = zeros(1,N);
            
            % Compute synthetic times
            if layered == 1 % 1D layered model
                for st = 1 : N
                    DistX = sqrt( (gridX(x)-xyz(st,1))^2 + (gridY(y)-xyz(st,2))^2 );
                    [~,tp_synth(st)] = TT(VelMod.depth*1000,VelMod.Vp*1000,...
                                          gridZ(z)*1000,-xyz(st,3)*1000,DistX*1000,10);
                    [~,ts_synth(st)] = TT(VelMod.depth*1000,VelMod.Vs*1000,...
                                          gridZ(z)*1000,-xyz(st,3)*1000,DistX*1000,10);
                end
            else  % homogeneous model
                for st = 1 : N
                    tp_synth(st) = (1/vp)*sqrt( (gridX(x)-xyz(st,1))^2 + ...
                                   (gridY(y)-xyz(st,2))^2 + (gridZ(z)+xyz(st,3))^2 );
                    ts_synth(st) = (1/vs)*sqrt( (gridX(x)-xyz(st,1))^2 + ...
                                   (gridY(y)-xyz(st,2))^2 + (gridZ(z)+xyz(st,3))^2 );
                end
            end

            % Average origin time
            orig_t = mean([tp(st_p_used) - tp_synth(st_p_used), ...
                           ts(st_s_used) - ts_synth(st_s_used)]);
            
            sp = zeros(1,N);
            ss = zeros(1,N);
            for st = 1 : N
                % Find sigma
                DistX = sqrt( (gridX(x)-xyz(st,1))^2 + (gridY(y)-xyz(st,2))^2 );
                sp(st) = polyval(pCoeff(1,:), DistX);
                ss(st) = polyval(sCoeff(1,:), DistX);
                
                % P-waves
                if tp(st) >= 0
                    diff(st) = ((tp(st) - tp_synth(st) - orig_t)^2)/(sp(st)^2);
                end
                
                % S-waves
                if ts(st) >= 0
                    diff(st) = diff(st) + ((ts(st)-ts_synth(st)-orig_t)^2)/(ss(st)^2);
                end
            end
            
            PDF(x,y,z) = exp((-1/2) * sum(diff) );
            
        end
    end
end
% normalization to 1
PDF = PDF/sum(PDF(:));


%% -------------------------------------------------------------------
% Evaluate results
max_PDF = max(PDF(:));
for x = 1 : length(gridX)
    for y = 1 : length(gridY)
        for z = 1 : length(gridZ)
            if max_PDF == PDF(x,y,z)
                loc_index(1) = x;
                loc_index(2) = y;
                loc_index(3) = z;
            end
        end
    end
end

% The most likelyhood solution as max PDF
loc_res(1) = gridX(loc_index(1));
loc_res(2) = gridY(loc_index(2));
loc_res(3) = gridZ(loc_index(3));

% Convert solution back into Lat Lon
[Loc] = locWGS84(loc_res(1), loc_res(2), Ref_LatLonEle(1), Ref_LatLonEle(2));
loc_res_WGS(1:2) = Loc;
loc_res_WGS(3) = loc_res(3);

disp('----------------------------------')
disp('The maximum likelihood solution:')
disp('(Easting, Northing, Depth (km))')
disp(num2str(loc_res,'%5.1f'))
disp('(Lat, Lon, Depth (km))')
disp(num2str(loc_res_WGS,'%10.4f'))
disp('----------------------------------')


%% -------------------------------------------------------------------
% Evaluate misfit
st_dist = zeros(1,N);
tp_synth = zeros(1,N);
ts_synth = zeros(1,N);
tp_missfit = zeros(1,N);
ts_missfit = zeros(1,N);

% Compute synthetic times
if layered == 1 % 1D layered model
    for st = 1 : N
        DistX = sqrt( (loc_res(1)-xyz(st,1))^2 + (loc_res(2)-xyz(st,2))^2 );
        [~,tp_synth(st)] = TT(VelMod.depth*1000,VelMod.Vp*1000,...
                              loc_res(3)*1000,-xyz(st,3)*1000,DistX*1000,10);
        [~,ts_synth(st)] = TT(VelMod.depth*1000,VelMod.Vs*1000,...
                              loc_res(3)*1000,-xyz(st,3)*1000,DistX*1000,10);
    end
else % homogeneous model
    for st = 1 : N
        st_dist = sqrt( (loc_res(1)-xyz(st,1))^2 + ...
                        (loc_res(2)-xyz(st,2))^2 + (loc_res(3)+xyz(st,3))^2 );
        tp_synth(st) = (1/vp)*st_dist;
        ts_synth(st) = (1/vs)*st_dist;
    end
end

% Average origin time
orig_t = mean([tp(st_p_used) - tp_synth(st_p_used), ts(st_s_used) - ts_synth(st_s_used)]);

for st = 1 : N
    % P-waves
    if tp(st) >= 0
        tp_missfit(st) = (tp(st) - tp_synth(st) - orig_t);
    end
    
    % S-waves
    if ts(st) >= 0
        ts_missfit(st) = (ts(st) - ts_synth(st) - orig_t);
    end
end


%% -------------------------------------------------------------------
% Solution uncertainty

% Marginal PDFs
marginal_z_PDF = sum(PDF,3);
marginal_zy_PDF = sum(marginal_z_PDF,2); % marginal PDF for X direction
marginal_zx_PDF = sum(marginal_z_PDF,1); % marginal PDF for Y direction

marginal_y_PDF = squeeze(sum(PDF,2));
marginal_yx_PDF = sum(marginal_y_PDF,1); % marginal PDF for Z direction

% Fit the marginal by Gauss (by log, polyfit) and find sigma for X direction
lny = log(marginal_zy_PDF);
coeffs = polyfit(gridX, lny', 2);
sigma = sqrt(-1/coeffs(1));
loc_xyz_sigma(1) = sigma;

% Fit the marginal by Gauss (by log, polyfit) and find sigma for Y direction
lny = log(marginal_zx_PDF);
coeffs = polyfit(gridY, lny, 2);
sigma = sqrt(-1/coeffs(1));
loc_xyz_sigma(2) = sigma;

% Fit the marginal by Gauss (by log, polyfit) and find sigma for Z direction
lny = log(marginal_yx_PDF);
coeffs = polyfit(gridZ, lny, 2);
sigma = sqrt(-1/coeffs(1));
loc_xyz_sigma(3) = sigma;

disp('Loc XYZ 2sigma uncertainty (km)')
disp(num2str(loc_xyz_sigma*2,'%6.2f'))
disp('----------------------------------')

% -------------------------------------------------------------------
% Save results into text file
fid = fopen(fullfile(resDir, [outfile,'.txt']),'w');
fprintf(fid,'%s\r\n','# The maximum likelihood solution for the hypocentre location');
fprintf(fid,'%s\r\n','# Latitude, Longitude, Depth[km], E[km], N[km], E_sigma[km], N_sigma[km], Z_sigma[km], E_2sigma[km], N_2sigma[km], Z_2sigma[km]');
fprintf(fid,'%10.5f %10.5f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\r\n', loc_res_WGS, loc_res(1:2), loc_xyz_sigma, loc_xyz_sigma.*2);
fclose(fid);
disp(['Results successfully saved in: ', outfile,'.txt']);


%% -------------------------------------------------------------------
% plot cross-section
close all
% Plot map of stations
fih(1) = figure('color','w');
hold on
for i = 1 : N
    if (st_p_used(i)==0) && (st_s_used(i)==0)
        plot(xyz(i,1),xyz(i,2),'^','color',[0.8 0.8 0.8],'MarkerSize',8,'LineWidth',1);
    else
        hli(1) = plot(xyz(i,1),xyz(i,2),'^k','MarkerSize',8,'LineWidth',1);
    end
end
hli(2) = plot(loc_res(1),loc_res(2),'*','Color',[0.8 0.2 0.2],'MarkerSize',8,'LineWidth',1);
text(gridX(end),loc_res(2), {[' ', num2str(loc_res_WGS(1),'%10.4f')], ...
                             [' ', num2str(loc_res_WGS(2),'%10.4f')], ...
                             [' ', num2str(loc_res_WGS(3),'%6.1f'),' km']},'Color',[0.8 0.2 0.2])
hli(3) = plot([gridX(1), gridX(end), gridX(end), gridX(1), gridX(1)],...
    [gridY(1), gridY(1), gridY(end), gridY(end), gridY(1)],'Color',[0.2 0.8 0.2],'LineWidth',1);
hold off
title('Station and the maximum likelihood solution', 'FontWeight', 'normal')
legend(hli,'Stations','Epicentre','Search area','Location','northeast')
xlabel('Easting (km)')
ylabel('Northing (km)')
box on;
axis equal

% plot horizontal slice
fih(2) = figure('color','w');
imagesc(gridX,gridY,squeeze(PDF(:,:,loc_index(3)))');
set(gca,'YDir','normal')
axis image;
colormap(dusk)
colorbar
title(['PDF cross-section at depth ',num2str(gridZ(loc_index(3)),'%6.1f'),' km'],...
       'FontWeight', 'normal')
xlabel('Easting (km)')
ylabel('Northing (km)')
box on;

% plot vertical slice W-E
fih(3) = figure('color','w');
imagesc(gridX,gridZ,squeeze(PDF(:,loc_index(2),:))');
set(gca,'YDir','normal')
axis image;
colormap(dusk)
colorbar
set(gca,'YDir','reverse')
title('PDF cross-section (W-E)', 'FontWeight', 'normal')
xlabel('Easting (km)')
ylabel('Depth (km)')
box on;

% plot vertical slice S-N
fih(4) = figure('color','w');
imagesc(gridY,gridZ,squeeze(PDF(loc_index(1),:,:))');
set(gca,'YDir','normal')
axis image;
colormap(dusk)
colorbar
set(gca,'YDir','reverse')
title('PDF cross-section (S-N)', 'FontWeight', 'normal')
xlabel('Northing (km)')
ylabel('Depth (km)')
box on;

% Plot misfits
fih(5) = figure('color','w');
hold on
for st = 1 : N
    % P-waves
    if tp(st) >= 0
        plot([tp(st) tp(st)],[st-1 st],'Color',[0.8 0.2 0.2])
        text(tp(st)+0.05,st-0.5,num2str(tp_missfit(st),'%5.2f'),'Color',[0.8 0.2 0.2])
    end
    plot([tp_synth(st)+orig_t tp_synth(st)+orig_t],[st-1 st],':','Color',[0.8 0.2 0.2])
    % S-waves
    if ts(st) >= 0
        plot([ts(st) ts(st)],[st-1 st],'Color',[0.2 0.8 0.2])
        text(ts(st)+0.05,st-0.5,num2str(ts_missfit(st),'%5.2f'),'Color',[0.2 0.8 0.2])
    end
    plot([ts_synth(st)+orig_t ts_synth(st)+orig_t],[st-1 st],':','Color',[0.2 0.8 0.2])
    % Station common
    plot([0 10],[st st],'color',[0.8 0.8 0.8])
    text(0,st-0.5,[' ST-',num2str(st)])
end
hold off
set(gca,'YDir','reverse')
set(gca, 'Layer', 'top')
set(gca,'Ylim',[0 N])
set(gca,'YTickLabel',{})
title('Data misfits (observed - synthetic)', 'FontWeight', 'normal')
xlabel('Time (s)')
box on;


%% -------------------------------------------------------------------
% Save figures
outfile_tmp = [outfile,'_map'];
exportgraphics(fih(1), fullfile(resDir, [outfile_tmp,'.png']), 'Resolution', 300);
fprintf('Figure successfully saved as: %s.png\n', outfile_tmp);

outfile_tmp = [outfile,'_pdf_hor'];
exportgraphics(fih(2), fullfile(resDir, [outfile_tmp,'.png']), 'Resolution', 300);
fprintf('Figure successfully saved as: %s.png\n', outfile_tmp);

outfile_tmp = [outfile,'_pdf_we'];
exportgraphics(fih(3), fullfile(resDir, [outfile_tmp,'.png']), 'Resolution', 300);
fprintf('Figure successfully saved as: %s.png\n', outfile_tmp);

outfile_tmp = [outfile,'_pdf_sn'];
exportgraphics(fih(4), fullfile(resDir, [outfile_tmp,'.png']), 'Resolution', 300);
fprintf('Figure successfully saved as: %s.png\n', outfile_tmp);

outfile_tmp = [outfile,'_misfit'];
exportgraphics(fih(5), fullfile(resDir, [outfile_tmp,'.png']), 'Resolution', 300);
fprintf('Figure successfully saved as: %s.png\n', outfile_tmp);

