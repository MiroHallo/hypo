%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Grid search of the earthquake hypocenter in homogenous medium
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Method by Hallo et al. (2019)
% Hallo,M., Oprsal,I., Asano,K., Gallovic,F. (2019): Seismotectonics of the 2018
%      Northern Osaka M6.1 earthquake and its aftershocks: joint
%      movements on strike-slip and reverse faults in inland Japan, Earth,
%      Planets and Space, 71:34.
%
% Probability density function (PDF) following theory by Tarantola (2005, Chapter 7.1)
% Tarantola, A., 2005. Inverse Problem Theory and Methods for Model Parameter Estimation,
%      Society for Industrial and Applied Mathematics, Philadelphia.
%
% Code author: Miroslav Hallo
% Charles University in Prague, Faculty of Mathematics and Physics
% Web: http://geo.mff.cuni.cz/~hallo/
% E-mail: hallo@karel.troja.mff.cuni.cz
% Revision 6/2017: The first version of the function.
% Revision 12/2018: Enhanced version.
%
% Copyright (C) 2017,2018  Miroslav Hallo
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

% location of stations in kilometers (x = easting; y = northing; z = elevation)
xyz = [
    9.9036   18.4123   -0.8140
   -9.3280  -14.2601   -1.1420
   11.3380    6.0447   -0.9310
  -21.2324  -20.1712   -2.1330
    7.6013   -9.0169   -0.1000
  -12.5531   22.5868   -0.0670
  -14.0971    9.4177    0.5000
  -10.3147   25.3030    0.0400
    6.8521   -7.2915    0.0700
   -4.8056    2.0149   -0.0020
   12.8820   24.3234    0.0400
   25.8370    0.9017    0.1500
    5.5000  -20.8391    0.1200
   -4.6960    2.1707   -0.0020
  -26.4035   -2.5715   -0.3510
  -10.0680  -20.5830   -0.6660
  -22.6211    5.1207   -0.8760
];

% location of reference point in WGS84 Latitude Longitude (in degrees) Elevation (in kilometers)
Ref_LatLonEle = [34.844 135.622 0.140];

% P-wave arrival times in seconds, if missing use -1234
tp = [3.9983 2.9318 2.7002 4.6666 2.5732 4.6214 3.2787 4.8322 -1234 -1234 5.1884 4.7214 -1234 1.9000 4.5943 3.9344 3.9426];

% S-wave arrival times in seconds, if missing use -1234
ts = [6.6714 5.1230 4.8450 8.1520 -1234 8.0706 -1234 8.3990 -1234 -1234 -1234 8.2709 -1234 3.4341 8.0368 -1234 6.8460];

% Coefficients for a polynomial of 3rd degree for arrival times uncertainties
pCoeff = [-1.63394549384411e-06 0.000124393138622651 0.00033456939971057 0.0318835318434683]; % P - 10km depth
sCoeff = [-3.83041158964040e-06 0.000293120791393656 0.00077613429952614 0.0744937975783552]; % S - 10km depth

% Set grid for the grid search in kilometers (Z direction is depth here = positive sign down)
gridX = -5:0.5:5;
gridY = -5:0.5:5;
gridZ = 5:2:25;

% Seismic waves velocities in km/s (the ray tracing is in homogenous medium)
vp = 6.4;
vs = vp/1.73;

% END INPUT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%% Computing

% number of stations
N = length(xyz(:,1));

PDF = zeros(length(gridX),length(gridY),length(gridZ));
tp_synth = zeros(1,N);
ts_synth = zeros(1,N);

st_p_used = tp ~= -1234;
st_s_used = ts ~= -1234;

for x = 1 : length(gridX)
    for y = 1 : length(gridY)
        for z = 1 : length(gridZ)
            
            diff = zeros(1,N);
            
            % Compute synthetic times
            for st = 1 : N
                tp_synth(st) = (1/vp)*sqrt( (gridX(x)-xyz(st,1))^2 + (gridY(y)-xyz(st,2))^2 + (gridZ(z)+xyz(st,3))^2 );
                ts_synth(st) = (1/vs)*sqrt( (gridX(x)-xyz(st,1))^2 + (gridY(y)-xyz(st,2))^2 + (gridZ(z)+xyz(st,3))^2 );
            end
            
            % Average origin time
            orig_t = mean([tp(st_p_used) - tp_synth(st_p_used), ts(st_s_used) - ts_synth(st_s_used)]);
            
            for st = 1 : N
                % Find sigmas
                DistX = sqrt( (gridX(x)-xyz(st,1))^2 + (gridY(y)-xyz(st,2))^2 );
                sp(st) = polyval(pCoeff(1,:), DistX);
                ss(st) = polyval(sCoeff(1,:), DistX);
                
                % P-waves
                if tp(st) ~= -1234
                    diff(st) = ((tp(st) - tp_synth(st) - orig_t)^2)/(sp(st)^2);
                end
                
                % S-waves
                if ts(st) ~= -1234
                    diff(st) = diff(st) + ((ts(st) - ts_synth(st) - orig_t)^2)/(ss(st)^2);
                end
            end
            
            PDF(x,y,z) = exp((-1/2) * sum(diff) );
            
        end
    end
end
% normalization to 1
PDF = PDF/sum(PDF(:));
    
%% Results

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

% The most likelihood solution as max PDF
loc_res(1) = gridX(loc_index(1));
loc_res(2) = gridY(loc_index(2));
loc_res(3) = gridZ(loc_index(3));

% Convert solution back into Lat Lon
[Loc] = locWGS84(loc_res(1),loc_res(2),Ref_LatLonEle(1),Ref_LatLonEle(2));
loc_res_WGS(1:2) = Loc;
loc_res_WGS(3) = loc_res(3);

disp('----------------------------------')
disp('The maximum likelihood solution:')
disp('(Easting, Northing, Depth (km))')
disp(num2str(loc_res,'%5.1f'))
disp('(Lat, Lon, Depth (km))')
disp(num2str(loc_res_WGS,'%10.4f'))
disp('----------------------------------')

%% Misfit

st_dist = zeros(1,N);
tp_synth = zeros(1,N);
ts_synth = zeros(1,N);
tp_missfit = zeros(1,N);
ts_missfit = zeros(1,N);

% Compute synthetic times
for st = 1 : N
    st_dist = sqrt( (loc_res(1)-xyz(st,1))^2 + (loc_res(2)-xyz(st,2))^2 + (loc_res(3)+xyz(st,3))^2 );
    tp_synth(st) = (1/vp)*st_dist;
    ts_synth(st) = (1/vs)*st_dist;
end

% Average origin time
orig_t = mean([tp(st_p_used) - tp_synth(st_p_used), ts(st_s_used) - ts_synth(st_s_used)]);

for st = 1 : N
    % P-waves
    if tp(st) ~= -1234
        tp_missfit(st) = (tp(st) - tp_synth(st) - orig_t);
    end
    
    % S-waves
    if ts(st) ~= -1234
        ts_missfit(st) = (ts(st) - ts_synth(st) - orig_t);
    end
end

%% Solution uncertainty

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


%% Cross-section plots

% Plot map of stations
fih(1) = figure('color','w');
hold on
for i = 1 : N
    if (st_p_used(i)==0) && (st_s_used(i)==0) % station not used
        plot(xyz(i,1),xyz(i,2),'^','color',[0.8 0.8 0.8]);
    else
        plot(xyz(i,1),xyz(i,2),'^k');
    end
end
plot(loc_res(1),loc_res(2),'*r');
text(gridX(end),loc_res(2), {num2str(loc_res_WGS(1),'%10.4f'),num2str(loc_res_WGS(2),'%10.4f'),[num2str(loc_res_WGS(3),'%6.2f'),' km']},'color','r')
plot([gridX(1),gridX(end),gridX(end),gridX(1),gridX(1)],...
    [gridY(1),gridY(1),gridY(end),gridY(end),gridY(1)],'g');
hold off
title('Stations network with the maximum likelihood solution')
xlabel('Easting (km)')
ylabel('Northing (km)')
box on;
axis equal

% plot horizontal slice
fih(2) = figure('color','w');
h = pcolor(gridX,gridY,squeeze(PDF(:,:,loc_index(3)))');
set(h, 'EdgeColor', 'none');
colormap(flipud(hot))
colorbar
set(gca,'Color',[0.8 0.8 0.8])
title(['PDF: cross-section at depth ',num2str(gridZ(loc_index(3)),'%6.2f'),' km'])
xlabel('Easting (km)')
ylabel('Northing (km)')
box on;
axis equal

% plot vertical slice W-E
fih(3) = figure('color','w');
h = pcolor(gridX,gridZ,squeeze(PDF(:,loc_index(2),:))');
set(h, 'EdgeColor', 'none');
colormap(flipud(hot))
colorbar
set(gca,'Color',[0.8 0.8 0.8])
set(gca,'YDir','reverse')
title('PDF: W-E vertical cross-section')
xlabel('Easting (km)')
ylabel('Depth (km)')
box on;
axis equal

% plot vertical slice S-N
fih(4) = figure('color','w');
h = pcolor(gridY,gridZ,squeeze(PDF(loc_index(1),:,:))');
set(h, 'EdgeColor', 'none');
colormap(flipud(hot))
colorbar
set(gca,'Color',[0.8 0.8 0.8])
set(gca,'YDir','reverse')
title('PDF: S-N vertical cross-section')
xlabel('Northing (km)')
ylabel('Depth (km)')
box on;
axis equal

% Plot misfits
fih(5) = figure('color','w');
hold on
for st = 1 : N
    % P-waves
    if tp(st) ~= -1234
        plot([tp(st) tp(st)],[st-1 st],'r')
        text(tp(st)+0.05,st-0.5,num2str(tp_missfit(st),'%5.2f'),'color','r')
    end
    plot([tp_synth(st)+orig_t tp_synth(st)+orig_t],[st-1 st],':r')
    % S-waves
    if ts(st) ~= -1234
        plot([ts(st) ts(st)],[st-1 st],'g')
        text(ts(st)+0.05,st-0.5,num2str(ts_missfit(st),'%5.2f'),'color','g')
    end
    plot([ts_synth(st)+orig_t ts_synth(st)+orig_t],[st-1 st],':g')
    % Station common
    plot([0 10],[st st],'color',[0.8 0.8 0.8])
    text(0,st-0.5,['ST',num2str(st)])
end
hold off
set(gca,'YDir','reverse')
set(gca,'Ylim',[0 N])
set(gca,'YTickLabel',{})
title('Picked and synthetic arrivals misfits')
xlabel('Time (s)')
box on;





return