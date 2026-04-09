% RUN1_PERT_MODEL Perturb velocity model by Monte Carlo simulations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Perturb velocity model by Monte Carlo simulations to get dependence 
% of uncertainty (1 sigma) on distance for P and S wave arrival times.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Miroslav HALLO
% Charles University in Prague, Faculty of Mathematics and Physics
% E-mail: hallo@karel.troja.mff.cuni.cz
% Revision 2018/07: First version
% Revision 2019/03: Enhanced version
% Revision 2026/04: New Matlab version
% Tested in Matlab R2025b
% Method:
% Hallo,M., Oprsal,I., Asano,K., Gallovic,F. (2019): Seismotectonics of the 
%      2018 Northern Osaka M6.1 earthquake and its aftershocks: joint
%      movements on strike-slip and reverse faults in inland Japan, Earth,
%      Planets and Space, 71:34. https://doi.org/10.1186/s40623-019-1016-8
%
% Copyright (C) 2018,2019 Miroslav Hallo
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

% Velocity model
crustName = 'example_crustal.dat';

Nmod = 1000; % Number of perturbed models
PrMod = 10; % P/S wave velocity perturbation in [%] (1 sigma)
PrDep = 10; % 1D layer depth perturbation in [%] (1 sigma)

Pmod = 300; % Number of models to plot
histXlim = [1 20]; % Limiths in seconds for hodochrones-histogram plot

DepthS = 10000; % Depth of source [m]
DepthR = 0; % Depth of receivers [m]
DistX = 0:2000:30000; % Source-receiver distances to compute rays
DX = 10; % Tolerance of ray for source-receiver distances [m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create directory for results
resDir = fullfile(projRoot, 'results');
if ~exist(resDir, 'dir')
    mkdir(resDir);
end

% Prepare timestamp
timestamp = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
outfile = [char(timestamp),'_pert_model'];

% -------------------------------------------------------------------
% Prepare perturbed velocity models

% Read crustal velocity model (top km format)
[veloc_mod0] = rCrustal(crustName);
dept = veloc_mod0(:,1) * 1000;
alfa0 = veloc_mod0(:,2) * 1000;
beta0 = veloc_mod0(:,3) * 1000;
dens = veloc_mod0(:,4) * 1000;

% Number of layers in the velocity model
nLayers = length(veloc_mod0(:,1));

% Number of distances
nDist = length(DistX);

% The first (correct) model
vm_changeP(1,:) = alfa0;
vm_changeS(1,:) = beta0;
vm_changeZ(1,:) = dept;

% Tarantola (2005) formula 6.34
G0 = (beta0.^2).*dens;
K0 = (alfa0.^2-(4/3)*beta0.^2).*dens;

% -------------------------------------------------------------------
% Create perturbed velocity models
for mod = 2:Nmod
    
    % Perturb bulk modulus
    G = G0 .* exp(randn(nLayers,1)*(PrMod/100));
    K = K0 .* exp(randn(nLayers,1)*(PrMod/100));
    
    % Tarantola (2005) formula 6.34
    alfa = sqrt((K+(4/3).*G)./dens);
    beta = sqrt(G./dens);
    
    vm_changeP(mod,1:nLayers) = alfa;
    vm_changeS(mod,1:nLayers) = beta;
    
    % Change depths of layers
    vm_changeZ(mod,:) = dept;
    for l = 2:nLayers
        SQ0 = vm_changeZ(mod,l)^2;
        SQ = SQ0*exp(randn(1,1)*(PrDep/100));
        Z = sqrt(SQ);
        if Z <= vm_changeZ(mod,l-1)
            Z = vm_changeZ(mod,l);
        elseif (l<nLayers) && (Z>=vm_changeZ(mod,l+1))
            Z = vm_changeZ(mod,l);
        end
        vm_changeZ(mod,l) = Z;
    end
end


%% -------------------------------------------------------------------
% Plot ensemble of perturbed models
fig = figure('Color','w');
cc = colormap(['copper(',num2str(Pmod),')']);
randcolol = max(1,(1:Pmod)' - randi(floor(Pmod/3),Pmod,1));
cc = cc(randcolol,:);

% P-wave models
subplot(1,2,1)
box on;
botDmin = 1e9;
for mod = 2:Pmod
    for l = 1:nLayers
        v = vm_changeP(mod,l);
        topD = vm_changeZ(mod,l);
        if l==nLayers
            botD = vm_changeZ(mod,l) + (vm_changeZ(mod,l)-vm_changeZ(mod,l-1));
            if botD<botDmin
                botDmin = botD;
            end
        else
            botD = vm_changeZ(mod,l+1);
            v2 = vm_changeP(mod,l+1);
            line([v v2]/1000,[botD botD]/1000,'color',cc(mod,:))
        end
        hln(2) = line([v v]/1000,[topD botD]/1000,'color',cc(mod,:));
    end
end
mod = 1;
for l = 1:nLayers
    v = vm_changeP(mod,l);
    topD = vm_changeZ(mod,l);
    if l==nLayers
        botD = vm_changeZ(mod,l) + (vm_changeZ(mod,l)-vm_changeZ(mod,l-1));
        if botD<botDmin
            botDmin = botD;
        end
    else
        botD = vm_changeZ(mod,l+1);
        v2 = vm_changeP(mod,l+1);
        line([v v2]/1000,[botD botD]/1000,'color','k','LineStyle','--','LineWidth',1.5)
    end
    hln(1) = line([v v]/1000,[topD botD]/1000,'color','k','LineStyle','--','LineWidth',1.5);
end
set(gca,'YLim',[0 botDmin/1000]);
set(gca,'YDir','Reverse')
legend(hln,'Input','Ensemble','Location','northeast');
xlabel('P-wave velocity (km/s)')
ylabel('Depth (km)')

% S-wave models
subplot(1,2,2)
box on;
for mod = 2:Pmod
    for l = 1:nLayers
        v = vm_changeS(mod,l);
        topD = vm_changeZ(mod,l);
        if l==nLayers
            botD = vm_changeZ(mod,l) + (vm_changeZ(mod,l)-vm_changeZ(mod,l-1));
            if botD<botDmin
                botDmin = botD;
            end
        else
            botD = vm_changeZ(mod,l+1);
            v2 = vm_changeS(mod,l+1);
            line([v/1000 v2/1000],[botD/1000 botD/1000],'color',cc(mod,:))
        end
        line([v/1000 v/1000],[topD/1000 botD/1000],'color',cc(mod,:))
    end
end
mod = 1; % average model
for l = 1:nLayers
    v = vm_changeS(mod,l);
    topD = vm_changeZ(mod,l);
    if l==nLayers
        botD = vm_changeZ(mod,l) + (vm_changeZ(mod,l)-vm_changeZ(mod,l-1));
        if botD<botDmin
            botDmin = botD;
        end
    else
        botD = vm_changeZ(mod,l+1);
        v2 = vm_changeS(mod,l+1);
        line([v/1000 v2/1000],[botD/1000 botD/1000],'color','k','LineStyle','--','LineWidth',1.5)
    end
    line([v/1000 v/1000],[topD/1000 botD/1000],'color','k','LineStyle','--','LineWidth',1.5)
end
set(gca,'YLim',[0 botDmin/1000]);
set(gca,'YDir','Reverse')
xlabel('S-wave velocity (km/s)')
ylabel('Depth (km)')

% Save figure
outfile_tmp = [outfile,'_ensemble'];
exportgraphics(fig, fullfile(resDir, [outfile_tmp,'.png']), 'Resolution', 300);
fprintf('Figure successfully saved as: %s.png\n', outfile_tmp);


%% -------------------------------------------------------------------
% Ray tracing
TTp = zeros(Nmod,nDist);
TTs = zeros(Nmod,nDist);
for mod = 1:Nmod
    DeptB = [vm_changeZ(mod,2:nLayers), vm_changeZ(mod,nLayers)+10000];
    [~,time] = TT(DeptB,vm_changeP(mod,1:end),DepthS,DepthR,DistX,DX);
    TTp(mod,1:nDist) = time';
    [~,time] = TT(DeptB,vm_changeS(mod,1:end),DepthS,DepthR,DistX,DX);
    TTs(mod,1:nDist) = time';
end


%% -------------------------------------------------------------------
% Plot hodochrones-histograms
fig = figure('Color','w');
for st = 1:nDist
    % P-waves
    subaxis(nDist,2,(st-1)*2+1,'Spacing',0.005,'Padding',0.0);
    histogram(TTp(:,st),'FaceColor',[0.2 0.2 0.2],'EdgeColor',[0.2 0.2 0.2],'FaceAlpha',1)
    if st<nDist
        set(gca,'xTickLabel',{})
    end
    set(gca,'yTickLabel',{})
    ylabel([num2str(DistX(st)/1000),' km'])
    set(get(gca,'ylabel'),'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right')
    if st==nDist
        xlabel('P-wave arrival (s)')
    end
    xlim(histXlim)
    % S-waves
    subaxis(nDist,2,(st-1)*2+2,'Spacing',0.005,'Padding',0.0);
    histogram(TTs(:,st),'FaceColor',[0.2 0.2 0.2],'EdgeColor',[0.2 0.2 0.2],'FaceAlpha',1)
    if st<nDist
        set(gca,'xTickLabel',{})
    end
    set(gca,'yTickLabel',{})
    if st==nDist
        xlabel('S-wave arrival (s)')
    end
    xlim(histXlim)
end

% Save figure
outfile_tmp = [outfile,'_hodochrones'];
exportgraphics(fig, fullfile(resDir, [outfile_tmp,'.png']), 'Resolution', 300);
fprintf('Figure successfully saved as: %s.png\n', outfile_tmp);


%% -------------------------------------------------------------------
% Compute standart deviations sigma 1 of arrival times
sigma1p = zeros(1,nDist);
sigma1s = zeros(1,nDist);
for mod = 1:PrMod
    for st = 1:nDist
        % P-waves
        tmp = TTp(:,st);
        smodCH = std(tmp);
        sigma1p(st) = smodCH;
        
        % S-waves
        tmp = TTs(:,st);
        smodCH = std(tmp);
        sigma1s(st) = smodCH;
    end
end

% Compute roots of uncertainty curves
csi = 1;
bp = polyfit(DistX(csi:end)/1000, sigma1p(csi:end), 3);
synt1p = polyval(bp, DistX/1000);
bs = polyfit(DistX(csi:end)/1000, sigma1s(csi:end), 3);
synt1s = polyval(bs, DistX/1000);

% Plot figure
fig = figure('Color','w');
subplot(1,2,1)
plot(DistX/1000,sigma1p,'Color',[0.8 0.2 0.2]); hold on
plot(DistX/1000,synt1p,'Color','k','LineStyle',':')
legend('Simulated','Polynomial')
xlabel('Distance (km)')
ylabel('P-wave arrival time uncertainty \sigma_P (s)')
subplot(1,2,2)
plot(DistX/1000,sigma1s,'Color',[0.2 0.8 0.2]); hold on
plot(DistX/1000,synt1s,'Color','k','LineStyle',':');
legend('Simulated','Polynomial')
xlabel('Distance (km)')
ylabel('S-wave arrival time uncertainty \sigma_S (s)')

% Save figure
outfile_tmp = [outfile,'_uncertainty'];
exportgraphics(fig, fullfile(resDir, [outfile_tmp,'.png']), 'Resolution', 300);
fprintf('Figure successfully saved as: %s.png\n', outfile_tmp);

% Display results
disp('--------------------------------')
disp('Resultant uncertainty dependence (polynomial of 3rd degree):')
disp(['sigma_P[s] = D[km]^3 *',num2str(bp(1)),' + D[km]^2 *',num2str(bp(2)),' + D[km] *',num2str(bp(3)),' + ',num2str(bp(4))])
disp(['sigma_S[s] = D[km]^3 *',num2str(bs(1)),' + D[km]^2 *',num2str(bs(2)),' + D[km] *',num2str(bs(3)),' + ',num2str(bs(4))])
disp('--------------------------------')

% Open text file for results
fid = fopen(fullfile(resDir, [outfile_tmp,'.txt']),'w');
fprintf(fid,'%s\r\n','# P-wave uncertainty (polynomial of 3rd degree):');
fprintf(fid,'%9.6e  %9.6e  %9.6e  %9.6e\r\n',bp);
fprintf(fid,'%s\r\n','# S-wave uncertainty (polynomial of 3rd degree):');
fprintf(fid,'%9.6e  %9.6e  %9.6e  %9.6e\r\n',bs);
fclose(fid);
disp(['Results successfully saved in: ', outfile_tmp,'.txt']);

