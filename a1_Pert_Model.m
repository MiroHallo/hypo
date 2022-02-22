%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Perturb velocity model to get dependence of sigma on distance for P and
% S wave arrival times.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Method by Hallo et al. (2019)
% Hallo,M., Oprsal,I., Asano,K., Gallovic,F. (2019): Seismotectonics of the 2018
%      Northern Osaka M6.1 earthquake and its aftershocks: joint
%      movements on strike-slip and reverse faults in inland Japan, Earth,
%      Planets and Space, 71:34.
%
% Code author: Miroslav Hallo
% Charles University in Prague, Faculty of Mathematics and Physics
% Web: http://geo.mff.cuni.cz/~hallo/
% E-mail: hallo@karel.troja.mff.cuni.cz
% Revision 7/2018: The first version of the function.
%
% Copyright (C) 2018  Miroslav Hallo
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

% Velocity model
crustName = 'crustal.dat';

Nmod = 1000; % Number of perturbed models
PrMod = 10; % P/S wave velocity perturbation in % (1 sigma)
PrDep = 10; % 1D layer depth perturbation in % (1 sigma)

Pmod = 300; % Number of models to plot
histXlim = [1 10]; % Limiths in seconds for hodochrones-histogram plot

DepthS = 10000; % Depth of source [m]
DepthR = 0; % Depth of receivers [m]
DistX = 0:2000:30000; % Source-receiver distances to compute rays
DX = 10; %Tolerance of ray for source-receiver distances [m]

% END INPUT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prepare perturbed velocity models

% Read crustal velocity model (top km format)
[veloc_mod0] = rCrustal(crustName);
dept = veloc_mod0(:,1)*1000;
alfa0 = veloc_mod0(:,2)*1000;
beta0 = veloc_mod0(:,3)*1000;
dens = veloc_mod0(:,4)*1000;

% Number of layers in the velocity model
nLayers = length(veloc_mod0(:,1));

% Number of distances
nDist = length(DistX);

% The first (correct) model
vm_changeP(1,:) = alfa0;
vm_changeS(1,:) = beta0;
vm_changeZ(1,:) = dept;

% Tarantola formula 6.34
G0 = (beta0.^2).*dens;
K0 = (alfa0.^2-(4/3)*beta0.^2).*dens;

% Create perturbed velocity models
for mod = 2:Nmod
    
    % Perturb bulk modulus
    G = G0.*exp(randn(nLayers,1)*(PrMod/100));
    K = K0.*exp(randn(nLayers,1)*(PrMod/100));
    
    % Tarantola formula 6.34
    alfa = sqrt((K+(4/3).*G)./dens);
    beta = sqrt(G./dens);
    
    vm_changeP(mod,1:nLayers) = alfa;
    vm_changeS(mod,1:nLayers) = beta;
    
    % change depths of layers
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

% P-wave models
figure('units','normalized','color',[1 1 1],'OuterPosition',[0.1,0.1,0.8,0.8]);
subplot(1,2,1)
title('Perturbed velocity models')
box on;
cc = colormap(['summer(',num2str(Pmod),')']);
botDmin = 9999999999;
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
            line([v/1000 v2/1000],[botD/1000 botD/1000],'color',cc(mod,:))
        end
        line([v/1000 v/1000],[topD/1000 botD/1000],'color',cc(mod,:))
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
        line([v/1000 v2/1000],[botD/1000 botD/1000],'color','k','LineWidth',2)
    end
    line([v/1000 v/1000],[topD/1000 botD/1000],'color','k','LineWidth',2)
end
set(gca,'YLim',[0 botDmin/1000]);
set(gca,'YDir','Reverse')
xlabel('P-wave velocity [km/s]')
ylabel('Depth [km]')

% S-wave models
subplot(1,2,2)
title('Perturbed velocity models')
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
        line([v/1000 v2/1000],[botD/1000 botD/1000],'color','k','LineWidth',2)
    end
    line([v/1000 v/1000],[topD/1000 botD/1000],'color','k','LineWidth',2)
end
set(gca,'YLim',[0 botDmin/1000]);
set(gca,'YDir','Reverse')
xlabel('S-wave velocity [km/s]')
ylabel('Depth [km]')

% Ray tracing
TTp = zeros(Nmod,nDist);
TTs = zeros(Nmod,nDist);
for mod = 1:Nmod
    DeptB = [vm_changeZ(mod,2:nLayers), vm_changeZ(mod,nLayers)+10000];
    [~,time] = ray1d(DeptB,vm_changeP(mod,1:end),DepthS,DepthR,DistX,DX);
    TTp(mod,1:nDist) = time';
    [~,time] = ray1d(DeptB,vm_changeS(mod,1:end),DepthS,DepthR,DistX,DX);
    TTs(mod,1:nDist) = time';
end

% Plot hodochrones-histograms
figure('units','normalized','color',[1 1 1],'OuterPosition',[0.1,0.1,0.8,0.8]);
for st = 1:nDist
    % P-waves
    subaxis(nDist,2,(st-1)*2+1,'Spacing',0.005,'Padding',0.0);
    hist(TTp(:,st))
    set(get(gca,'child'),'FaceColor',[0.5 0.5 0.5],'EdgeColor','k');
    if st<nDist
        set(gca,'xTickLabel',{})
    end
    set(gca,'yTickLabel',{})
    ylabel([num2str(DistX(st)/1000),' km'])
    set(get(gca,'ylabel'),'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right')
    if st==nDist
        xlabel('P-wave arrival [s]')
    elseif st==1
        title('Hodochrones by histograms')
    end
    xlim(histXlim)
    % S-waves
    subaxis(nDist,2,(st-1)*2+2,'Spacing',0.005,'Padding',0.0);
    hist(TTs(:,st))
    set(get(gca,'child'),'FaceColor',[0.5 0.5 0.5],'EdgeColor','k');
    if st<nDist
        set(gca,'xTickLabel',{})
    end
    set(gca,'yTickLabel',{})
    if st==nDist
        xlabel('S-wave arrival [s]')
    elseif st==1
        title('Hodochrones by histograms')
    end
    xlim(histXlim)
end

% Compute standart deviations sigma1 and plot
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
figure('units','normalized','color',[1 1 1],'OuterPosition',[0.1,0.1,0.8,0.8]);
subplot(1,2,1)
plot(DistX/1000,sigma1p,'r')
xlabel('Distance [km]')
ylabel('\sigma_P [s]')
title('Wave arrival time uncertainty')
subplot(1,2,2)
plot(DistX/1000,sigma1s,'g')
xlabel('Distance [km]')
ylabel('\sigma_S [s]')
title('Wave arrival time uncertainty')

% Compute roots and plot fits
csi = 1;
bp = polyfit(DistX(csi:end)/1000,sigma1p(csi:end), 3);
synt1p = polyval(bp, DistX/1000);
bs = polyfit(DistX(csi:end)/1000,sigma1s(csi:end), 3);
synt1s = polyval(bs, DistX/1000);
% Display results
disp(['sigma_P[s] = D[km]^3 *',num2str(bp(1)),' + D[km]^2 *',num2str(bp(2)),' + D[km] *',num2str(bp(3)),' + ',num2str(bp(4))])
disp(['sigma_S[s] = D[km]^3 *',num2str(bs(1)),' + D[km]^2 *',num2str(bs(2)),' + D[km] *',num2str(bs(3)),' + ',num2str(bs(4))])

figure('units','normalized','color',[1 1 1],'OuterPosition',[0.1,0.1,0.8,0.8]);
subplot(1,2,1)
plot(DistX/1000,sigma1p,'r',DistX/1000,synt1p,'k')
xlabel('Distance [km]')
ylabel('\sigma_P [s]')
legend('Simulated','Polyfit')
title('Fitted arrival time uncertainty')
subplot(1,2,2)
plot(DistX/1000,sigma1s,'g',DistX/1000,synt1s,'k')
xlabel('Distance [km]')
ylabel('\sigma_S [s]')
legend('Simulated','Polyfit')
title('Fitted arrival time uncertainty')





return