function [AI, TT] = TT(Depth, Vel, DepthS, DepthR, DistX, DX)
% TT Fastest travel time in the gradient layered model.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Fastest travel time in the gradient layered model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Author: Miroslav HALLO
% Revision 2017/01: First version
% 
% Copyright (C) 2017 Miroslav Hallo
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
% 
% INPUT:
%        Depth = depths of layers (bottom) [meters]
%        Vel = wave velocities [m/s]
%        DepthS = depth of the source [meters] (depth positive down)
%        DepthR = depth of the receiver [meters] (depth positive down)
%        DistX = vector of horizontal source-receiver distances [meters]
%        DX = Tolerance of final ray for source-receiver distances [meters]
%
% OUTPUT:
%        AI = Source take-off angles (0-180 [deg], where 0=UP, 180=DOWN)
%        TT = Travel time for all distances [seconds]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants

frac_ang = 4; % Lowering factor of the angle step search (4 is the best)
MaxLoop = int16(100); % Maximal number of loop cycles for searching of the ray

% number of layers in model
nLayers = int16(length(Depth));

% number of distances
nDist = int16(length(DistX));

% Optimalization for speed
if nDist <5
    irt_N = 4; % Number of Initial ray tracing angles (irt_step*(irt_N-1) must give 90)
    irt_step = 30; % Step of Initial ray tracing [deg]
elseif nDist <50
    irt_N = 19; % Number of Initial ray tracing angles (irt_step*(irt_N-1) must give 90)
    irt_step = 5; % Step of Initial ray tracing [deg]
elseif  nDist <500
    irt_N = 37; % Number of Initial ray tracing angles (irt_step*(irt_N-1) must give 90)
    irt_step = 2.5; % Step of Initial ray tracing [deg]
else
    irt_N = 91; % Number of Initial ray tracing angles (irt_step*(irt_N-1) must give 90)
    irt_step = 1; % Step of Initial ray tracing [deg]
end


%% -------------------------------------------------------------------
% Prepare

% Allocation of res
AI = zeros(nDist,1);
TT = zeros(nDist,1) + Inf;

% IS index of source Layer
IS = int16(find(Depth >= DepthS,1,'first'));
if isempty(IS)
    IS = nLayers;
end

% IR index of receiver Layer
IR = int16(find(Depth >= DepthR,1,'first'));
if isempty(IR)
    IR = nLayers;
end

% Layer thicknesses
Thick = zeros(nLayers,1);
Thick(1,1) = Inf;
Thick(2:end,1) = Depth(2:end) - Depth(1:end-1);

% find vertical distance to a boundary
if IS == 1
    ISresDepU = Inf; % vertical dist. to boundary at source layer
else
    ISresDepU =  DepthS - Depth(IS-1); % vertical dist. to boundary at source layer
end

if IR == 1
    IRresDepU = Inf; % vertical dist. to boundary at receiver layer
else
    IRresDepU = DepthR - Depth(IR-1); % vertical dist. to boundary at receiver layer
end

ISresDepD = Depth(IS) - DepthS; % vertical dist. to boundary at source layer
IRresDepD = Depth(IR) - DepthR; % vertical dist. to boundary at receiver layer

DepDel =  DepthS - DepthR; % delDepths (positive for deaper source)


%% -------------------------------------------------------------------
% Horizontal ray

if DepDel == 0
    TT = DistX./Vel(IS);
    AI = AI + 90;
    return
end


%% -------------------------------------------------------------------
% Initial ray tracing (blind shotting)

if DepDel > 0 % upgoing
    irtA = (0:irt_N-1)*irt_step;
else % downgoing
    irtA = 90 + (0:irt_N-1)*irt_step;
end

irtX = zeros(1,irt_N);
for i = 1 : irt_N
    [~,irtX(i),~] = TX(Vel,Thick,irtA(i),DepthS,DepthR,IS,IR,ISresDepD,ISresDepU,IRresDepD,IRresDepU);
end

%% -------------------------------------------------------------------
% Finishing ray tracing for all distances

for d = int16(1) : nDist
    
    % current distance
    di = DistX(d);
    
    % a-priory angle ---------------------------------
    if abs(irtX(1) - di) < DX % direct hit by first
        AI(d) = irtA(1);
        [TT(d),~,~] = TX(Vel,Thick,AI(d),DepthS,DepthR,IS,IR,ISresDepD,ISresDepU,IRresDepD,IRresDepU);
        continue
    elseif irtX(1) < di % typical for UP-going
        grad = 1; % normal state (higher distance = higher angle)
        chd = int16(1); % change search direction (1 = too small distance)
        AngI = find(irtX > di,1,'first');
        A = irtA(AngI-1);
    else % typical for DOWN-going
        grad = -1; % reverse state (higher distance = smaller angle)
        chd = int16(1); % change search direction (1 = too small distance)
        AngI = find(irtX <= di,1,'first');
        A = irtA(AngI);
    end
    
    AstepA = irt_step/frac_ang; % initial step for precisse ray search
    count = int16(0); % loop number
    

    while count <= MaxLoop
        count = count + int16(1);
        
        [Tn,Xn,~] = TX(Vel,Thick,A,DepthS,DepthR,IS,IR,ISresDepD,ISresDepU,IRresDepD,IRresDepU);
        
        % compare horizontal distance
        if abs(Xn - di) < DX % ray founded
            break
            
        else % set new Ang for next iteraction
            if chd == int16(1) % too small distance
                if Xn < di % still too small distance
                    A = A + grad*AstepA;
                else % too big distance
                    AstepA = AstepA/frac_ang;
                    chd = int16(-1);
                    A = A - grad*AstepA;
                end
                
            else % too big distance
                if Xn > di % still too big distance
                    A = A - grad*AstepA;
                else % too big distance
                    AstepA = AstepA/frac_ang;
                    chd = int16(1);
                    A = A + grad*AstepA;
                end
            end
        end
    end
    
    % Results
    AI(d) = A;
    TT(d) = Tn;
    
end

end
