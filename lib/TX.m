function [T,X,A] = TX(Vel,Thick,A,DepthS,DepthR,IS,IR,ISresDepD,ISresDepU,IRresDepD,IRresDepU)
% XT Ray tracer in the layered model.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Ray tracer in the layered velocity model
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
%        Vel = wave velocities [m/s]
%        Thick = layer thicknesses [meters]
%        A  = Take-off dips from vertical (0-180 [deg], where 0=UP, 180=DOWN)
%        DepthS = depth of the source [meters] (depth positive down)
%        DepthR = depth of the receiver [meters] (depth positive down)
%        IS = index of source Layer
%        IR = index of receiver Layer
%        ISresDepD = vertical dist. to boundary at source layer [meters]
%        ISresDepU = vertical dist. to boundary at source layer [meters]
%        IRresDepD = vertical dist. to boundary at receiver layer [meters]
%        IRresDepU = vertical dist. to boundary at receiver layer [meters]
%
% OUTPUT:
%        T = Travel time from the source to the depth of the receiver [seconds]
%        X = Horizontal distance between the source and position where the ray 
%            intercept the depth of receiver [meters]
%        A = Receiver incidence angle (0-180 [deg], where 0=UP, 180=DOWN)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Allocation
DepDel =  DepthS - DepthR; % delDepths (positive for deaper source)
p = sind(A) / Vel(IS); % ray parameter

% Horizontal ray
if (A == 90) %% Horizontal ray to Inf
    T = Inf;
    X = Inf;
    return

% Normal ray
else 
    
    % Invalid direction
    if (A < 90) && (DepDel <= 0)
        T = Inf;
        X = Inf;
        return
    end
    if (A > 90) && (DepDel >= 0)
        T = Inf;
        X = Inf;
        return
    end
    
    % ---------------------------------
    % First layer
    
    if (IS == IR) % Source and receive in the same layer (1st layer)
        X = abs(tand(A)*DepDel);
        T = sqrt(DepDel^2 + X^2) / Vel(IS);
        return
        
    else % Source and receive in diferent layers (1st layer)
        
        if A < 90 % Upgoing
            
            % Next layer
            IN = IS - 1;
            
            X = tand(A)*ISresDepU;
            T = sqrt(ISresDepU^2 + X^2) / Vel(IS);
            A = asind(p*Vel(IN));
            
        else  % Downgoing
            
            % Next layer
            IN = IS + 1;
            
            X = abs(tand(A))*ISresDepD;
            T = sqrt(ISresDepD^2 + X^2) / Vel(IS);
            A = 180 - asind(p*Vel(IN));
        end        
        
        % Return if a reflection
        if ~isreal(A) % Complex angle == reflection
            T = Inf;
            X = Inf;
            A = 0;
            return
        end
    end
    
    % ---------------------------------
    % Nth layer

    INi = IN;
    
    while INi ~= IR
        
        if A < 90 % Upgoing
            
            % Next layer
            IN = INi - 1;
            
            Xn = tand(A)*Thick(INi);
            X = X + Xn;
            T = T + (sqrt(Thick(INi)^2 + Xn^2) / Vel(INi));
            A = asind(p*Vel(IN));
            
        else  % Downgoing
            
            % Next layer
            IN = INi + 1;
            
            Xn = abs(tand(A))*Thick(INi);
            X = X + Xn;
            T = T + sqrt(Thick(INi)^2 + Xn^2) / Vel(INi);
            A = 180 - asind(p*Vel(IN));
            
        end
        
        % Return if a reflection
        if ~isreal(A) % Complex angle == reflection
            T = Inf;
            X = Inf;
            A = 0;
            return
        end
        
        INi = IN;
    end
    
    % ---------------------------------
    % Last layer
    if A < 90 % Upgoing
        
        Xn = tand(A)*IRresDepD;
        X = X + Xn;
        T = T + sqrt(IRresDepD^2 + Xn^2) / Vel(IR);

    else  % Downgoing
        
        Xn = abs(tand(A))*IRresDepU;
        X = X + Xn;
        T = T + sqrt(IRresDepU^2 + Xn^2) / Vel(IR);
        
    end
    
end

end
