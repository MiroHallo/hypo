function [AI,TT] = ray1d(Depth,Vel,DepthS,DepthR,DistX,DX)

frac_ang = 4; % Lowering factor of the angle step search
MaxLoop = int16(100); % Maximal number of loop cycles for searching of the ray
irt_N = 4; % Number of Initial ray tracing angles (irt_step*(irt_N-1) must give 90)
irt_step = 30; % Step of Initial ray tracing [deg]

%% Allocation
nLayers = int16(length(Depth));
nDist = int16(length(DistX));
AI = zeros(nDist,1);
TT = zeros(nDist,1) + Inf;
IS = int16(find(Depth >= DepthS,1,'first'));
if isempty(IS)
    IS = nLayers;
end
IR = int16(find(Depth >= DepthR,1,'first'));
if isempty(IR)
    IR = nLayers;
end
Thick = zeros(nLayers,1);
Thick(1,1) = Inf;
Thick(2:end,1) = Depth(2:end) - Depth(1:end-1);
if IS == 1
    ISresDepU = Inf;
else
    ISresDepU =  DepthS - Depth(IS-1);
end
if IR == 1
    IRresDepU = Inf;
else
    IRresDepU = DepthR - Depth(IR-1);
end
ISresDepD = Depth(IS) - DepthS;
IRresDepD = Depth(IR) - DepthR;
DepDel =  DepthS - DepthR;

%% Horizontal ray
if DepDel == 0
    TT = DistX./Vel(IS);
    AI = AI + 90;
    return
end

%% Initial ray tracing
if DepDel > 0
    irtA = (0:irt_N-1)*irt_step;
else
    irtA = 90 + (0:irt_N-1)*irt_step;
end
irtX = zeros(1,irt_N);
for i = 1 : irt_N
    [~,irtX(i),~] = TX(Vel,Thick,irtA(i),DepthS,DepthR,IS,IR,ISresDepD,ISresDepU,IRresDepD,IRresDepU);
end

%% Finishing ray tracing
for d = int16(1) : nDist
    di = DistX(d);
    if abs(irtX(1) - di) < DX
        AI(d) = irtA(1);
        [TT(d),~,~] = TX(Vel,Thick,AI(d),DepthS,DepthR,IS,IR,ISresDepD,ISresDepU,IRresDepD,IRresDepU);
        continue
    elseif irtX(1) < di
        grad = 1;
        chd = int16(1);
        AngI = find(irtX > di,1,'first');
        A = irtA(AngI-1);
    else
        grad = -1;
        chd = int16(1);
        AngI = find(irtX <= di,1,'first');
        A = irtA(AngI);
    end
    AstepA = irt_step/frac_ang;
    count = int16(0);
    
    while count <= MaxLoop
        count = count + int16(1);
        [Tn,Xn,~] = TX(Vel,Thick,A,DepthS,DepthR,IS,IR,ISresDepD,ISresDepU,IRresDepD,IRresDepU);
        if abs(Xn - di) < DX
            break
        else
            if chd == int16(1)
                if Xn < di
                    A = A + grad*AstepA;
                else
                    AstepA = AstepA/frac_ang;
                    chd = int16(-1);
                    A = A - grad*AstepA;
                end
                
            else
                if Xn > di
                    A = A - grad*AstepA;
                else
                    AstepA = AstepA/frac_ang;
                    chd = int16(1);
                    A = A + grad*AstepA;
                end
            end
        end
    end
    AI(d) = A;
    TT(d) = Tn;
end
return