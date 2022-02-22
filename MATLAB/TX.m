function [T,X,A] = TX(Vel,Thick,A,DepthS,DepthR,IS,IR,ISresDepD,ISresDepU,IRresDepD,IRresDepU)
DepDel =  DepthS - DepthR;
p = sind(A) / Vel(IS);
if (A == 90)
    T = Inf;
    X = Inf;
    return
else 
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
    if (IS == IR)
        X = abs(tand(A)*DepDel);
        T = sqrt(DepDel^2 + X^2) / Vel(IS);
        return
    else
        if A < 90
            IN = IS - 1;
            X = tand(A)*ISresDepU;
            T = sqrt(ISresDepU^2 + X^2) / Vel(IS);
            A = asind(p*Vel(IN));
        else
            IN = IS + 1;
            X = abs(tand(A))*ISresDepD;
            T = sqrt(ISresDepD^2 + X^2) / Vel(IS);
            A = 180 - asind(p*Vel(IN));
        end
        if ~isreal(A)
            T = Inf;
            X = Inf;
            A = 0;
            return
        end
    end
    INi = IN;
    while INi ~= IR
        if A < 90
            IN = INi - 1;
            Xn = tand(A)*Thick(INi);
            X = X + Xn;
            T = T + (sqrt(Thick(INi)^2 + Xn^2) / Vel(INi));
            A = asind(p*Vel(IN));
        else
            IN = INi + 1;
            Xn = abs(tand(A))*Thick(INi);
            X = X + Xn;
            T = T + sqrt(Thick(INi)^2 + Xn^2) / Vel(INi);
            A = 180 - asind(p*Vel(IN));
        end
        if ~isreal(A)
            T = Inf;
            X = Inf;
            A = 0;
            return
        end
        INi = IN;
    end
    if A < 90
        Xn = tand(A)*IRresDepD;
        X = X + Xn;
        T = T + sqrt(IRresDepD^2 + Xn^2) / Vel(IR);
    else
        Xn = abs(tand(A))*IRresDepU;
        X = X + Xn;
        T = T + sqrt(IRresDepU^2 + Xn^2) / Vel(IR);
    end
end
return