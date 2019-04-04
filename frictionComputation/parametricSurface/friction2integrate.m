% function [fx2intSyms,fy2intSyms,fz2intSyms,taux2intSyms,tauy2intSyms,tauz2intSyms,...
%     fxMinus2intSyms,fyMinus2intSyms,fzMinus2intSyms,tauxMinus2intSyms,tauyMinus2intSyms,tauzMinus2intSyms]...
%     = friction2integrate(projectedVp,r,origin,dA,symbols)
function [fx2intSyms,fy2intSyms,fz2intSyms,taux2intSyms,tauy2intSyms,tauz2intSyms]...
    = friction2integrate(pressureDist,projectedVp,r,origin,dA,symbols)

a = symbols.a;
b = symbols.b;
c = symbols.c;
x0 = symbols.x0;
y0 = symbols.y0;
z0 = symbols.z0;
u = symbols.u;
v = symbols.v;
h = symbols.h;

%% Compute friction based on correct form

[fxSyms,fySyms,fzSyms,tauxSyms,tauySyms,tauzSyms]...
    = computefrictionVppNormbased(pressureDist,projectedVp.normalizedVp,r,origin,dA);

% 
% fx2intSyms(a,b,c,x0,y0,z0,h,u,v) = simplify(fxSyms);
% fy2intSyms(a,b,c,x0,y0,z0,h,u,v) = simplify(fySyms);
% fz2intSyms(a,b,c,x0,y0,z0,h,u,v) = simplify(fzSyms);
% 
% taux2intSyms(a,b,c,x0,y0,z0,h,u,v) = simplify(tauxSyms);
% tauy2intSyms(a,b,c,x0,y0,z0,h,u,v) = simplify(tauySyms);
% tauz2intSyms(a,b,c,x0,y0,z0,h,u,v) = simplify(tauzSyms);


fx2intSyms(a,b,c,x0,y0,z0,h,u,v) = fxSyms;
fy2intSyms(a,b,c,x0,y0,z0,h,u,v) = fySyms;
fz2intSyms(a,b,c,x0,y0,z0,h,u,v) = fzSyms;

taux2intSyms(a,b,c,x0,y0,z0,h,u,v) = tauxSyms;
tauy2intSyms(a,b,c,x0,y0,z0,h,u,v) = tauySyms;
tauz2intSyms(a,b,c,x0,y0,z0,h,u,v) = tauzSyms;



end

