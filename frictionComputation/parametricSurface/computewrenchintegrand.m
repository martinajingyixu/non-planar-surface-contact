function [surfaceInfo,frictionParam] = computewrenchintegrand(pathIntegrand,usePrecomputedVars)
%% define parametric form of the suface
syms u v r real
syms rToSample(u,v) 

%load('unitSphere/varsParametricForm.mat');
% path = ['frictionResults/'  surfaceType  '/'  pressureType  '/'];
load([pathIntegrand '/varsParametricForm.mat']);
%%
lowerBoundU = surfaceInfo.rangeU(1);
upperBoundU = surfaceInfo.rangeU(2);

lowerBoundV = surfaceInfo.rangeV(1);
upperBoundV = surfaceInfo.rangeV(2);

assume(u >= lowerBoundU);
assumeAlso(u <= upperBoundU);
assumeAlso(v >= lowerBoundV);
assumeAlso(v <= upperBoundV);

% dA = cos(theta); % later compute surface integral, dA = radius^2 * cos(theta)*dTheta*dPhi;
% pressureDistWeighted = 1;
dA = computedA(surfaceInfo.r,u,v);
% if not integrate friction center, compute it with discretized values
if ~usePrecomputedVars
    [frictionCenter,surfaceArea,normalForceSum] = ...
        frictioncenter(surfaceInfo.r,surfaceInfo.pressureDist,dA,u,v,surfaceInfo.rangeU,surfaceInfo.rangeV);
else

    surfaceArea = surfaceInfo.surfaceAreaInt;
    if strcmp(surfaceInfo.pressureType,'hetzCenter')
        normalForceSum = surfaceInfo.forceSumIntHetz;
        frictionCenter = surfaceInfo.pwcIntHetz;
    elseif strcmp(surfaceInfo.pressureType,'uniform')
        normalForceSum = surfaceInfo.forceSumIntUni;
        frictionCenter = surfaceInfo.pwcIntUni;
    end

end

% sampled coordinates for visualization
% weight the pressureDist such that the sum of the force is always 1N.
p0 = 1/double(normalForceSum);
pressureDistWeighted = surfaceInfo.pressureDist .* p0;
% the Nt normalized is the norm of the tangent plane
NtNormalized = computesurfacenormal(surfaceInfo.r,u,v);

%% compute the 6d wrench of normal pressure
if ~usePrecomputedVars

    [FIntegrand,tauIntegrand] = computenormalforcetorque(pressureDistWeighted,NtNormalized,dA,surfaceInfo.r,frictionCenter);
    normalForce = [integralcustimized(FIntegrand(1),lowerBoundU,upperBoundU,lowerBoundV,upperBoundV);...
        integralcustimized(FIntegrand(2),lowerBoundU,upperBoundU,lowerBoundV,upperBoundV);...
        integralcustimized(FIntegrand(3),lowerBoundU,upperBoundU,lowerBoundV,upperBoundV)];
    normalTorque = [integralcustimized(tauIntegrand(1),lowerBoundU,upperBoundU,lowerBoundV,upperBoundV);...
        integralcustimized(tauIntegrand(2),lowerBoundU,upperBoundU,lowerBoundV,upperBoundV);...
        integralcustimized(tauIntegrand(3),lowerBoundU,upperBoundU,lowerBoundV,upperBoundV)];
    normalWrench = [normalForce;normalTorque];
    pressureSum = integralcustimized(pressureDistWeighted,lowerBoundU,upperBoundU,lowerBoundV,upperBoundV);
    normalForceSumWeighted = integralcustimized(dA*pressureDistWeighted, lowerBoundU,upperBoundU,lowerBoundV,upperBoundV);
    assert(abs(normalForceSumWeighted-1)<1e-5,'The normal Force is scale to be 1');
    
end

%% ToDo: solve discretized problem later. plot the O, compute surface area.
%% define COR
% COR axis defined as: [x0;y0;z0] + [a;b;c]*t, where x0,y0,z0 are the point
% it goes through and a,b,c define the direction of COR
% assume the COR is the new x axis, distY distZ is the distance to the 
% COR, because sampling along the COR axis does not make sense.
syms a b c x0 y0 z0 real

assume(a>=0);
assumeAlso(a<=1);
assumeAlso(b>=0);
assumeAlso(b<=1);
assumeAlso(c>=0);
assumeAlso(c<=1);

%% compute d: the vector from the point P to the COR axis, d should be perpendicular to COR
%% Case 2: x0,y0,z0 not 0
% COR = [x0;y0;z0] + [a;b;c]*t
% syms d
% d = perpendicular2vec(r,a,b,c,x0,y0,z0);
% d = computeDpoint2COR(r,a,b,c,x0,y0,z0);
%% solve v: the vector perpendicular to both COR and d
syms vp h real
vp = computevelocity(surfaceInfo.r,a,b,c,x0,y0,z0,h);
%% compute the velocity projected to the tangent plane of P
projectedVp= projectedvelocity(vp,NtNormalized);

%% put all velocity related to one struct. later use for integral

%% General Form of integral
symbols.a = a;
symbols.b = b;
symbols.c = c;
symbols.x0 = x0;
symbols.y0 = y0;
symbols.z0 = z0;
symbols.u = u;
symbols.v = v;
symbols.h = h;

[fx2intSyms,fy2intSyms,fz2intSyms,taux2intSyms,tauy2intSyms,tauz2intSyms]...
    = friction2integrate(pressureDistWeighted,projectedVp,surfaceInfo.r,frictionCenter,dA,symbols);

%% Save info to struct
surfaceInfo.dA = dA;
surfaceInfo.frictionCenter = double(frictionCenter);
surfaceInfo.surfaceArea = surfaceArea;
surfaceInfo.pressureDistWeighted = pressureDistWeighted;
surfaceInfo.NtNormalized = NtNormalized;
if ~usePrecomputedVars
surfaceInfo.normalWrench = normalWrench; % O: friction center
surfaceInfo.normalForce = normalForce;
surfaceInfo.pressureSum = pressureSum;
surfaceInfo.normalForceSumWeighted = normalForceSumWeighted;
end

frictionParam.u = u;
frictionParam.v = v;
% frictionParam.d = d;
frictionParam.fx2intSyms = fx2intSyms;
frictionParam.fy2intSyms = fy2intSyms;
frictionParam.fz2intSyms = fz2intSyms;
frictionParam.taux2intSyms = taux2intSyms;
frictionParam.tauy2intSyms = tauy2intSyms;
frictionParam.tauz2intSyms = tauz2intSyms;
save([pathIntegrand '/FrictionParametricForm.mat'],'surfaceInfo','frictionParam');

end

