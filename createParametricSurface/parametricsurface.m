function [surfaceInfo] = parametricsurface(pathSurface,surfaceType,pressureType,rangeU,rangeV,r1,r2)
% Unit sphere, unit pressure distribution
syms u v r x y z ox oy oz real % u is Theta, v is Phi
syms rToSample(u,v) 
coordCenterInt = [];
surfaceAreaInt = [];
pwcIntHetz = [];
forceSumIntHetz = [];
forceSumIntUni = [];
pwcIntUni = [];
r3 = 0;

%% Default range of parameters of the contact area
if nargin< 7
    if strcmp(surfaceType,'cylinder')
        rangeU = [0,pi];
        rangeV = [0,1];% 
        r1 = 1;
        r2 = r1;
        r3 = 0;
        surfaceAreaInt = 3.14159;
        coordCenterInt = [0;0.63662;0.5];
%         forceSumIntHetz = 2.44664;
%         pwcIntHetz = [0;0.696282;0.5];
        forceSumIntHetz = 2.83047;
        pwcIntHetz = [0;0.660038;0.5];
        forceSumIntUni = surfaceAreaInt;
        pwcIntUni = coordCenterInt;
    elseif strcmp(surfaceType,'ellipticalCylinder')
        rangeU = [0,pi];
        rangeV = [0,1];
        r1 = 1;
        r2 = 0.5;
        r3 = 0;
        surfaceAreaInt = 2.42211;
        coordCenterInt = [0;0.352832;0.5];       
%         forceSumIntHetz = 1.94606;
%         pwcIntHetz = [0;0.382197;0.5];
        forceSumIntHetz = 2.24669;
        pwcIntHetz = [0;0.362045;0.5];
        forceSumIntUni = surfaceAreaInt;
        pwcIntUni = coordCenterInt;
    elseif strcmp(surfaceType,'sphere')
        rangeU = [-0.5*pi,0.5*pi];
        rangeV = [0,pi];
        r1 = 1;
        r2 = r1;
        r3 = 0;
        surfaceAreaInt = 6.28319;
        coordCenterInt = [0;0.5;0];  
%         forceSumIntHetz = 3.81203;
%         pwcIntHetz = [-0.00787777;0.58622;0.00696208];
        forceSumIntHetz = 5.49884;
        pwcIntHetz = [0;0.516813;0];
        forceSumIntUni = surfaceAreaInt;
        pwcIntUni = coordCenterInt;
    elseif strcmp(surfaceType,'ellipsoid')
        rangeU = [-0.5*pi,0.5*pi];
        rangeV = [0,pi];
        r1 = 1;
        r2 = 0.5*r1;
        r3 = 0.6*r1;
        surfaceAreaInt = 2.99715;
        coordCenterInt = [0;0.26757;0];  
%         forceSumIntHetz = 2.27279;
%         pwcIntHetz = [0.00539069;0.289907;0.000871128];
        forceSumIntHetz = 2.787;
        pwcIntHetz = [0;0.272562;0];
        forceSumIntUni = surfaceAreaInt;
        pwcIntUni = coordCenterInt;
    elseif strcmp(surfaceType,'paraboloid')
        rangeU = [0,pi];
        rangeV = [0,1];
        r1 = 1;
        r2 = r1; 
        r3 = 0;
        surfaceAreaInt = 2.66521;
        coordCenterInt = [0;0.455002;0.558937]; 
        forceSumIntHetz = 2.47883;
        pwcIntHetz = [0;0.46087;0.551788];
        forceSumIntUni = surfaceAreaInt;
        pwcIntUni = coordCenterInt;
    elseif strcmp(surfaceType,'ellipticalParaboloid')
        rangeU = [0,pi];
        rangeV = [0,1];
        r1 = 1;
        r2 = 0.5;
        r3 = 0;
        surfaceAreaInt = 1.82807;
        coordCenterInt = [0;0.252747;0.575772]; 
        forceSumIntHetz = 1.72721;
        pwcIntHetz = [0;0.25583;0.571008];
        forceSumIntUni = surfaceAreaInt;
        pwcIntUni = coordCenterInt;
    else
        error('can not recognize the type of the parametric surface!');
    end
end


if strcmp(surfaceType,'cylinder')
    r  = [r1*cos(u); r1*sin(u); v];
elseif strcmp(surfaceType,'sphere')
    r  = [r1*cos(u)*cos(v); r1*cos(u)*sin(v); r1*sin(u)];
elseif strcmp(surfaceType,'ellipticalCylinder')
    r  = [r1*cos(u); r2*sin(u); v];
elseif strcmp(surfaceType,'ellipticalParaboloid')
    r  = [r1*cos(u)*v; r2*sin(u)*v;v*v];    
elseif strcmp(surfaceType,'paraboloid')
    r  = [r1*cos(u)*v; r1*sin(u)*v;v*v];        
elseif strcmp(surfaceType,'ellipsoid')
    r  = [r1*cos(u)*cos(v); r2*cos(u)*sin(v); r3*sin(u)];
else
    error('can not recognize the type of the parametric surface!');
end

lowerBoundU = rangeU(1);
upperBoundU = rangeU(2);

lowerBoundV = rangeV(1);
upperBoundV = rangeV(2);

%% provide an estimated coord center. only use for integrand computation when integral can not be solved.
rToSample(u,v) = r;

surfaceInfo.lowerBoundU = lowerBoundU;
surfaceInfo.upperBoundU = upperBoundU;
surfaceInfo.lowerBoundV = lowerBoundV;
surfaceInfo.upperBoundV = upperBoundV;
surfaceInfo.rangeU = rangeU;
surfaceInfo.rangeV = rangeV;
surfaceInfo.meanRadius = mean([r1,r2,r3]);
surfaceInfo.r = r;
surfaceInfo.rToSample = rToSample;

surfaceInfo.coordCenterInt = coordCenterInt; % the pre-integrated com in mathematica
surfaceInfo.surfaceAreaInt = surfaceAreaInt;
surfaceInfo.pressureType = pressureType;
surfaceInfo.surfaceType = surfaceType;

surfaceInfo.forceSumIntHetz = forceSumIntHetz;
surfaceInfo.forceSumIntUni = forceSumIntUni;
surfaceInfo.pwcIntUni = pwcIntUni;
surfaceInfo.pwcIntHetz = pwcIntHetz;
surfaceInfo.r1 = r1;
surfaceInfo.r2 = r2;
surfaceInfo.r3 = r3;

%% store variables 
if ~exist(pathSurface, 'dir')
   mkdir(pathSurface);
end
save([pathSurface  'varsParametricForm.mat'],'surfaceInfo');

end

