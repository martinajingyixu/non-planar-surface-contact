%% Visualize the friction here
% pathProcessingSurface = 'unitSphereSpecialCase';
%pathProcessingSurface = 'unitCylinderThin';
pathProcessingSurface = 'unitCylinderThinQuater';

load([pathProcessingSurface '/FrictionParametricForm.mat']);
dataCOR = load([pathProcessingSurface '/CORSamples.mat']);

%% Update samples
rangeU = [0,0.5*pi];

% samplingRangeU = 0.1*pi;
% samplingRangeV = 0.2*pi;

samplingRangeU = 0.05*pi;

lowerBoundU = rangeU(1);
upperBoundU = rangeU(2);



samplesU = lowerBoundU:samplingRangeU:upperBoundU;
%%

syms projectedVelocityDirection(a, b,c,x0,y0,z0,u)
projectedVelocityDirection(a, b,c,x0,y0,z0,u)...
    = projectedVpNormalized.correct;

syms projectedVelocityDirectionLB(a, b,c,x0,y0,z0,u)
projectedVelocityDirectionLB(a, b,c,x0,y0,z0,u)...
    = projectedVpNormalized.lowerBoundAbs;
iSampleV = 1;
%%
CORSamples4Analysis = dataCOR.CORSamples4Analysis;
CORSamples4Integral = dataCOR.CORSamples4Integral;

COR4Subs = CORSamples4Integral;
COR4Analyse = CORSamples4Analysis;
fx = integratedFriction.fxArray;
fy = integratedFriction.fyArray;
fz = integratedFriction.fzArray;
taux = integratedFriction.tauxArray;
tauy = integratedFriction.tauyArray;
tauz = integratedFriction.tauzArray;
%%
idx = find(abs(CORSamples4Integral(:,1))<1e-2 &...
CORSamples4Integral(:,4) ==0 & abs(CORSamples4Integral(:,6))<1e-1)
%%
idx = find(abs(CORSamples4Integral(:,1))<1e-2&...
abs(CORSamples4Integral(:,4))<1e-2&abs(CORSamples4Integral(:,6)-1)<1e-2)
%%
idx = find(abs(CORSamples4Integral(:,1))<1e-2)
%% visualize velocity and new COr sampling
iTest = 1;
%%
close all
curSample = COR4Subs(iSampleV,:);
aV= curSample(1);
bV = curSample(2);
cV = curSample(3);
x0V= curSample(4);
y0V= curSample(5);
z0V= curSample(6);

alphaSample = COR4Analyse(iSampleV,1);
betaSample = COR4Analyse(iSampleV,2);
distXSample = COR4Analyse(iSampleV,3);
distYSample = COR4Analyse(iSampleV,4);
distZSample = COR4Analyse(iSampleV,5);


frictionForce = [fx(iSampleV),fy(iSampleV),fz(iSampleV)];
frictionTourque = [taux(iSampleV),tauy(iSampleV),tauz(iSampleV)];



projectedVelocityArray = zeros(3,numel(samplesU));
projectedVelocityLBArray = zeros(3,numel(samplesU));



disp([ 'alpha: ' num2str(alphaSample)  ' beta: ' num2str(betaSample) ...
    ' distx: ' num2str(distXSample) ' disty: ' num2str(distYSample) ' distz: ' num2str(distZSample)]);

disp([ 'a: ' num2str(aV)  ' b: ' num2str(bV) ...
    ' c: ' num2str(cV) ]);
disp([ 'x0: ' num2str(x0V)  ' y0: ' num2str(y0V) ...
    ' z0: ' num2str(z0V) ]);

disp(['friction force: ' num2str(frictionForce)]);
disp(['friction tourque: ' num2str(frictionTourque)]);
% figure;

[FxDistrib,FyDistrib,FzDistrib,TauxDistrib,TauyDistrib,TauzDistrib]...
    =  deal(ones(numel(samplesU)*numel(samplesV),3));

iS = 1;
for iU = 1:numel(samplesU)
        uSample = samplesU(iU);
      
        % plot forces direction
        projectedVelocityArray(:,iU) = ...
    vpa(projectedVelocityDirection(aV,bV,cV,x0V,y0V,z0V,uSample));
        projectedVelocityLBArray(:,iU) = ...
    vpa(projectedVelocityDirectionLB(aV,bV,cV,x0V,y0V,z0V,uSample));


    pointOnSurface = rToSample(uSample,1);

    plot3(pointOnSurface(1),pointOnSurface(2),pointOnSurface(3),'.')

    P1 = pointOnSurface-0.5*projectedVelocityArray(:,iU);
    P2= pointOnSurface+0.5*projectedVelocityArray(:,iU);
    vectarrow(P1,P2,'b');
        hold on


    P3 = pointOnSurface-0.5*projectedVelocityLBArray(:,iU);
    P4= pointOnSurface+0.5*projectedVelocityLBArray(:,iU);
    vectarrow(P3,P4,'g');
    hold on
    axis equal
    
    
%     
%     FxDistrib(iS,1:2) = [uSample,fx2intSyms.general(aV,bV,cV,x0V,y0V,z0V,uSample)];
%     FyDistrib(iS,1:2) = [uSample,fy2intSyms.general(aV,bV,cV,x0V,y0V,z0V,uSample)];
%     FzDistrib(iS,1:2) = [uSample,fz2intSyms.general(aV,bV,cV,x0V,y0V,z0V,uSample)];
% 
%     TauxDistrib(iS,1:2) = [uSample,taux2intSyms.general(aV,bV,cV,x0V,y0V,z0V,uSample)];
%     TauyDistrib(iS,1:2) = [uSample,tauy2intSyms.general(aV,bV,cV,x0V,y0V,z0V,uSample)];
%     TauzDistrib(iS,1:2) = [uSample,tauz2intSyms.general(aV,bV,cV,x0V,y0V,z0V,uSample)];
%     iS = iS + 1;
    
end

hold on
%plotCOR(aV,bV,cV,x0V,y0V,z0V)
iSampleV = iSampleV+1;

%% plot u,v,fx before integral and fit the surface
plot3Dsurface(FxDistrib,FyDistrib,FzDistrib,...
    TauxDistrib,TauyDistrib,TauzDistrib);
%%
xdata = FxDistrib(:,1);
ydata = FxDistrib(:,2);
zdata = FxDistrib(:,3);
%%
fxFit = fit([FxDistrib(:,1),FxDistrib(:,2)],FxDistrib(:,3),'linearinterp');

plot(fxFit,[FxDistrib(:,1),FxDistrib(:,2)],FxDistrib(:,3))
%%
testValue = 0:0.01*pi:pi;
sinValues = sin(testValue);
syms x approc(x)
approc(x) = (16*x*(pi-x))/(5*pi*pi-4*x*(pi-x));
approcValues = double(vpa(approc(testValue)));
diff = sinValues - approcValues;
values = [sinValues',approcValues',diff'];
%%
% ui = linspace(0,0.5*pi,100);
% vi = linspace(0,2*pi-0.00001,100);
% 
% test = gridfit(FxDistrib(:,1),FxDistrib(:,2),FxDistrib(:,3),ui,vi);
% figure
% surf(ui,vi,test)

%%
syms aa bb cc
% only when a*a > b*b + c*c
fsqrt(aa,bb,cc) = sqrt(aa*aa+bb*bb+cc*cc);
fsqrtapproc(aa,bb,cc) = 0.96*aa+0.384*bb + 0.16*cc;
aaV = 0.2:0.52:10.2;
bbV = 0.1:0.51:9.9;
ccV = 0:0.5:9.8;
aaV = aaV*10;
sqrtValues = fsqrt(aaV,bbV,ccV);
fsqrtapprocValue = fsqrtapproc(aaV,bbV,ccV);

diff = double(sqrtValues)'-double(fsqrtapprocValue)';
values = [double(sqrtValues)',double(fsqrtapprocValue)',diff];
%% Old: test COR coord transform
angleY = alphaSample;
angleZ = betaSample;

Ry = rotateY(angleY);
Rz = rotateZ(angleZ);
rotationMatrixS = Ry*Rz; % first rotate around z then y, counter clock wise.

oldX = (rotationMatrixS)\[1;0;0];
oldY = (rotationMatrixS)\[0;1;0];
oldZ = (rotationMatrixS)\[0;0;1];

newCOR = [0;1;-1];

% newCOR = [0;distYSample;distZSample];
COR = (rotationMatrixS)\newCOR;
CORSym = simplify((rotationMatrix)\[0;distY;distZ]);
distO2COR(alpha,beta,distY,distZ)= frictionCenter-CORSym;
distO2COR(alphaSample,betaSample,distYSample,distZSample)

origin = [0;0;0];

P1CORGlobal = origin + oldX;
vectarrow(origin,P1CORGlobal,'r');
hold on

P1CORGlobal = origin + oldY;
vectarrow(origin,P1CORGlobal,'g');
hold on

P1CORGlobal = origin + oldZ;
vectarrow(origin,P1CORGlobal,'b');
hold on

P1CORGlobal = origin + [1;0;0];
vectarrow(origin,P1CORGlobal,'-.r');
hold on

P1CORGlobal = origin + [0;1;0];
vectarrow(origin,P1CORGlobal,'-.g');
hold on

P1CORGlobal = origin + [0;0;1];
vectarrow(origin,P1CORGlobal,'-.b');
hold on

plot3(COR(1),COR(2),COR(3),'bo','MarkerSize',30)

% plot3(newCOR(1),newCOR(2),newCOR(3),'ro','MarkerSize',30)

axis equal
