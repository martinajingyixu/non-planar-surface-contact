%% Visualize the friction here
% surfaceType = 'unitSphereNew/';
% surfaceType = 'unitCylinderNew/';
% surfaceType = 'unitSphereNewVel/';
surfaceType = 'unitCylinderPartialHetztianLeft/';


pathProcessingSurface = ['frictionResults/' surfaceType];
load([pathProcessingSurface 'wrenchTwist101871.mat']);

%% Update samples
com = computedParams.com;

samplingRangeU = 0.06*pi;
% samplingRangeV = 0.4*pi;
samplingRangeV = 0.003;

lowerBoundU = computedParams.rangeU(1);
upperBoundU = computedParams.rangeU(2);

lowerBoundV = computedParams.rangeV(1);
upperBoundV = computedParams.rangeV(2);

samplesU = lowerBoundU:samplingRangeU:upperBoundU;
samplesV = lowerBoundV:samplingRangeV:upperBoundV;

%%

syms projectedVelocityDirection(a, b,c,x0,y0,z0,u,v) dFromCOR(a,b,c,x0,y0,z0,u,v)...
    frictionforcex(a,b,c,x0,y0,z0,u,v) frictionforcey(a,b,c,x0,y0,z0,u,v) ...
    frictionforcez(a,b,c,x0,y0,z0,u,v) VelocityDirection(a, b,c,x0,y0,z0,u,v)
projectedVelocityDirection(a, b,c,x0,y0,z0,u,v)...
    = computedParams.projectedVp.normalizedVp;
% projectedVelocityDirection(a, b,c,x0,y0,z0,u,v)...
%     = computedParams.projectedVp.normalizedVp;
VelocityDirection(a, b,c,x0,y0,z0,u,v) = computedParams.vp;
dFromCOR(a,b,c,x0,y0,z0,u,v) = computedParams.d;
frictionforcex(a,b,c,x0,y0,z0,u,v) = computedParams.fx2intSyms.general.correct;
frictionforcey(a,b,c,x0,y0,z0,u,v) = computedParams.fy2intSyms.general.correct;
frictionforcez(a,b,c,x0,y0,z0,u,v) = computedParams.fz2intSyms.general.correct;

%% visualize velocity and new COr sampling
iSampleV = 1% has problem
% iSampleV = 46029;
% iSampleV = 45970;
% iSampleV = 70301;
% iSampleV = 28375;

curSample = CORSamples4Integral(iSampleV,:);
aV= curSample(1);
bV = curSample(2);
cV = curSample(3);
x0V= curSample(4);
y0V= curSample(5);
z0V= curSample(6);
alphaSample = CORSamples4Analysis(iSampleV,1);
betaSample = CORSamples4Analysis(iSampleV,2);
distXSample = CORSamples4Analysis(iSampleV,3);
distYSample = CORSamples4Analysis(iSampleV,4);
distZSample = CORSamples4Analysis(iSampleV,5);

fx = wrench(:,1);
fy = wrench(:,2);
fz = wrench(:,3);
taux = wrench(:,4);
tauy = wrench(:,5);
tauz = wrench(:,6);

frictionForce = [fx(iSampleV),fy(iSampleV),fz(iSampleV)];
frictionTorque = [taux(iSampleV),tauy(iSampleV),tauz(iSampleV)];

figure
projectedVelocityArray = zeros(3,numel(samplesU),numel(samplesV));
velocityArray = zeros(3,numel(samplesU),numel(samplesV));

dArray = zeros(3,numel(samplesU),numel(samplesV));
frictionForceArray = zeros(3,numel(samplesU),numel(samplesV));


disp([ 'alpha: ' num2str(alphaSample)  ' beta: ' num2str(betaSample) ...
    ' distx: ' num2str(distXSample) ' disty: ' num2str(distYSample) ' distz: ' num2str(distZSample)]);
disp(['friction force: ' num2str(frictionForce.*1000)]);
disp(['friction tourque: ' num2str(frictionTorque.*1000)]);
% figure;

[FxDistrib,FyDistrib,FzDistrib,TauxDistrib,TauyDistrib,TauzDistrib]...
    =  deal(zeros(numel(samplesU)*numel(samplesV),3));

% iS = 1;
for iU = 1:2:numel(samplesU)
    for iV = 1:2:numel(samplesV)
        uSample = samplesU(iU);
        vSample = samplesV(iV);
      
        % plot forces direction
        projectedVelocityArray(:,iU,iV) = ...
    projectedVelocityDirection(aV,bV,cV,x0V,y0V,z0V,uSample,vSample);

velocityArray(:,iU,iV) = ...
        VelocityDirection(aV,bV,cV,x0V,y0V,z0V,uSample,vSample);


    dArray(:,iU,iV) = dFromCOR(aV,bV,cV,x0V,y0V,z0V,uSample,vSample);
    frictionForceArray(:,iU,iV) = ...
        [frictionforcex(aV,bV,cV,x0V,y0V,z0V,uSample,vSample);...
        frictionforcey(aV,bV,cV,x0V,y0V,z0V,uSample,vSample);...
        frictionforcez(aV,bV,cV,x0V,y0V,z0V,uSample,vSample)];

    pointOnSurface = parametric2cartesiansamples...
    (computedParams.rToSample,samplesU(iU),samplesV(iV));

    plot3(pointOnSurface(1),pointOnSurface(2),pointOnSurface(3),'.')
    normVel = norm(projectedVelocityArray(:,iU,iV));
   % PLOT the velocity

    P1 = pointOnSurface-0.05*projectedVelocityArray(:,iU,iV);
    P2= pointOnSurface+0.05*projectedVelocityArray(:,iU,iV);
    vectarrow(P1,P2,'b');
    hold on

    axis equal
       % PLOT d
    P3= pointOnSurface-dArray(:,iU,iV);
    P4 = pointOnSurface;

%     vectarrow(P3,P4,'-.og');
    hold on


    end
end
hold on
%%
% velocity of COM
% plotvector(twistProjected(iSampleV,1),twistProjected(iSampleV,2),twistProjected(iSampleV,3),...
%     com(1),com(2),com(3),'m')
hold on
plotvector(twistVelBasedProjected(iSampleV,1),twistVelBasedProjected(iSampleV,2),twistVelBasedProjected(iSampleV,3),...
    com(1),com(2),com(3),'k')
hold on
plotvector(twistVelBasedOriginal(iSampleV,1),twistVelBasedOriginal(iSampleV,2),twistVelBasedOriginal(iSampleV,3),...
    com(1),com(2),com(3),'m')
hold on
plotvector(frictionForce(1),frictionForce(2),frictionForce(3),...
    com(1),com(2),com(3),'k')
hold on
% COR
plotvector(aV,bV,cV,x0V,y0V,z0V,'c',10);

hold on
plotvector(twistAngularVelBasedProjected(iSampleV,4),twistAngularVelBasedProjected(iSampleV,5),...
    twistAngularVelBasedProjected(iSampleV,6),x0V,y0V,z0V,'r');

