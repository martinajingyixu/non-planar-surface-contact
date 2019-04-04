function [fittedRes4th,fittedResEllip] = fitlimitsurfaces(pathLS,config)

lambda = config.lambda;
wrenchWeight = config.wrenchWeight;
twistWeight=config.twistWeight;
ifConvex=config.ifConvex;
ifVisualize=config.ifVisualize;
ifSampling=config.ifSampling;
samplingStepSize=config.samplingStepSize;

load([pathLS '/wrenchTwist.mat']);

if ~exist('twist2fit','var')
    twist2fit = zeros(size(wrench2fit));
end

%% fit 4th-order poly model
if config.ifFitLS4th == 1
    fittedRes4th = ...
        fit4thorderpoly6d(wrench2fit, twist2fit, lambda, twistWeight, wrenchWeight,ifConvex);
    if ifSampling
        [sampledWrench4th,shapeArray4th]=...
            samplefittedlimitsurface(wrench2fit,fittedRes4th.fittedCoeff,samplingStepSize,4);
        fittedRes4th.sampledWrench = sampledWrench4th;
        fittedRes4th.shapeArray = shapeArray4th;
    end
    if ifVisualize
        visualizedfitted6dgeometry(wrench2fit,sampledWrench4th,twist2fit);           
    end
else
    fittedRes4th = [];
end

%% fit ellipsoidal model
if  config.ifFitLSEllip == 1
    fittedResEllip = ...
        fitellipsoid6d(wrench2fit,twist2fit,lambda,wrenchWeight,twistWeight);
    if ifSampling 
        numPoints = 1;
        [~,sampledWrenchEllip] = ...
            sample6dellipsoid(fittedResEllip.fittedCoeff,numPoints,false,ifVisualize);

        sampledWrenchEllip = sampledWrenchEllip';    
        sampledWrenchEllip= sampledWrenchEllip...
                ( ~any( isnan( sampledWrenchEllip ) | isinf( sampledWrenchEllip ), 2 ),: );
        sampledWrenchEllip = sampledWrenchEllip';
        fittedResEllip.sampledWrench = sampledWrenchEllip;
    end
    if ifVisualize
        visualizedfitted6dgeometry(wrench2fit,sampledWrenchEllip,twist2fit);           
    end
else
    fittedResEllip = [];
end

%%
save([pathLS 'fittedLS.mat'], 'fittedRes4th', 'fittedResEllip','-v7.3');

end

