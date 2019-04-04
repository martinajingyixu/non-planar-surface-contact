%% 1. Create contact with parametric surface with different pressure distribution, compute wrench and fit limit surfaces
surfaceTypeArray = {'cylinder','ellipticalCylinder','sphere','ellipsoid','paraboloid','ellipticalParaboloid'};
pressureTypeArray = {'uniform','hetzCenter'};
%% settings
config.ifVisualize = 0;
config.ifSampling = 1;
config.ifConvex = 1;
config.lambda = 1;
config.wrenchWeight = 10;
config.twistWeight = 1;
config.samplingStepSize = 5;
config.fineLevel = 0;
config.ifFitLS4th = 1;
config.ifFitLSEllip = 1;
config.fineLevel = -1; % different level of IAR densities

usePrecomputedVars = 1; % integrated friction center for parametric surfaces
%% compute friction wrench for each parametric surface with the two pressure distributions
for iP = 1:2
    for iS = 1:6
        iS
%         try
        surfaceType = surfaceTypeArray{iS};
        pressureType = pressureTypeArray{iP};
        pathResults = ['frictionResults/parametricSurface/'  surfaceType  '/'  pressureType  '/'];

        % create parametric surface and compute integrands for wrench and twists
        surfaceInfo = parametricsurface(pathResults,surfaceType,pressureType);
        surfaceInfo = definepressureparametricsurface(pathResults,surfaceInfo,config.ifVisualize);
        [surfaceInfo,frictionParam] = computewrenchintegrand(pathResults,usePrecomputedVars);
        % sample Twist
        twistSampleVariables = ...
            generateTwistSamples(surfaceInfo.frictionCenter,surfaceInfo.lengthArray,config.fineLevel,pathResults);
        % integrate wrench
        [wrench,twist] = integratewrench(pathResults,config.ifVisualize);
        % fit LS
         [fittedRes4th,fittedResEllip] = fitlimitsurfaces(pathResults,config);            
    end
end
%% post-processing and visualization
%% fitting error
% compute normalzed wrench twist
[wrenchNormalized,twistNormalized] = normalizewrenchtwist(wrench,twist);
[meanWrenchError4th,meanTwistError4th] = ...
    computefittingerror4th(fittedRes4th.fittedCoeff,wrenchNormalized,twistNormalized);
[meanWrenchErrorEllip,meanTwistErrorEllip] = ...
    computefittingerrorellip(fittedResEllip.fittedCoeff,wrenchNormalized,twistNormalized);
%% visualize computed friction wrench
plot6dwrenches(wrenchNormalized)
%% visualize LS cross-sections
% use stepSize to adjust number of samples to be visualized. densely down-sample the wrenches.
stepSize = 10;
[~,idxArray2Plot] = downsamplewrench(wrenchNormalized,stepSize);
wrench2plotNormalized = wrenchNormalized(:,idxArray2Plot);
twist2plotNormalized = twistNormalized(:,idxArray2Plot);
visualizeLScrosssection(wrench2plotNormalized,twist2plotNormalized,fittedResEllip.fittedCoeff,2)
visualizeLScrosssection(wrench2plotNormalized,twist2plotNormalized,fittedRes4th.fittedCoeff,4)
