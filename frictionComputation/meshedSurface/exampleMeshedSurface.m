%% 1.parametric surface with different pressure distribution
surfaceTypeArray = {'cylinder','ellipticalCylinder','sphere','ellipsoid','paraboloid','ellipticalParaboloid'};
pressureTypeArray = {'uniform','hetzCenter'};
% current surface and pressure type
iS = 4; % ellipsoid
iP = 1; % uniform
%% 
config.ifNewFigure = true;
config.ifSampling = 1;
config.ifConvex = 1;
config.lambda = 1;
config.wrenchWeight = 10;
config.twistWeight = 1;
config.samplingStepSize = 6;
config.fineLevel = 1;
config.elemSz = 0.2;
config.ifFitLS4th = 1;
config.ifFitLSEllip = 1;
%%
surfaceType = surfaceTypeArray{iS};
pressureType = pressureTypeArray{iP};
pathResults = ['frictionResults/discretizedSurface/'  surfaceType  '/'  pressureType  '/'];

surfaceInfo = parametricsurface(pathResults,surfaceType,pressureType);
%     descritize the elements. Use center coord of the nodes to determine the pressure of the current element. Then scale the sum of the force to one
[elementsInfo,surfaceInfo] = getdiscretizationresults...
    (surfaceInfo,pathResults,config.ifVisualize,config.elemSz);
surfaceInfo.elementsInfo = elementsInfo;

%%
[fittingResults,wrenchTwist] = ...
    computewrenchtwistfitls(elementsInfo,pathResults,config);
displayfittingerror(fittingResults.fittingErr);

save([pathResults 'frictionContactResults.mat'],...
    'elementsInfo','surfaceInfo','fittingResults','wrenchTwist');

%%
wrenchNormalized = wrenchTwist.wrenchNormalized;
twistNormalized =  wrenchTwist.twistNormalized;
%% visualize computed friction wrench
if config.ifVisualize
    plot6dwrenches(wrenchTwist.wrenchNormalized)
end
%% visualize LS 3D cross-sections
% use stepSize to adjust number of samples to be visualized. densely down-sample the wrenches.
% stepSize = 10;
% [~,idxArray2Plot] = downsamplewrench(wrenchNormalized,stepSize);
% wrench2plotNormalized = wrenchNormalized(:,idxArray2Plot);
% twist2plotNormalized = twistNormalized(:,idxArray2Plot);
% visualizeLScrosssection(wrench2plotNormalized,twist2plotNormalized,fittingResults.fittedResEllip.fittedCoeff,2)
% visualizeLScrosssection(wrench2plotNormalized,twist2plotNormalized,fittingResults.fittedRes4th.fittedCoeff,4)
