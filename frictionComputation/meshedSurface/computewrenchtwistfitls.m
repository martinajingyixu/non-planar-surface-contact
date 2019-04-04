function [fittingResults,wrenchTwist] = computewrenchtwistfitls(contactInfo,pathResults,config)
%% given contactInfo (surface, pressure, etc), compute wrench twist samples then fit LS models
if ~isfield(config,'ifFitLS4th')
    ifFitLS4th = 0;
end
if ~isfield(config,'ifFitLSEllip')
    ifFitLSEllip = 0;
end
twistSampleVariables = generateTwistSamples(contactInfo.pwc,contactInfo.lengthArray,config.fineLevel);
% tic
wrench = computewrenchdiscrete(contactInfo,twistSampleVariables,config.ifVisualize);
% wrenchComputationTime = toc;
[row, col] = find(isnan(wrench));
wrench(:,unique(col)) = [];
twistSampleVariables(:,unique(col)) = [];
twist = computetwist(contactInfo.pwc,twistSampleVariables);
[wrenchNormalized,twistNormalized] = normalizewrenchtwist(wrench,twist);
    %%
[~,idxArray] = downsamplewrench(wrenchNormalized);
wrench2fit = wrenchNormalized(:,idxArray);
twist2fit = twistNormalized(:,idxArray);
twist2fit = bsxfun(@rdivide, twist2fit, sqrt(sum(twist2fit.^2)));
if ~exist(pathResults,'dir')
    mkdir(pathResults);
end
save([pathResults, 'wrenchTwist.mat'],'twistSampleVariables','wrench','twist','idxArray',...
    'wrench2fit','twist2fit');
wrenchTwist.wrench = wrench;
wrenchTwist.twist = twist;
wrenchTwist.idxArrayForLSFitting = idxArray;
wrenchTwist.wrenchNormalized = wrenchNormalized;
wrenchTwist.twistNormalized = twistNormalized;

%% fit
fittingResults = [];
fittingError = [];
if config.ifFitLS4th == 0 && config.ifFitLSEllip == 0
    fittingResults = [];
else
    [fittedRes4th,fittedResEllip] = fitlimitsurfaces(pathResults,config);
    if config.ifFitLS4th == 1
        [meanWrenchError4th,meanTwistError4th] = ...
        computefittingerror4th(fittedRes4th.fittedCoeff,wrenchNormalized,twistNormalized);
        fittingError.meanWrenchError4th = meanWrenchError4th;
        fittingError.meanTwistError4th = meanTwistError4th;
        fittingResults.fittedRes4th = fittedRes4th;

    end
    if config.ifFitLSEllip == 1
        [meanWrenchErrorEllip,meanTwistErrorEllip] = ...
        computefittingerrorellip(fittedResEllip.fittedCoeff,wrenchNormalized,twistNormalized);
        fittingError.meanWrenchErrorEllip = meanWrenchErrorEllip;
        fittingError.meanTwistErrorEllip = meanTwistErrorEllip;
        fittingResults.fittedResEllip = fittedResEllip;
    end
    fittingResults.fittingErr = fittingError;
    displayfittingerror(fittingError)
end

%%

if config.ifVisualize == 1
    stepSize = 8;
    [~,idxArray2Plot] = downsamplewrench(wrenchNormalized,stepSize);
    wrench2plotNormalized = wrenchNormalized(:,idxArray2Plot);
    twist2plotNormalized = twistNormalized(:,idxArray2Plot);
    if config.ifFitLS4th == 1
        visualizeLScrosssection(wrench2plotNormalized,twist2plotNormalized,fittedRes4th.fittedCoeff,4)
    end

    %%
    if config.ifFitLSEllip == 1

        visualizeLScrosssection(wrench2plotNormalized,twist2plotNormalized,fittedResEllip.fittedCoeff,2)
    end
end
end

