function [wrench,twist] = integratewrench(pathResults,ifShow)

% pathResults = ['frictionResults/' surfaceType '/' pressureType '/'];
pathTwistSamples = [pathResults '/twistSampleVariables.mat'];
computedParams = load([pathResults '/FrictionParametricForm.mat']);

%% compute the integral of friction based on Twist samples
% get the Twist Samples. load or generate
dataTwist = load(pathTwistSamples);

twistSampleVariables = dataTwist.twistSampleVariables;
numSamples = size(twistSampleVariables,2);
%%
for idxSampleStart = 1:1000:numSamples

    idxSampleEnd = min(idxSampleStart + 999,numSamples);

    computefrictionparfor(pathResults,pathTwistSamples,...
        idxSampleStart,idxSampleEnd);
end

%% Merge the friction from multiple .mat files into a single variable/File
% merge the friction from multiple .mat files into 1

[integratedFriction.fxArray,integratedFriction.fyArray,integratedFriction.fzArray,...
    integratedFriction.tauxArray,integratedFriction.tauyArray,integratedFriction.tauzArray]...
    = deal(zeros(1,numSamples));

for idxSampleStart = 1:1000:numSamples
    idxSampleEnd = min(idxSampleStart + 999,numSamples);
    integratedFriction = mergeintegraldata...
        (pathResults,idxSampleStart,idxSampleEnd,integratedFriction);
end


% save([pathProcessingSurface '/integratedFriction' num2str(numSamples) '.mat'],'integratedFriction','dataTwist');

%% NEW
wrench =  [integratedFriction.fxArray;integratedFriction.fyArray;...
    integratedFriction.fzArray;integratedFriction.tauxArray;...
    integratedFriction.tauyArray;integratedFriction.tauzArray];
% twist = computetwist(computedParams.surfaceInfo.com,twistSampleVariables);
twist = computetwist(computedParams.surfaceInfo.frictionCenter,twistSampleVariables);

[wrenchNormalized,twistNormalized] = normalizewrenchtwist(wrench,twist);

[~,idxArray] = downsamplewrench(wrenchNormalized);
numel(idxArray)
wrench2fit = wrenchNormalized(:,idxArray);
twist2fit = twistNormalized(:,idxArray);
twist2fit = bsxfun(@rdivide, twist2fit, sqrt(sum(twist2fit.^2)));
save([pathResults, 'wrenchTwist.mat'],'twistSampleVariables','wrench','twist','idxArray',...
    'wrench2fit','twist2fit');
%%
if ifShow
wrenchPlot = wrench./max(wrench);
plot6dwrenches(wrenchPlot)
end
end


