function [wrench] = computewrenchdiscrete(elementArray,twistSamples,ifVisualize)
% COR data preparation
if nargin <3
    ifVisualize = false;    
end

numElements = size(elementArray.elementIdxArray,2);
elementInfoArray = {1,numElements} ;
for iElement = 1:numElements
    elementInfo.idx = elementArray.elementIdxArray(:,iElement);
    elementInfo.area = elementArray.areaArray(:,iElement);
    elementInfo.center = elementArray.centerArray(:,iElement);
    elementInfo.normal = elementArray.normalArray(:,iElement);
    elementInfo.pressure = elementArray.pressureArray(:,iElement);
    elementInfoArray{iElement} = elementInfo;
end

mu = 1; % not dealing with mu now

wrench = zeros(6,length(twistSamples));
% wrench = [];
% pitchValues = [0];

% for idxTwist = 1:length(twistSamples)
parfor idxTwist = 1:length(twistSamples)
% for idxTwist = 1

    idxTwist;
    [forceArray,torqueArray] = frictiondiscretizedsurface...
        (twistSamples(:,idxTwist),elementInfoArray,elementArray.pwc,numElements,mu);
    wrench(:,idxTwist) = [sum(forceArray'),sum(torqueArray')]';
%     wrench(:,(idxCOR-1)*numel(pitchValues)+1:idxCOR*numel(pitchValues)) = wrenchOnePitch;

end  
if ifVisualize
    wrenchDiscreteSurfacePlot = wrench';

    plot6dwrenches(wrenchDiscreteSurfacePlot./max(wrenchDiscreteSurfacePlot));
    
%         wrenchDiscreteSurfacePlot = wrench';
%     plot6dwrenches(wrenchDiscreteSurfacePlot./max(wrenchDiscreteSurfacePlot),...
%         wrenchDiscreteSurfacePlot(end-124:end,:)./max(wrenchDiscreteSurfacePlot));
end
% save([pathWrench '/wrenchTwist.mat'],'wrench','twistSamples');
end
%%

