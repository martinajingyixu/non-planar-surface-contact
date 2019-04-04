function [downSampledPoints,idxDownsampledPoints] = downsample3dpoints(pointArray,numSteps)
if nargin == 1
%     stepSize = 3;
    numSteps = 10; %77 % 60 used for parametric surface
%       stepSize = 2.2; % 60 used for FEM
%     stepSize = 2; % 52 samples


end
idxDownsampledPoints = [];

pointArrayScaled = pointArray./(0.5*(max(pointArray') - min(pointArray'))');
rangeMin = min(pointArrayScaled');
rangeMax = max(pointArrayScaled');
stepSz = ( max(pointArrayScaled') - min(pointArrayScaled'))./numSteps; % very fine downsample
% step = (rangeMax - rangeMin)./6; % default for showing the results of
% LS 
% step = (rangeMax - rangeMin)./5;

% step = (rangeMax - rangeMin)./3;


x1V = rangeMin(1):stepSz(1):rangeMax(1);
x2V = rangeMin(2):stepSz(2):rangeMax(2);
x3V = rangeMin(3):stepSz(3):rangeMax(3);

parfor ix1 = 1:numel(x1V)
    idx1 = find(pointArrayScaled(1,:)>= x1V(ix1) & pointArrayScaled(1,:)<x1V(ix1)+stepSz(1));
    if numel(idx1) == 0
        continue
    end
    for ix2 = 1:numel(x2V)
%    for ix2 = 1:numel(x2V)

        idx2 = intersect(idx1,find(pointArrayScaled(2,:)>= x2V(ix2) & pointArrayScaled(2,:)<x2V(ix2)+stepSz(2)));
        if numel(idx2) == 0
            continue
        end
        for ix3 = 1:numel(x3V)
            idx=  intersect(idx2,find(pointArrayScaled(3,:)>= x3V(ix3) & pointArrayScaled(3,:)<x3V(ix3)+stepSz(3)));
            if ~isempty(idx)
                % Find the sample which is nearst to the
                % centerid of the samples in grid
                [~,idxMinDist] = min(sqrt(sum((pointArrayScaled(:,idx) -  mean(pointArrayScaled(:,idx),2)).^2)));
                idxDownsampledPoints = [idxDownsampledPoints,idx(idxMinDist)];
            end              
        end
    end
end

idxDownsampledPoints = unique(idxDownsampledPoints);
downSampledPoints = pointArray(:,idxDownsampledPoints);

end


