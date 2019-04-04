% function [elementsInfo,meshFineLevel] = getdiscretizationresults(computedParams,sizeElem,ifConnectFirstLastNode,ifVisualize)
function [elementsInfo,surfaceInfo] = ...
    getdiscretizationresults(surfaceInfo,pathResults,ifVisualize,sizeElem)
    syms x y z real

    if nargin<3
        ifVisualize = false;
    end
    if nargin<4
%         ifVisualize = true;
%         sizeElem = mean([surfaceInfo.r1,surfaceInfo.r2,surfaceInfo.r3])/3;
        sizeElem = 0.3;
    end
    
    [nodesCoord,nodesCoordMatrix,idxNodeRows] = ...
        samplesurface(surfaceInfo.rToSample,surfaceInfo.rangeU,surfaceInfo.rangeV,sizeElem,0);

    % element info
    elementIdxArray = meshsurface(nodesCoordMatrix,idxNodeRows);

    numElements = length(elementIdxArray);
    areaArray = double(computareaelements(nodesCoordMatrix,elementIdxArray));
    normalArray = computenormalelements(nodesCoordMatrix,elementIdxArray);
    centerArray = computecenterelements(nodesCoord,elementIdxArray);

%     meshFineLevel = computedParams.surfaceInfo.surfaceArea - sum(areaArray);
    % compute pressure for each element. 
    % uniform pressure distribution: Its uniform pressure distribution,
    % and the sum of the force is 1
    % Hetzian pressure distribution
    com = mean(centerArray,2); % com is the pressure center

 
    if strcmp(surfaceInfo.pressureType,'hetzCenter') % check if pressure is uniform or hetztian
%         radiusPressure = max(sqrt(sum(bsxfun(@minus, centerArray,...
%         repmat(com,[1,length(centerArray)])).^2,1)))*1.01;
        radiusPressure = max([surfaceInfo.r1,surfaceInfo.r2,surfaceInfo.r3])*1.8; 
        % the radius, same as in integrated friction

        pressure2Sample(x,y,z) = hetzianpressure([x;y;z],radiusPressure,[com(1:2);z]); 
        pressureArray = double(pressure2Sample(centerArray(1,:),centerArray(2,:),centerArray(3,:)));
        
    elseif strcmp(surfaceInfo.pressureType,'hetzCenterMaxDist') % check if pressure is uniform or hetztian
        radiusPressure = max(sqrt(sum(bsxfun(@minus, centerArray,...
        repmat(com,[1,length(centerArray)])).^2,1)))*1.01;
%         radiusPressure = max([surfaceInfo.r1,surfaceInfo.r2,surfaceInfo.r3])*1.8; 
        % the radius, same as in integrated friction

        pressure2Sample(x,y,z) = hetzianpressure([x;y;z],radiusPressure,[com(1:2);z]); 
        pressureArray = double(pressure2Sample(centerArray(1,:),centerArray(2,:),centerArray(3,:)));
    elseif strcmp(surfaceInfo.pressureType,'uniform') 
        pressureArray = ones([1,numElements]);     
    else
        error('unknown pressure type! only accept uniform and hetztian distribution')
    end

    % scale the pressure such that the sum of the normal force is 1N.
%     pressureOneElement = 1/sum(areaArray) ;
%     pressureArray = double(repmat(pressureOneElement,[1,numElements]));
    % p0 is the scaling factor, such that is normal force is 1. 
%     normalForceMagnitude = norm(sum(normalArray .* areaArray .* pressureArray,2));
    p0 = 1/sum(areaArray .* pressureArray); 
    pressureArray = p0 .* pressureArray;
    % compute the pressure weighted friction center
    pwc = double(computefrictionweightedpressurecenter(centerArray,pressureArray,areaArray));    
    
    elementsInfo.elementIdxArray = elementIdxArray;
    elementsInfo.nodesCoordArray = nodesCoordMatrix;
    elementsInfo.pressureArray = pressureArray;
    elementsInfo.areaArray = areaArray;
    elementsInfo.normalArray = normalArray;
    elementsInfo.com = com;
    elementsInfo.centerArray = centerArray;
    elementsInfo.lengthArray = max(centerArray')' - min(centerArray')';
    elementsInfo.coord.x = nodesCoordMatrix(1,:);
    elementsInfo.coord.y = nodesCoordMatrix(2,:);
    elementsInfo.coord.z = nodesCoordMatrix(3,:);
    elementsInfo.pwc = pwc;
    
    surfaceInfo.frictionCenter = pwc; 
    surfaceInfo.pwc = pwc;
    surfaceInfo.normalForce = sum(normalArray .* areaArray .*pressureArray,2);

if ifVisualize

    showtriangularelements(elementsInfo.coord,elementsInfo.elementIdxArray,...
        elementsInfo.centerArray,elementsInfo.normalArray);
    hold on
    axis equal
    plot3(elementsInfo.pwc(1),elementsInfo.pwc(2),elementsInfo.pwc(3),'bo');
    hold on
    plot3(elementsInfo.com(1),elementsInfo.com(2),elementsInfo.com(3),'ro');
    hold on
    C = repmat([100],[1,numel(centerArray(1,:))]);
        scatter3(centerArray(1,:),centerArray(2,:),...
        centerArray(3,:),C,pressureArray,'filled')
    hold on
    axis equal
    colorbar
    hold on
    xlabel('x')
    ylabel('y')
    zlabel('z')

end


if nargin >1
    save([pathResults  'varsParametricForm.mat'],'surfaceInfo','elementsInfo');
end  
end

