function [] = showtriangularelements(nodesCoord,elementArray,elementCenterArray,normalArray,ifNewfigure,color)

idxNodesElementArray = elementArray;
if ~exist('color', 'var') || isempty(color)
    color = 'b';
end

if ~exist('ifNewfigure', 'var') || isempty(ifNewfigure)
    ifNewfigure = true;
end

if ~exist('normalArray', 'var') || isempty(normalArray)
    showNormal = false;
else
    showNormal = true;
end

if ifNewfigure
    figure;
end
% figure;
% color = 'b';
for iElement = 1:size(idxNodesElementArray,2)
    
    idxArray = idxNodesElementArray(:,iElement);
%     iElement
%     idxArray
    elementNodesCoords.x = [...
        nodesCoord.x(idxArray(1));...
        nodesCoord.x(idxArray(2));...
        nodesCoord.x(idxArray(3));...
        nodesCoord.x(idxArray(1))];
    
    elementNodesCoords.y = [...
        nodesCoord.y(idxArray(1));...
        nodesCoord.y(idxArray(2));...
        nodesCoord.y(idxArray(3));...
        nodesCoord.y(idxArray(1))];
    
    
    elementNodesCoords.z = [...
        nodesCoord.z(idxArray(1));...
        nodesCoord.z(idxArray(2));...
        nodesCoord.z(idxArray(3));...
        nodesCoord.z(idxArray(1))];
    
    if showNormal
        elementNormal = normalArray(:,iElement).*0.1;
        vectarrow(elementCenterArray(:,iElement),elementCenterArray(:,iElement)+elementNormal,'g');
    end
    hold on

    plot3(elementNodesCoords.x,elementNodesCoords.y,elementNodesCoords.z,color,'LineWidth',1.5);
    hold on

%     plot3(mean(elementNodesCoords.x),mean(elementNodesCoords.y),mean(elementNodesCoords.z),[color 'o'])
%     axis equal;
    hold on;
end

