function  showfrictiondirection...
    (nodesCoord,elementInfoArray,forceArray,projectedVelArray,velArray,CORAxis)
% figure
    hold on

for iElement = 1:length(elementInfoArray)
        hold on

    idxArray = elementInfoArray{iElement}.idx;
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
    hold on;

    vectarrow(elementInfoArray{iElement}.center,forceArray(:,iElement).*5+elementInfoArray{iElement}.center,'b');
    hold on
%     vectarrow(elementInfoArray{iElement}.center,...
%         projectedVelArray(:,iElement).*0.3+elementInfoArray{iElement}.center,'g');
%     hold on
    vectarrow(elementInfoArray{iElement}.center,...
        velArray(:,iElement).*0.3+elementInfoArray{iElement}.center,'c');
    hold on

    plot3(elementNodesCoords.x,elementNodesCoords.y,elementNodesCoords.z,'black');
    hold on
%     elementNormal = elementInfoArray{iElement}.normal.*0.1;
%     vectarrow(elementInfoArray{iElement}.center,elementInfoArray{iElement}.center+elementNormal,'g');

%     axis equal;
    hold on;
    
 
end
    CORP0 = CORAxis.location - 0.5 .*CORAxis.direction;
    CORP1 = CORAxis.location + 0.5 .*CORAxis.direction;
    hold on

%     vectarrow(CORP0,CORP1,'g');
    hold on

end

