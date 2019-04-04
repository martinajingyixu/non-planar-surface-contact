function [elementsInfo] = convertelementinfo(contact)
elementIdxArray = [contact.elements.nodeComponents(:,[1,4,3]);...
    contact.elements.nodeComponents(:,[1,3,2])]';

nodesCoordArray = zeros(max(contact.nodal.idxArray),3);
nodesCoordArray(contact.nodal.idxArray,:) = ...
    [contact.nodal.coord.x,contact.nodal.coord.y,contact.nodal.coord.z];

pressureArray = [contact.elements.pressureArray;contact.elements.pressureArray]';

coord.x = nodesCoordArray(:,1)';
coord.y = nodesCoordArray(:,2)';
coord.z = nodesCoordArray(:,3)';

areaArray = computareaelements(nodesCoordArray',elementIdxArray);
normalArray = computenormalelements(nodesCoordArray',elementIdxArray);
centerArray = computecenterelements(coord,elementIdxArray);
pwc = computefrictionweightedpressurecenter(centerArray,pressureArray,areaArray);
lengthArray = max(nodesCoordArray) - min(nodesCoordArray);

elementsInfo.elementIdxArray = elementIdxArray;
elementsInfo.nodesCoordArray = nodesCoordArray;
elementsInfo.pressureArray = pressureArray;
elementsInfo.coord = coord;
elementsInfo.areaArray = areaArray;
elementsInfo.normalArray = normalArray;
elementsInfo.pwc = pwc;
elementsInfo.centerArray = centerArray;
elementsInfo.lengthArray = lengthArray;
elementsInfo.com = contact.contactCenter;


end

