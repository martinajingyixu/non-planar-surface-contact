function [left,right] = separatecontacts(contact,isSqueezingX)
% separate the contact for the left and right fingers
% Unit: mm and mPa
if isSqueezingX
    elementsIdxLeft = contact.nodesCoord.x(:,1) < 0; 
    % one node of the four is enough for distinguish left and right contact
    elementsIdxRight = contact.nodesCoord.x(:,1) >= 0;
    nodesIdxLeft = contact.nodal.coord.x < 0;
    nodesIdxRight = contact.nodal.coord.x >= 0;
   
else
    elementsIdxLeft = contact.nodesCoord.y(:,1) < 0;
    elementsIdxRight = contact.nodesCoord.y(:,1) >= 0;
    nodesIdxLeft = contact.nodal.coord.y < 0;
    nodesIdxRight = contact.nodal.coord.y >= 0;
        
end
%% left contact
    left.nodal.coord = extractrowsfieldsxyz(contact.nodal.coord,nodesIdxLeft);   
    left.nodal.idxArray = contact.nodal.idxArray(nodesIdxLeft);    

    left.elements.nodesCoord = extractrowsfieldsxyz(contact.nodesCoord,elementsIdxLeft);
    left.elements.nodeComponents = contact.nodeComponents(elementsIdxLeft,:);
    left.elements.pressureArray = contact.pressureArray(elementsIdxLeft);
    left.elements.areaArray = contact.elementsAreaArray(elementsIdxLeft);

    coord = left.elements.nodesCoord;
   left.elements.centerArray = structfun(@(coord) mean(coord,2),coord,'UniformOutput', false);
%     
    left.elements.nodesCoord = extractrowsfieldsxyz(contact.nodesCoord,elementsIdxLeft);

    
    left.contactCenter = [mean(left.elements.centerArray.x);...
        mean(left.elements.centerArray.y);mean(left.elements.centerArray.z)];
     left.com = [mean(left.elements.centerArray.x);...
        mean(left.elements.centerArray.y);mean(left.elements.centerArray.z)]; 
    
    
    % compute sum of the force and normalize contact areaArray
    left.sumForce = computesumforce(left.elements.areaArray, left.elements.pressureArray);
    left.elements.pressureArrayNormalized = left.elements.pressureArray./left.sumForce;
    
%% right contact
  
    right.nodal.coord = extractrowsfieldsxyz(contact.nodal.coord,nodesIdxRight);
    right.nodal.idxArray = contact.nodal.idxArray(nodesIdxRight);    

    right.elements.nodesCoord = extractrowsfieldsxyz(contact.nodesCoord,elementsIdxRight);

    right.elements.nodeComponents = contact.nodeComponents(elementsIdxRight,:);
    right.elements.pressureArray = contact.pressureArray(elementsIdxRight);
    right.elements.areaArray = contact.elementsAreaArray(elementsIdxRight);
    
    coord = right.elements.nodesCoord;
    right.elements.centerArray = ...
        structfun(@(coord) mean(coord,2),coord,'UniformOutput', false);

    right.contactCenter = [mean(right.elements.centerArray.x);...
        mean(right.elements.centerArray.y);mean(right.elements.centerArray.z)];
    
    right.sumForce = computesumforce(right.elements.areaArray, right.elements.pressureArray);
    right.elements.pressureArrayNormalized = right.elements.pressureArray./right.sumForce;


end

function [sumForce] = computesumforce(elementsareaArray,pressureArray)
    sumForce = sum(mm2tom2(elementsareaArray).* mpa2pa(pressureArray));
end

function [pressureArrayInPa] = mpa2pa(pressureArrayInMPa)
pressureArrayInPa = pressureArrayInMPa.*1e6;
end

function [areaArrayInm2] = mm2tom2(areaArrayInmm2)
areaArrayInm2 = areaArrayInmm2.*1e-6;
end


