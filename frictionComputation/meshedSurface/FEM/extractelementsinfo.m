function [elements,elementsParts,contact] = extractelementsinfo(contactInfo,nodal)
 
% elements contains the information of all elements
% elementsParts contains the information of different simulation components
% e.g. object, steel, finger
    % get element info
    % each element consists of 4 nodes. Node index here

    elements.nodeComponents = contactInfo.element_info(:,2:5);

    elements.pressure = contactInfo.element_pressure(:,2);
     
    for iNodeOfOneElement = 1:4
        
        elements.nodesCoord.x(:,iNodeOfOneElement) = nodal.deformedCoord...
            .x(elements.nodeComponents(:,iNodeOfOneElement));
        elements.nodesCoord.y(:,iNodeOfOneElement) = nodal.deformedCoord...
            .y(elements.nodeComponents(:,iNodeOfOneElement));
        elements.nodesCoord.z(:,iNodeOfOneElement) = nodal.deformedCoord...
            .z(elements.nodeComponents(:,iNodeOfOneElement));
        
        elements.nodesCoordNoDeformation.x(:,iNodeOfOneElement) = nodal.coord...
            .x(elements.nodeComponents(:,iNodeOfOneElement));
        elements.nodesCoordNoDeformation.y(:,iNodeOfOneElement) = nodal.coord...
            .y(elements.nodeComponents(:,iNodeOfOneElement));
        elements.nodesCoordNoDeformation.z(:,iNodeOfOneElement) = nodal.coord...
            .z(elements.nodeComponents(:,iNodeOfOneElement));

    end
    
%     elements.contactPairs.idx1 = contactInfo.element_real(contactInfo.element_real(:,2) == 1,1);
%     elements.contactPairs.idx2 = contactInfo.element_real(contactInfo.element_real(:,2) == 2,1);
%     elements.contactPairs.idx3 = contactInfo.element_real(contactInfo.element_real(:,2) == 3,1);
%     elements.contactPairs.idx4 = contactInfo.element_real(contactInfo.element_real(:,2) == 4,1);

    
    % remove the overlap of the object and the finger 
    idxFinger = find(contactInfo.element_mat(:,2) == 3);
    idxObj = find(contactInfo.element_mat(:,2) == 1);

    elementsParts.object.nodeComponents =  elements.nodeComponents(idxObj,:);
    elementsParts.finger.nodeComponents =  elements.nodeComponents(idxFinger,:);
    elementsParts.finger.elementsAreaArray = contactInfo.element_area_rough(idxFinger);


    %% check overlap of nodeComponents from object and finger. 
    % And remove it from the elementParts.object
    elementsParts.object.nodesCoord = ...
    extractrowsfieldsxyz(elements.nodesCoord,idxObj);

%     idxMatObject = contactInfo.element_mat(:,2) == 1;
    elementsParts.object.nodesCoordNoDeformation = ...
        extractrowsfieldsxyz(elements.nodesCoordNoDeformation,idxObj);
  
    elementsParts.finger.nodesCoord = ...
        extractrowsfieldsxyz(elements.nodesCoord,idxFinger);

    elementsParts.finger.nodesCoordNoDeformation = ...
        extractrowsfieldsxyz(elements.nodesCoordNoDeformation,idxFinger);

    elementsParts.object.pressureArray = elements.pressure(idxObj);
    elementsParts.finger.pressureArray = elements.pressure(idxFinger);
    
    elementsParts.finger.nodal.idxArray = unique(elementsParts.finger.nodeComponents);    
    
    elementsParts.finger.nodal.coord = ...
        extractrowsfieldsxyz(nodal.deformedCoord,elementsParts.finger.nodal.idxArray);
    
    
%     elementsParts.finger.nodal.coord = ...
%         extractrowsfieldsxyz(nodal.deformedCoord,idxFinger);
    
%     elementsParts.finger.nodal.idxArray =1:length(elementsParts.finger.nodal.coord.x);    

%% Finger pressure counts
    % extract contact area of the finger, where pressure larger than zero
    % use absolute pressure because sometimes its negative
    elementsParts.finger.pressureArray = abs(elementsParts.finger.pressureArray);
    isNonZeroPressure = elementsParts.finger.pressureArray>0;
 
    
    contact.nodeComponents = ...
        elementsParts.finger.nodeComponents(isNonZeroPressure,:);
    
    contact.nodesCoord = ...
        extractrowsfieldsxyz(elementsParts.finger.nodesCoord,isNonZeroPressure);
        
    contact.pressureArray = elementsParts.finger.pressureArray(isNonZeroPressure);
    contact.elementsAreaArray = elementsParts.finger.elementsAreaArray(isNonZeroPressure);
    
    contact.nodal.idxArray = unique(contact.nodeComponents);
    contact.nodal.coord = ...
        extractrowsfieldsxyz(nodal.deformedCoord,contact.nodal.idxArray);
    

end

function [conditions] = elementidxfoamremoveoverlap(contactInfo,elements)
% get elements which belong to both finger and no overlap
% the overlap is defined by the FEM contact info
% contact 3 and 4 are between the finger and steel
% contact 1 and 2 are between object and the finger. But half of the object
% is defind within the contact so can not use this for overlap removal for
% object and the finger.
    isfoam = contactInfo.element_mat(:,2) == 3;
    [isInContact1, isInContact2]= deal(true(length(elements.nodeComponents),1));
    isInContact1(elements.contactPairs.idx3) = false;
    isInContact2(elements.contactPairs.idx4) = false;
    conditions= isfoam == 1 & isInContact1 == 1 & isInContact2 == 1;
end

function [overlapElementIdxObj] = overlapnodesidx(elementsParts)

% find the overlap 

    [~,iax,~] = intersect(elementsParts.object.nodesCoord.x,elementsParts.finger.nodesCoord.x,'rows');
    [~,iay,~] = intersect(elementsParts.object.nodesCoord.y,elementsParts.finger.nodesCoord.y,'rows');
    [~,iaz,~] = intersect(elementsParts.object.nodesCoord.z,elementsParts.finger.nodesCoord.z,'rows');
%     [~,idx,~] = intersect(intersect(ia,ia2),ia3);

%     [overlap,iasteelx,ibx] = intersect(int64(elementsParts.steel.nodesCoord.x*100),int64(elementsParts.finger.nodesCoord.x*100),'rows');
%     [overlap,iasteelx,ibx] = intersect(elementsParts.steel.nodesCoord.x,elementsParts.finger.nodesCoord.x,'rows');
% 
%     [~,~,iby] = intersect(elementsParts.steel.nodesCoord.y,elementsParts.finger.nodesCoord.y,'rows');
%     [~,~,ibz] = intersect(elementsParts.steel.nodesCoord.z,elementsParts.finger.nodesCoord.z,'rows');

    overlapElementIdxObj = unique([iax;iay;iaz]);  
%     overlapElementIdxFoam = unique([ibx;iby;ibz]);    

%     [~,iatest,~] = intersect(testx,floor(elementsParts.finger.nodesCoord.x*10),'rows');

%     overlapElementIdx = ia(idx);
end