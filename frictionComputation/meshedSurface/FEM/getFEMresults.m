function [contactResult] = getFEMresults(contactInfo,isSqueezingX,ifVisualize)
% iObject = 4;
% iSqueezeLocation = 9;
% iPressure = 6;
if nargin<3
    ifVisualize = false;
end

nodal = extractnodalinfo(contactInfo);

[elements,elementsParts,contactNotSeparated] = extractelementsinfo(contactInfo,nodal);
% left: x<0 if is squeezing along the x axis. y<0 if squezing along y axis.
[contact.left,contact.right] = separatecontacts(contactNotSeparated,isSqueezingX); % contact only has elements where pressure>0
[contactAllElements.left,contactAllElements.right] = separatecontacts(elementsParts.finger,isSqueezingX);

%% get elements info
contactResult.left = convertelementinfo(contact.left);
contactResult.right = convertelementinfo(contact.right);
contactResult.left.allElementsInfo= convertelementinfo(contactAllElements.left);
contactResult.right.allElementsInfo = convertelementinfo(contactAllElements.right);


if ifVisualize
    %%
%     figure
    %% all object and finger info
    h=figure;
    elementsLeft = contactResult.left;
    elementsRight = contactResult.right;

% %     sz = repmat([10],[1,numel(elementsParts.object.nodesCoord.x(:,1))]);
% %     scatter3(elementsParts.object.nodesCoord.x(:,1),elementsParts.object.nodesCoord.y(:,1),...
% %         elementsParts.object.nodesCoord.z(:,1),sz,sz.*0.1,'filled','MarkerFaceColor','g','MarkerFaceAlpha',.2)
% %     hold on
% 
%     sz = repmat([1],[1,numel(elementsLeft.allElementsInfo.nodesCoordArray(:,1))]);
%     scatter3(elementsLeft.allElementsInfo.nodesCoordArray(:,1),elementsLeft.allElementsInfo.nodesCoordArray(:,2),...
%         elementsLeft.allElementsInfo.nodesCoordArray(:,3),sz,sz.*0.1,'filled','MarkerFaceColor','b','MarkerFaceAlpha',.1)
%    
%     hold on
%     sz = repmat([1],[1,numel(elementsRight.allElementsInfo.nodesCoordArray(:,1))]);
%     scatter3(elementsRight.allElementsInfo.nodesCoordArray(:,1),elementsRight.allElementsInfo.nodesCoordArray(:,2),...
%         elementsRight.allElementsInfo.nodesCoordArray(:,3),sz,sz.*0.1,'filled','MarkerFaceColor','b','MarkerFaceAlpha',.1)
%     hold on
%%
    showtriangularelements(elementsLeft.coord,elementsLeft.elementIdxArray,...
        elementsLeft.centerArray,elementsLeft.normalArray,0,'b');
    hold on
    axis equal
    hold on
%     colorbar
    C = repmat([60],[1,numel(elementsLeft.centerArray(1,:))]);
    scatter3(elementsLeft.centerArray(1,:),elementsLeft.centerArray(2,:),...
        elementsLeft.centerArray(3,:),C,abs(elementsLeft.pressureArray),'filled')
%     hold on
    %%

%     plot3(elementsLeft.pwc(1),elementsLeft.pwc(2),elementsLeft.pwc(3),'bo');
%     hold on

%     plot3(elementsLeft.allElementsInfo.com(1),...
%         elementsLeft.allElementsInfo.com(2),elementsLeft.allElementsInfo.com(3),'bx');
%     hold on

%     plot3(elementsRight.pwc(1),elementsRight.pwc(2),elementsRight.pwc(3),'ro');
%     hold on
% 
%     plot3(elementsRight.allElementsInfo.com(1),...
%         elementsRight.allElementsInfo.com(2),elementsRight.allElementsInfo.com(3),'rx');
%%
    showtriangularelements(elementsRight.coord,elementsRight.elementIdxArray,...
        elementsRight.centerArray,elementsRight.normalArray,0,'b');
    hold on

    C = repmat([60],[1,numel(elementsRight.centerArray(1,:))]);
% 
    scatter3(elementsRight.centerArray(1,:),elementsRight.centerArray(2,:),...
    elementsRight.centerArray(3,:),C,abs(elementsRight.pressureArray),'filled')   
    axis equal
%     colorbar
    hold on
    exportfigure(h,'figures/FEM/')
    hold on
    %%
    
end

    
end

