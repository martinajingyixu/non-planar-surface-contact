function [normalArray] = computenormalelements(nodesCoordMatrix,elementArray)


numElement = size(elementArray,2);
normalArray = zeros(3,numElement);

    for iElement = 1:numElement

        p1 = nodesCoordMatrix(:,elementArray(1,iElement));
        p2 = nodesCoordMatrix(:,elementArray(2,iElement));
        p3 = nodesCoordMatrix(:,elementArray(3,iElement));
        normalArray(:,iElement)=plannormaltriangular(p1,p2,p3);

    end
end

function [normal] = plannormaltriangular(p1,p2,p3)
%plannormaltriangular Compute plan normal of 3 points
%order of the nodes must be clockwise!
    normal = cross(p1-p2, p1-p3);
    normal = normal./norm(normal); % normalize the normal
   
end

