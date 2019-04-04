function [areaArray] = computareaelements(nodesCoordMatrix,elementArray)
%Given nodes coordinate and idx of element, compute the area of each
%element

numElement = size(elementArray,2);
areaArray = zeros(1,numElement);

    for iElement = 1:numElement

        p1 = nodesCoordMatrix(:,elementArray(1,iElement));
        p2 = nodesCoordMatrix(:,elementArray(2,iElement));
        p3 = nodesCoordMatrix(:,elementArray(3,iElement));
        areaArray(1,iElement)=getTriangleArea(p1,p2,p3);

    end

end

function [area]=getTriangleArea(p1,p2,p3)
% This function gives the area of a triangle given 3 vertics coordinates
%http://math.stackexchange.com/questions/128991/how-to-calculate-area-of-3d-triangle

    AB = p1-p2;
    AC = p1-p3;

    x1 = AB(1);
    x2 = AB(2);
    x3 = AB(3);

    y1 = AC(1);
    y2 = AC(2);
    y3 = AC(3);

    area = 0.5*sqrt((x2*y3-x3*y2)^2 + (x3*y1-x1*y3)^2 + (x1*y2-x2*y1)^2);

end

