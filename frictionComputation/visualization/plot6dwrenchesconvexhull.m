function plot6dwrenchesconvexhull(wrench1,idxConvexHull, color)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% subplot(3,3,1)
if nargin == 2
    color = 'b.';
elseif nargin == 1
    idxConvexHull = [];
    color = 'b.';

end
for iAxis1 = 1:6
    for iAxis2 = 1:6
        for iAxis3 = 1:6
            
            if iAxis1 < iAxis2 && iAxis2 < iAxis3 
                figure
                plotwrench(wrench1, iAxis1,iAxis2,iAxis3,color)
                hold on
                if ~isempty(idxConvexHull)
                    trimesh(idxConvexHull(:,[iAxis1,iAxis2,3]),wrench1(:,iAxis1),...
                    wrench1(:,iAxis2),wrench1(:,iAxis3)) 
                end
                title([int2str(iAxis1) int2str(iAxis2) int2str(iAxis3)]);
                xlabel(getlabel(iAxis1))
                ylabel(getlabel(iAxis2))
                zlabel(getlabel(iAxis3))
%                 hold on
%                 savefig(['figures/unitCylinder/CORZAxisCompare' int2str(iAxis1) int2str(iAxis2) int2str(iAxis3) '.fig'])
            else
                continue
            end
        end
    end
end
end

function plotwrench(wrench1, axis1,axis2,axis3,color)
    plot3(wrench1(:,axis1),...
        wrench1(:,axis2),...
        wrench1(:,axis3),color);
    hold on
    axis equal
end

function label = getlabel(iAxis)
switch iAxis
    case 1
        label = 'fx';
    case 2
        label = 'fy';   
    case 3
        label = 'fz';           
    case 4
        label = 'taux';           
    case 5
        label = 'tauy';           
    case 6
        label = 'tauz'; 
end
        
end