function plot6dwrenches(wrench1,wrench2, color)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% subplot(3,3,1)
if size(wrench1,1) == 6
    wrench1 = wrench1';
end


if nargin == 2
    color = [0.8500,0.3250,0.0980];
elseif nargin == 1
    wrench2 = [];
    color = [0.8500,0.3250,0.0980];    
end

if size(wrench2,1) == 6
    wrench2 = wrench2';
end

for iAxis1 = 1:6
    for iAxis2 = 1:6
        for iAxis3 = 1:6
            
            if iAxis1 < iAxis2 && iAxis2 < iAxis3 
                figure
                plotwrench(wrench1, wrench2, iAxis1,iAxis2,iAxis3,color)
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

function plotwrench(wrench1, wrench2, axis1,axis2,axis3,color)
    plot3(wrench1(:,axis1),...
        wrench1(:,axis2),...
        wrench1(:,axis3),'Marker', 'o', 'MarkerFaceColor', color, 'Markersize', 3, 'LineStyle', 'none');
    hold on
    axis equal
    if numel(wrench2) > 0
        plot3(wrench2(:,axis1),...
        wrench2(:,axis2),...
        wrench2(:,axis3), 'Marker', 'o', 'MarkerFaceColor', [0 .75 .75], 'Markersize', 10, 'LineStyle', 'none');
%         wrench2(:,axis3),'ro');
%     scatter3(wrench2(:,axis1),wrench2(:,axis2),...
%         wrench2(:,axis3),'MarkerEdgeColor','none','MarkerFaceColor',[0 .75 .75],'MarkerFaceAlpha',.1)
    end
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