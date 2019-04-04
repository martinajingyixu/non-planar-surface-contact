function plot3dconvechullof6dwrenches(wrench1,wrench2, color)

if nargin == 2
    color = 'b.';
elseif nargin == 1
    wrench2 = [];
    color = 'b.';

end
for iAxis1 = 1:6
    for iAxis2 = 1:6
        for iAxis3 = 1:6
            
            if iAxis1 < iAxis2 && iAxis2 < iAxis3 
                figure
                plotwrenchandhull(wrench1, wrench2, iAxis1,iAxis2,iAxis3,color)
                title([int2str(iAxis1) int2str(iAxis2) int2str(iAxis3)]);
%                 hold on
%                 savefig(['figures/unitCylinder/CORZaxisConvexHullComp' int2str(iAxis1) int2str(iAxis2) int2str(iAxis3) '.fig'])
            else
                continue
            end
        end
    end
end
end

function plotwrenchandhull(wrench1, wrench2, axis1,axis2,axis3,color)
    K1 = convhulln([wrench1(axis1,:);wrench1(axis2,:);wrench1(axis3,:)]');
    K1unique = unique(K1);

%     plot3(wrench1(axis1,K1unique),...
%         wrench1(axis2,K1unique),...
%         wrench1(axis3,K1unique),color);
    
    trisurf(K1,wrench1(axis1,:)',wrench1(axis2,:)',wrench1(axis3,:)','FaceColor','cyan');

    hold on
    if numel(wrench2) > 0
%         K2 = convhulln([wrench2(axis1,:);wrench2(axis2,:);wrench2(axis3,:)]');
%         K2unique = unique(K2);
        plot3(wrench2(axis1,:),...
        wrench2(axis2,:),...
        wrench2(axis3,:),'rx');
    end
    axis equal
end