function [contactL,contactR,transMat] = ...
    createcontactparallelgripper(radius,surfaceInfo,wrenchSampled,ifVisualize,normalForceLeftMag,normalForceRightMag)
if nargin<4
    ifVisualize = false;
end
if nargin<6
    normalForceLeftMag = 1;
    normalForceRightMag = 1;
end

% 
transMat.R_L = getRotationMatrix( 0, 0, -pi/2);
transMat.R_R = getRotationMatrix( 0, 0, +pi/2);

transMat.p_L = [0;0;0];
transMat.p_R = [0;0;0];
%% use this transform to compute the translation vector
contactL = ...
    transformcontact(surfaceInfo,wrenchSampled,transMat.R_L,transMat.p_L); % check max coord after translation

contactL.normalForce =  contactL.normalForce .*normalForceLeftMag;
contactR = ...
    transformcontact(surfaceInfo,wrenchSampled,transMat.R_R,transMat.p_R);

contactR.normalForce =  contactR.normalForce .*normalForceRightMag;
% 
% % 
transL =  max((contactL.nodes2plot)');
transR =  min((contactR.nodes2plot)');
%% move obj center to [0;0;0]

% transform the contact
transMat.p_L = [radius-transL(1);0;0];
transMat.p_R = [-radius-transR(1);0;0];
% 
% % transform the contact
contactL = ...
    transformcontact(surfaceInfo,wrenchSampled,transMat.R_L,transMat.p_L);
contactR = ...
    transformcontact(surfaceInfo,wrenchSampled,transMat.R_R,transMat.p_R);



if ifVisualize
    figure
    plot3(surfaceInfo.elementsInfo.nodesCoordArray(1,:),...
    surfaceInfo.elementsInfo.nodesCoordArray(2,:),surfaceInfo.elementsInfo.nodesCoordArray(3,:),'g.');
    hold on
    vectarrow(surfaceInfo.frictionCenter,surfaceInfo.frictionCenter+surfaceInfo.normalForce,'g');
    hold on
    plot3(contactL.nodes2plot(1,:),...
    contactL.nodes2plot(2,:),contactL.nodes2plot(3,:),'r.');
    hold on
    vectarrow(contactL.frictionCenter,contactL.frictionCenter+contactL.normalForce,'r');
    hold on
    hold on
    plot3(contactL.frictionCenter(1),contactL.frictionCenter(2),contactL.frictionCenter(3),'r.','MarkerSize',20);
    hold on
    plot3(contactR.frictionCenter(1),contactR.frictionCenter(2),contactR.frictionCenter(3),'r.','MarkerSize',20);
    hold on
    vectarrow(contactR.frictionCenter,contactR.frictionCenter+contactR.normalForce,'b');
    hold on
    plot3(contactR.nodes2plot(1,:),...
        contactR.nodes2plot(2,:),contactR.nodes2plot(3,:),'b.');
    hold on
    axis equal;
    hold on
    xlabel('x') % x-axis label
    ylabel('y') % y-axis label
    zlabel('z') % z-axis label
end
%%
end

