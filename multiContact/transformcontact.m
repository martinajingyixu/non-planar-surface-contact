function [contact] = transformcontact(surface,wrench,RRot,pTrans)
% We have the contact area and friction in at the origin.
% Now we rotate and translate the contact into another pose, so we have a grasp of
% parallel gripper

% here we rotate the friction center, the surface, and the friction wrench
% to the new pose

%% rotate and translate the contact

contact.frictionCenter = ...
    RRot * surface.frictionCenter + pTrans;
if isfield(surface,'rSampledCatesianCoord')
    contact.nodes2plot = ...
        RRot * surface.rSampledCatesianCoord + pTrans;
elseif isfield(surface,'elementsInfo')
   contact.nodes2plot = ...
    RRot *surface.elementsInfo.nodesCoordArray + pTrans; 
end

%% rotate and translate the contact
% the friction center is moved as well, so the friction torque is still
% with respect to the same point
contact.normalForce = RRot * surface.normalForce;

pwave = [0;0;0]; % pwave is the translation of the origin or different coordinate system

RWrench = [RRot,zeros(size(RRot));...
    RRot* crossproductmatrix(pwave),RRot];
contact.frictionWrench = RWrench * wrench;
       
end

