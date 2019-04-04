function [graspWrenchSpaceMink,graspWrenchSpaceUnion,wrenchL,wrenchR] =...
    computegws(contactL,contactR,mu,objOrigin)
% compute the grasp wrench space of the left right finger contact.
% Now is the torque with respect to the origin of the object.

if nargin <4
    objOrigin = [0;0;0];
end

% Compute the wrench with respect to the object origin
%MB = MA + rAB x F

wrenchL = computewrenchobj(contactL,objOrigin,mu);
wrenchR = computewrenchobj(contactR,objOrigin,mu);
wrenchL = [wrenchL,[0;0;0;0;0;0]];
wrenchR = [wrenchR,[0;0;0;0;0;0]];

graspWrenchSpaceMink = minkowskisum(wrenchL,wrenchR);
% graspWrenchSpaceMink = [];
% figure
% k = convhulln(graspWrenchSpaceMink(4:6,:)');
% trisurf(k,graspWrenchSpaceMink(4,:)',graspWrenchSpaceMink(5,:)',graspWrenchSpaceMink(6,:)')
graspWrenchSpaceUnion =  union(wrenchL',wrenchR','Rows')';

end

function wrenchOriginObj = computewrenchobj(contact,objOrigin,mu)
    contact.frictionWrench = contact.frictionWrench .* norm(contact.normalForce).*mu;
    torqueArm = contact.frictionCenter - objOrigin;
    force = contact.frictionWrench(1:3,:) + contact.normalForce;
    torque = cross(repmat(torqueArm,[1,length(force)]),force) + contact.frictionWrench(4:6,:);
    wrenchOriginObj = [force;torque];
end

