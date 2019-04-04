function [contact] = computecontactareainfo(contact)
%compute the contact center, contact normal, pressure-weighted fricion
%center of one contact area
contact.pwc = computefrictionweightedpressurecenter(contact);
[contact.centerNode, contact.centerNodeIdx] = centernode(contact);
[contact.normal,contact.normalAverageElementsNormal] ...
    = normalcontactarea(contact);
% contact.projectedPWC = projectedpwc(contact);

%% in local coordinate system
contact.localCoordSystem = localcoordinatesystem(contact);
% [contact.transformationMatrix, contact.localCoord.normal] ...
%     = transformationmatrix(contact);
[contact.transformationMatrix] = transformationmatrix(contact);
end

