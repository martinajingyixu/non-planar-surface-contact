function [nodal] = extractnodalinfo(contactInfo)
 
    % node coordinates before deformation
    nodal.coord.x = contactInfo.nodal_coord(:,2);
    nodal.coord.y = contactInfo.nodal_coord(:,3);
    nodal.coord.z = contactInfo.nodal_coord(:,4);
    
    % nodal displacement: 
    % data in contactInfo.nodal_displacement idx, UX, UY,UZ, USUM

    nodal.dispacement.x = contactInfo.nodal_displacement(:,2);
    nodal.dispacement.y = contactInfo.nodal_displacement(:,3);
    nodal.dispacement.z = contactInfo.nodal_displacement(:,4);

    % nodal_displacement
    nodal.deformedCoord.x = nodal.coord.x + nodal.dispacement.x;
    nodal.deformedCoord.y = nodal.coord.y + nodal.dispacement.y;
    nodal.deformedCoord.z = nodal.coord.z + nodal.dispacement.z;

end