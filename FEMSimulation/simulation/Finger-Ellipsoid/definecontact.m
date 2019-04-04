function [] = definecontact(filename,MODEL_PARAMETERS,contactDist1,contactDist2)
%definecontact define contact of ansys models


    global CONSTANTS
    
    fileID = fopen(filename,'a');
    fprintf(fileID,'\n');
    
     fprintf(fileID,'ET,4,targe170\n'); % target 170 is the target surface to be contacted
    fprintf(fileID,'ET,5,conta174\n');

%     fprintf(fileID,'keyopt,5,12,3\n');
% %     fprintf(fileID,'keyopt,5,5,1\n');
%     fprintf(fileID,'keyopt,5,5,3\n');
% 
%     fprintf(fileID,'keyopt,5,7,1\n');
%     fprintf(fileID,'keyopt,5,9,1\n');
%     fprintf(fileID,'keyopt,5,11,1\n');

    fprintf(fileID,'keyopt,5,4,2\n');
    fprintf(fileID,'keyopt,5,5,3\n');

    fprintf(fileID,'keyopt,5,6,1\n');
    fprintf(fileID,'keyopt,5,7,1\n');
    fprintf(fileID,'keyopt,5,8,2\n');

    fprintf(fileID,'keyopt,5,9,1\n');
    fprintf(fileID,'keyopt,5,10,2\n');

    fprintf(fileID,'keyopt,5,11,1\n');
    fprintf(fileID,'keyopt,5,12,3\n');



    fprintf(fileID,'r,1,,1,1.0,0.1\n');
    fprintf(fileID,'r,2,,1,1.0,0.1\n');
    fprintf(fileID,'r,3,,1,1.0,0.1\n');
    fprintf(fileID,'r,4,,1,1.0,0.1\n');
    

    fprintf(fileID,'vsel,s,volu,,1,8\n');
    % comment the for a automatic mesh because of the complex geometry of
    % the finger
%     fprintf(fileID,'ESIZE,%f\n',CONSTANTS.SIMULATION.ELEMENT_SZ_FINGER); % Mesh

    fprintf(fileID,'TYPE,10\n');
    fprintf(fileID,'MAT,3\n');
    fprintf(fileID,'VMESH,1,8\n');

    fprintf(fileID,'csys,0\n');
    fprintf(fileID,[MODEL_PARAMETERS.SELECT_POSITIVE_HALF num2str(contactDist1-0.01) '\n']); %Define material and object part in contact, dist from center to finger pos -0.01
    fprintf(fileID,'nsla,s,1\n');
    fprintf(fileID,'TYPE,4\n');
    fprintf(fileID,'MAT,1\n');
    fprintf(fileID,'real,1\n');
    fprintf(fileID,'ESURF\n');
    % each econtact is 1/8 of ellipsoid, 4 times contact each finger
    fprintf(fileID,'ASEL,s,area,,1001,1010,3\n');
    fprintf(fileID,'nsla,s,1\n');
    fprintf(fileID,'TYPE,5\n');
    fprintf(fileID,'MAT,3\n');
    fprintf(fileID,'ESURF\n');

    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!contact pair 2
    
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!target
    fprintf(fileID,[MODEL_PARAMETERS.SELECT_NEGATIVE_HALF num2str(abs(contactDist2)-0.01) '\n']); % make num2str(cp.r1) the radius of contact
    fprintf(fileID,'nsla,s,1\n');
    fprintf(fileID,'TYPE,4\n');
    fprintf(fileID,'MAT,1\n');
    fprintf(fileID,'real,2\n');
    fprintf(fileID,'ESURF\n');

    fprintf(fileID,'ASEL,s,area,,2001,2010,3\n');
    fprintf(fileID,'nsla,s,1\n');
    fprintf(fileID,'TYPE,5\n');
    fprintf(fileID,'MAT,3\n');
    fprintf(fileID,'ESURF\n');


end

