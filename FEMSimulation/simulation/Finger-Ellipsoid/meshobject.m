function [] = meshobject(filename,total_areas_wall,area_idx_bottem)
    fileID = fopen(filename,'a'); 
    fprintf(fileID,'\n');
    
    global CONSTANTS


    % mesh the elements on the wall without Rim part
    fprintf(fileID,'\n');

    sel_wall_mat1 = ['ASEL,s,area,,1,' num2str(total_areas_wall-4) '\n'];
    fprintf(fileID,sel_wall_mat1);
    fprintf(fileID,'AATT,1,,1,,1\n');
    fprintf(fileID,'ESIZE,%f\n',CONSTANTS.SIMULATION.ELEMENT_SZ_OBJECT); % Mesh

    amesh_wall_mat1 = ['AMESH,1,' num2str(total_areas_wall-4) '\n'];
    fprintf(fileID,amesh_wall_mat1);
    
    % mesh the bottem
    fprintf(fileID,'\n');
    sel_bottem = ['ASEL,s,area,,' num2str(area_idx_bottem) ',' num2str(area_idx_bottem) '\n'];
    fprintf(fileID,sel_bottem);
    fprintf(fileID,'AATT,1,,1,,1\n');
    fprintf(fileID,'ESIZE,%f\n',CONSTANTS.SIMULATION.ELEMENT_SZ_OBJECT); % Mesh

    amesh_bottem = ['AMESH,'  num2str(area_idx_bottem) ',' num2str(area_idx_bottem) '\n'];
    fprintf(fileID,amesh_bottem);
    
    % mesh the rim structure
    fprintf(fileID,'\n');
    sel_wall_mat2 = ['ASEL,s,area,,' num2str(total_areas_wall-3) ',' num2str(total_areas_wall) '\n'];
    fprintf(fileID,sel_wall_mat2);
    fprintf(fileID,'AATT,1,,2,,2\n');
%     fprintf(fileID,'ESIZE,7.000000\n');
    fprintf(fileID,'ESIZE,%f\n',CONSTANTS.SIMULATION.ELEMENT_SZ_OBJECT); % Mesh

    amesh_wall_mat2 = ['AMESH,' num2str(total_areas_wall-3) ',' num2str(total_areas_wall) '\n'];
    fprintf(fileID,amesh_wall_mat2);

end