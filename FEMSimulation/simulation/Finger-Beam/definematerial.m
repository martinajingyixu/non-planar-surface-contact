function [] = definematerial(filename)
%definematerial Define material for ANSYS for Object material, sensor and
%steel to prevent undisired deformation

    global CONSTANTS

    fileID = fopen(filename,'a');
    fprintf(fileID,'\n');
      
    fprintf(fileID,'MP,EX,1,%f\n',CONSTANTS.MATERIAL.YOUNGS_MODULUS); % Mesh
    fprintf(fileID,'MP,NUXY,1,%f\n',CONSTANTS.MATERIAL.NUXY);
    fprintf(fileID,'MP,MU,1,%f\n',CONSTANTS.MATERIAL.MU);
    fprintf(fileID,'MP,DENS,1,%f\n',CONSTANTS.MATERIAL.DENSITY);
    fprintf(fileID,'\n'); % MP, EX 1 is the material for object

%     fprintf(fileID,'MP,EX,2,%f\n',CONSTANTS.SIMULATION.STEEL_STIFFNESS); % Mesh
%     fprintf(fileID,'MP,NUXY,2,%f\n',CONSTANTS.SIMULATION.STEEL_NUXY);
%     fprintf(fileID,'MP,MU,2,%f\n',CONSTANTS.SIMULATION.STEEL_MU);
%     fprintf(fileID,'MP,DENS,2,%f\n',CONSTANTS.SIMULATION.STEEL_DENSITY);
    fprintf(fileID,'\n'); % MP, EX 1 is the material for object

    % MP, EX 3 is the material for sensor
    fprintf(fileID,'MP,EX,3,%f\n',CONSTANTS.SIMULATION.SENSOR_STIFFNESS);
    fprintf(fileID,'MP,NUXY,3,%f\n',CONSTANTS.SIMULATION.SENSOR_NUXY);
    fprintf(fileID,'MP,MU,3,%f\n',CONSTANTS.SIMULATION.SENSOR_MU);
    fprintf(fileID,'MP,DENS,3,%f\n',CONSTANTS.SIMULATION.SENSOR_DENSITY);
    

    fprintf(fileID,'\n');
end

