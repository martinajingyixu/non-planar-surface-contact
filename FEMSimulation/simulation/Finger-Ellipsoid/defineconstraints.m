function [] = defineconstraints(filename,MODEL_PARAMETERS,area_idx_bottem)
%defineconstraints Define the constraints of FEM simulation
    %% Constraints X axis   
    fileID = fopen(filename,'a');
    fprintf(fileID,'\n');
    fprintf(fileID,'/solu\n');
    fprintf(fileID,'allsel\n');
    fprintf(fileID,'csys,0\n');
    fprintf(fileID,'nlgeom,on\n');
    fprintf(fileID,'solc\n');

    fprintf(fileID,'EQSLV,SPARSE\n');
    fprintf(fileID,'bcsoption,,incoere\n');
    fprintf(fileID,'time,0.85\n');
    fprintf(fileID,'autot,on\n');

    fprintf(fileID,'nsubst,25,500,10\n');
    fprintf(fileID,'NEQIT,100\n');
    comman_constraints = ['DA,' num2str(area_idx_bottem) ',all,0\n'];
    fprintf(fileID,comman_constraints);

    fprintf(fileID,MODEL_PARAMETERS.CONSTRAINTS);

    % common constraints
    fprintf(fileID, 'DA,1002,uz,0\n DA,1005,uz,0\n DA,1008,uz,0\n DA,1011,uz,0\n DA,2002,uz,0\n DA,2005,uz,0\n DA,2008,uz,0\n DA,2011,uz,0\n ');


    fprintf(fileID,'/psf,pres,norm,2\n');


end

