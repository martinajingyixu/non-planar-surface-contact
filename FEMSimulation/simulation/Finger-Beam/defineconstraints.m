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

    fprintf(fileID,MODEL_PARAMETERS.CONSTRAINTS1);
    fprintf(fileID,MODEL_PARAMETERS.CONSTRAINTS2);
    fprintf(fileID,MODEL_PARAMETERS.CONSTRAINTS3);
    fprintf(fileID,MODEL_PARAMETERS.CONSTRAINTS4);
    fprintf(fileID,'/psf,pres,norm,2\n');


end

