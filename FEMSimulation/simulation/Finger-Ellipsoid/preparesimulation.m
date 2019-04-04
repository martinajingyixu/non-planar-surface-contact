function [] = preparesimulation(filename)
    fileID = fopen(filename,'a');
    fprintf(fileID,'/UIS,MSGPOP,2 !! Disable error message pop up\n');
    fprintf(fileID,'!/UIS,MSGPOP,3 \n');
    fprintf(fileID,'KEYW,PR_SGVOF,1 \n');
    fprintf(fileID,'/NERR,5,10000, ,1,5, \n');
    fprintf(fileID,'!/UIS,MSGPOP,3  \n');
    fprintf(fileID,'KEYW,PR_SGVOF,1 \n');
    fprintf(fileID,'/NERR,5,10000, ,1,5, \n');

    fprintf(fileID,'FINISH\n');
    fprintf(fileID,'/CLEAR\n');
    fprintf(fileID,'/PREP7\n');

    fprintf(fileID,'/AUX15\n');
    fprintf(fileID,'IOPTN,IGES,SMOOTH\n');
    fprintf(fileID,'IOPTN,MERGE,YES\n');
    fprintf(fileID,'IOPTN,SOLID,YES\n');
    fprintf(fileID,'IOPTN,SMALL,YES\n');
    fprintf(fileID,'IOPTN,GTOLER, DEFA\n');
end