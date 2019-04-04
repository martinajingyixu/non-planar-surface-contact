function [] = defineelementtype(filename)
    fileID = fopen(filename,'a');
    fprintf(fileID,'\n');
    

    
    fprintf(fileID,'ET,3,solid187\n');
    fprintf(fileID,'ET,10,hyper86\n');
    fprintf(fileID,'keyopt,10,2,1\n');

    fprintf(fileID,'\n'); % MP, EX 2 is the material for steel

    fprintf(fileID,'\n');
    fprintf(fileID,'TB,HYPER,3,,,BLATZ\n');
    fprintf(fileID,'TBDATA,1,mu\n');




end