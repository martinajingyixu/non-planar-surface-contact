function [ifSuccuess] = runansysfingerellipsoid(dirIGESModel,simulationResultPath,objectInfo,MODEL,locationPos,locationNeg,pressure)
%runansys Run simulation with ansys. For one squeeze location and one
%pressure
global CONSTANTS

filename = [simulationResultPath '/ansys_code.txt'];
fileID = fopen(filename,'w');
fileID = fopen(filename,'a');

%% write ansys code
preparesimulation(filename);

commandLoadIGES = ['IGESIN,' char(39) objectInfo.filename char(39) ',' char(39) 'igs' char(39) ',' char(39) dirIGESModel char(39) '\n'];
fprintf(fileID,commandLoadIGES);

%Merge keypoints for AUX15 data (iges etc)

fprintf(fileID,'/PREP7\n');
fprintf(fileID,'alls\n');
fprintf(fileID,'nummrg,kp,,0.01\n');

%Create Bottem of the object
fprintf(fileID,'LSEL,s,loc,z,0\n');
fprintf(fileID,'AL,all\n');
fprintf(fileID,'ALLSEL\n');
fprintf(fileID,'NUMCMP,ALL\n');

%% Create VS sensor 
% specify two keypoints
formatSpec_KP = 'K,%d,%f,%f,%f\n';
sqeezeLocation1 = locationPos - 0.5; 
sqeezeLocation2 = locationNeg - 0.5;

fprintf(fileID,formatSpec_KP,1000,sqeezeLocation1(1),sqeezeLocation1(2),sqeezeLocation1(3));
fprintf(fileID,formatSpec_KP,1001,sqeezeLocation2(1),sqeezeLocation2(2),sqeezeLocation2(3));

% move workspace to squeeze location 1
fprintf(fileID,'KWPAVE,1000\n');
fprintf(fileID,MODEL.ROTATE_WORKSPACE);
fprintf(fileID,MODEL.MOVE_WORKSPACE_FINGER1,CONSTANTS.SIMULATION.FINGER_THICK,CONSTANTS.SIMULATION.SENSOR_HEIGHT);
fprintf(fileID,'CSWPLA,11,0,1,1 \n');

% create the finger 1, consists of 4 1/8 ellpsoid
fprintf(fileID,'NUMSTR,AREA,1000 \n'); %Area Rectangle for foam 1
fprintf(fileID,'NUMSTR,KP,2000 \n'); %Area Rectangle for foam 1
fprintf(fileID,'cyl4,,,,,%f,90\n',CONSTANTS.SIMULATION.FINGER_THICK);
fprintf(fileID,MODEL.VROTATE_FINGER1);
fprintf(fileID,MODEL.VSCALE_FINGER1);


fprintf(fileID,'KWPAVE,1001\n');
fprintf(fileID,MODEL.ROTATE_WORKSPACE);
% fprintf(fileID,MODEL.MOVE_WORKSPACE_NEGATIVE,CONSTANTS.SIMULATION.FINGER_THICK);
fprintf(fileID,MODEL.MOVE_WORKSPACE_FINGER2,CONSTANTS.SIMULATION.FINGER_THICK,CONSTANTS.SIMULATION.SENSOR_HEIGHT);

fprintf(fileID,'CSWPLA,11,0,1,1 \n');
fprintf(fileID,'NUMSTR,AREA,2000 \n'); %Area Rectangle for foam 1
fprintf(fileID,'NUMSTR,KP,3000 \n'); %Area Rectangle for foam 1
fprintf(fileID,'cyl4,,,,,%f,90\n',CONSTANTS.SIMULATION.FINGER_THICK);
fprintf(fileID,MODEL.VROTATE_FINGER2);
fprintf(fileID,MODEL.VSCALE_FINGER2);

fprintf(fileID,'VGLUE,ALL\n');


%%
if objectInfo.rim == 0
    thickRim = objectInfo.thick;
    objHeight = objectInfo.h;
else
    thickRim = CONSTANTS.SIMULATION.RIM_THICKNESS_FACTOR * objectInfo.thick;
    objHeight = objectInfo.h+objectInfo.res;
end

fprintf(fileID,'ET,1,SHELL181\n');
fprintf(fileID,'SECTYPE,1,SHELL\n');
command_thickness = ['SECDATA,' num2str(objectInfo.thick) '\n'];
fprintf(fileID,command_thickness);

fprintf(fileID,'ET,2,SHELL181\n');
fprintf(fileID,'SECTYPE,2,SHELL\n');
command_thickness2 = ['SECDATA,' num2str(thickRim) '\n'];
fprintf(fileID,command_thickness2);


%% Define material for object, sensor, steel (to prevent undisired deformation)
definematerial(filename);

%% Mesh
numWallArea = ceil(objHeight/objectInfo.res)*4; 
idxBottemArea = numWallArea+1;
meshobject(filename,numWallArea,idxBottemArea);

%% Define Element Type
defineelementtype(filename);

%% Define contacts and fem constraints
% sqeezeLocation1(1),sqeezeLocation1(2)
contactDist1 = MODEL.transToCenter * sqeezeLocation1;
contactDist2 = MODEL.transToCenter * sqeezeLocation2;

definecontact(filename,MODEL,contactDist1,contactDist2);
defineconstraints(filename,MODEL,idxBottemArea)

%% Apply pressure

fprintf(fileID,['SFA,1002,,pres,' num2str(pressure) '\n']); % Finger 1
fprintf(fileID,['SFA,1005,,pres,' num2str(pressure) '\n']);
fprintf(fileID,['SFA,1008,,pres,' num2str(pressure) '\n']);
fprintf(fileID,['SFA,1011,,pres,' num2str(pressure) '\n']);

fprintf(fileID,['SFA,2002,,pres,' num2str(pressure) '\n']); % Finger 2
fprintf(fileID,['SFA,2005,,pres,' num2str(pressure) '\n']);
fprintf(fileID,['SFA,2008,,pres,' num2str(pressure) '\n']);
fprintf(fileID,['SFA,2011,,pres,' num2str(pressure) '\n']);


fprintf(fileID,'solve\n');
fprintf(fileID,'finish/post1\n');

fprintf(fileID,'set, last\n');
fprintf(fileID,'/post1\n');
fprintf(fileID,'set, last\n');

deleteredundantfiles();

%% Store files and delete undesireble files
% will create useful .mac file here
storefiles(filename,simulationResultPath);

if CONSTANTS.SIMULATION.SAVE_DB_FILE
    save_command = ['save,' objectInfo.filename ',db,' simulationResultPath '\n'];
    fprintf(fileID,save_command);
end

commandRunAnsys = ['KMP_STACKSIZE=' CONSTANTS.SIMULATION.MAX_STACKSIZE ';find ~/ -name “*.lock” -exec rm {} \;  ANSYS_LOCK=OFF; ' CONSTANTS.DIR.ANSYS ' -np 6 -b -i ' simulationResultPath '/ansys_code.txt ' ' -o ' simulationResultPath '/simulationoutput.txt;'];
% commandRunAnsys = ['KMP_STACKSIZE=' CONSTANTS.SIMULATION.MAX_STACKSIZE ';  ANSYS_LOCK=OFF;' CONSTANTS.DIR.ANSYS ' -np 6 -b -i ' simulationResultPath '/ansys_code.txt ' ' -o ' simulationResultPath '/simulationoutput.txt;'];

% prevent system output
[ status, cmdout ] = system(commandRunAnsys);

ifSuccuess = issimulationsuccess(simulationResultPath);

fclose('all');
deleteredundantfiles();

end


function deleteredundantfiles()
    if exist([pwd '/file.err'],'file')
        delete ('file.*');
    end
    if ~isempty(dir('*.mac'))   
        delete ('*.mac');
    end
end



