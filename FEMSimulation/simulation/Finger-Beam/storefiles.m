function [] = storefiles(filename,simulationResultPath)
%storefiles Store simulation results

    global CONSTANTS
    fileID = fopen(filename,'a');
    fprintf(fileID,'/post1\n');
    fprintf(fileID,'set, last \n');
    fprintf(fileID,'/post1\n');
    fprintf(fileID,'set, last\n ');
    
%        
%     fprintf(fileID,['*create, '   'cubefingerelemcont,mac ' '\n']);
%     fprintf(fileID,'/post1\n');
%     fprintf(fileID,['/output, elementComponents,txt,',simulationResultPath,'\n']);
%     fprintf(fileID,'presol,cont\n');
%     fprintf(fileID,'/output\n');
%     fprintf(fileID,'*end\n');
%     fprintf(fileID,'cubefingerelemcont\n');
%     
    fprintf(fileID,['*create,  '   'cubefingerelempres,mac'  '\n ']);
    fprintf(fileID,'/post1\n');
    fprintf(fileID,['/output, elementPressure,txt,',simulationResultPath,'\n']);
    fprintf(fileID,'etable,ElePres,cont,pres\n');
    fprintf(fileID,'pretab,ELEPRES\n');
    fprintf(fileID,'/output\n');
    fprintf(fileID,'*end\n');
    fprintf(fileID, 'cubefingerelempres\n');
    
    createmacfile('exportelementareatemp', [CONSTANTS.DIR.SIMULATION_ROOT '/exportelementarea'], simulationResultPath, 'element_area');
    fprintf(fileID,['/INPUT, exportelementareatemp,mac,' pwd ' \n']);

    createmacfile('exportelementinfotemp', [CONSTANTS.DIR.SIMULATION_ROOT '/exportelementinfo'], simulationResultPath, 'element_info');
    fprintf(fileID,['/INPUT, exportelementinfotemp,mac,' pwd ' \n']);
    
    createmacfile('exportelementkinematicenergytemp', [CONSTANTS.DIR.SIMULATION_ROOT '/exportelementkinematicenergy'], simulationResultPath, 'element_kinematic_energy');
    fprintf(fileID,['/INPUT, exportelementkinematicenergytemp,mac,' pwd ' \n']);
    
    createmacfile('exportelementmattemp', [CONSTANTS.DIR.SIMULATION_ROOT '/exportelementmat'], simulationResultPath, 'element_mat');
    fprintf(fileID,['/INPUT, exportelementmattemp,mac,' pwd ' \n']);
    
    createmacfile('exportelementrealtemp', [CONSTANTS.DIR.SIMULATION_ROOT '/exportelementreal'], simulationResultPath, 'element_real');
    fprintf(fileID,['/INPUT, exportelementrealtemp,mac,' pwd ' \n']);
    
    createmacfile('exportelementstiffnessenergytemp', [CONSTANTS.DIR.SIMULATION_ROOT '/exportelementkinematicenergy'], simulationResultPath, 'element_stiffness_energy');
    fprintf(fileID,['/INPUT, exportelementstiffnessenergytemp,mac,' pwd ' \n']);
    
    createmacfile('exportnodalcoordstemp', [CONSTANTS.DIR.SIMULATION_ROOT '/exportnodalcoords'], simulationResultPath, 'nodal_coord');
    fprintf(fileID,['/INPUT, exportnodalcoordstemp,mac,' pwd ' \n']);
    
    createmacfile('exportnodaldisplacementtemp', [CONSTANTS.DIR.SIMULATION_ROOT '/exportnodaldisplacement'], simulationResultPath, 'nodal_displacement');
    fprintf(fileID,['/INPUT, exportnodaldisplacementtemp,mac,' pwd ' \n']);
    
%     fprintf(fileID,['*create, '   'cubefingernodecont,mac ' '\n ']);
%     fprintf(fileID,'/post1\n');
%     fprintf(fileID,['/output, nodalResults,txt,',simulationResultPath,'\n']);
%     fprintf(fileID,'prnsol,cont\n');
%     fprintf(fileID,'/output\n');
%     fprintf(fileID,'*end\n');
%     fprintf(fileID,'cubefingernodecont\n');
%     
%     fprintf(fileID,['*create, '   'cubefingernodecoord,mac'  ' \n']);
%     fprintf(fileID,'/post1\n');
%     fprintf(fileID,['/output, nodeCoordinates,txt,',simulationResultPath,'\n']);
%     fprintf(fileID,'NLIST,all,COORD,X,Y,Z\n');
%     fprintf(fileID,'/output\n');
%     fprintf(fileID,'*end\n');
%     fprintf(fileID,'cubefingernodecoord\n');
%     
%     fprintf(fileID,['*create, '   'cubefingerux,mac '  ' \n']);
%     fprintf(fileID,'/post1\n');
%     fprintf(fileID,['/output, deformationUx,txt,',simulationResultPath,'\n']);
%     fprintf(fileID,'prnsol,u,x\n');
%     fprintf(fileID,'/output\n');
%     fprintf(fileID,'*end\n');
%     fprintf(fileID,'cubefingerux\n');
%     
%     fprintf(fileID,['*create, '   'cubefingeruy,mac '  '\n ']);
%     fprintf(fileID,'/post1\n');
%     fprintf(fileID,['/output, deformationUy,txt,',simulationResultPath,'\n']);
%     fprintf(fileID,'prnsol,u,y\n');
%     fprintf(fileID,'/output\n');
%     fprintf(fileID,'*end\n');
%     fprintf(fileID,'cubefingeruy\n');
% % 
%     fprintf(fileID,['*create,'   'cubefingeruz,mac' '\n ']);
%     fprintf(fileID,'/post1\n');
%     fprintf(fileID,['/output, deformationUz,txt,',simulationResultPath,'\n']);
%     fprintf(fileID,'prnsol,u,z\n');
%     fprintf(fileID,'/output\n');
%     fprintf(fileID,'*end\n');
%     fprintf(fileID,'cubefingeruz\n');  
end

function [] = createmacfile(macFileName, txtBaseFileName, simulationResultPath, targetFileName)
    filename = [ macFileName '.mac'];
    command_merge_mac = ['touch ' filename];
    [~,~] = system(command_merge_mac);

    fileID = fopen(filename,'w');
    fprintf(fileID,'alls\n');
    command_add_target_file = ['*CFOPEN, ' targetFileName ',csv,'  simulationResultPath ' \n'];
    fprintf(fileID,command_add_target_file);
    command_add_base_file = ['cat ' txtBaseFileName '.txt  >> ' macFileName '.mac'];
    [~,~] = system(command_add_base_file);
end
