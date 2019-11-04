% This script create object and finger with NURBS description and export to IGES files. 
% It calls ANSYS for simulation. ANSYS results are written in .csv files. 

% write the absolute path of this package, because ANSYS requires an
% absolute path
% PACKAGE_ROOT = '/home/xxx/xxx/xxx/3d-surface-contact/';
PACKAGE_ROOT = pwd; 
% A struct CONSTS constists of constants for:
% GRASP, CONTACT, DEBUG, GEOMETRY, MATERIAL

% setelct a finger type for simulation
% loadconstantsFingerBeam; 
loadconstantsFingerEllipsoid; % setup ANSYS directory in the config file
addpackagespaths(PACKAGE_ROOT);
[dbDir,objectParameters] = db_parameters(PACKAGE_ROOT);
objectRange = 1:height(objectParameters);
%% model generation
[dbDir,objectParameters] = db_parameters(PACKAGE_ROOT);
db_build(dbDir, objectParameters);
%% FEM simulation
db_simulation(dbDir,objectParameters,objectRange);
%% Check simulation success
dbSimulationSuccess = checksimulationsuccess(dbDir,objectParameters,objectRange);

%% read ansys data
% convert .csv files to .mat
FEMSimulationResults = exportdatatodb(dbDir,objectParameters,objectRange);
