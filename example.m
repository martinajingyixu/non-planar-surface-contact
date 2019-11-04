%% example for FEM simulations
% exampleFEMSimulation

%% Example code to model non-planar surface contacts and its application in grasping. 
% The pipeline includes:
% 1. Contact profile 
% 2. Friction computation
% 3. Fit limit surfaces
% 4. Multi-contacts and applications in grasping 

% add paths
addpath(genpath('createParametricSurface'));
addpath(genpath('fitLimitSurface'));
addpath(genpath('frictionComputation'));
addpath(genpath('IO'));

%%  Step 1-3
config.ifVisualize = 0;
% examle for a meshed surface
exampleMeshedSurface; 

% example for a parametric surface
% exampleParametricSurface

% example for a contact from FEM simulation
% exampleFEM

%% Step 4: Example for multi contacts
addpath(genpath('multiContact'));
exampleMultiContact;

