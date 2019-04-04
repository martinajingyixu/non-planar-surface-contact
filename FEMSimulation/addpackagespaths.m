function [] = addpackagespaths(package_root)
%addlibrariespaths Add path of each relavent libraties and packages
global CONSTANTS
addpath(package_root);
%% Add path of dataset and its subfolers
addpath(genpath([package_root '/dataset/']));
    
%% Add path of libs  
addpath(genpath([package_root '/libs/']));

%% Add path of self-implemented algorthims 
addpath(genpath([package_root '/frictionAnalysis/parametricSurface/']));
addpath([package_root '/simulation/' CONSTANTS.FINGER_TYPE '/']);
%% Add path of general functions
% addpath(genpath([package_root '/common/']));
addpath(genpath([package_root '/IO/']));

%% Add path of results and ouputs
addpath(genpath([package_root '/results/figures']));
addpath(genpath([package_root '/results/files']));

end

