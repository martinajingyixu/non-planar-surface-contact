function [outputArg1,outputArg2] = processansysdata(dirFiles,inputArg2)
%processansysdata Read the simulation results in txt to .mat files

%% for tests
dirFiles = '/home/jxu/promotion/programming/softContactGrasping//simulation/results/Ademo5//0001_bottle-gen_w5edel_150x035x020-03/xaxis//pressure0.0001/loc1_x35_y0_z0';
%% specify file names
fileDeformationUx = [dirFiles '/deformationUx.txt'];
fileDeformationUx = [dirFiles '/deformationUy.txt'];
fileDeformationUx = [dirFiles '/deformationUz.txt'];
fileElementComponents = [dirFiles '/elementComponents.txt'];
fileElementpressure = [dirFiles '/elementPressure.txt'];
fileNodalResults = [dirFiles '/nodalResults.txt'];
fileNodalCoordinates= [dirFiles '/nodalCoordinates.txt'];

end

