%% create a multi-contact from parametric surface
surfaceTypeArray = {'cylinder','ellipticalCylinder','sphere','ellipsoid','paraboloid','ellipticalParaboloid'};
pressureTypeArray = {'uniform','hetzCenter'};
%% config
mu = 0.333; % only included in GWS computation
config.ifVisualize = 0;
contact.dist = 1;
iS = 4;
iP = 1; % a elliptic cylinder with a uniform pressure distribution

surfaceType = surfaceTypeArray{iS};
pressureType = pressureTypeArray{iP};
pathResults = ['frictionResults/discretizedSurface/'  surfaceType  '/'  pressureType  '/'];
load([pathResults 'frictionContactResults.mat'],...
    'elementsInfo','surfaceInfo','fittingResults','wrenchTwist');
%% create left and right contacts from parametric surfaces
% transform the surface and friction wrench. But friction torque is still
% with respect to the friction center
wrench = wrenchTwist.wrench;
fittedResEllip = fittingResults.fittedResEllip;
[contactL,contactR,transMat] = ...
    createcontactparallelgripper(contact.dist,surfaceInfo,...
    (fittedResEllip.sampledWrench' .* max(wrench'))',config.ifVisualize);
%% scale the magnitude of the normal force to one. Normal wrench and friction wrench w.r.t to the origin (COM)
[GWSMink,GWSUnion,wrenchL,wrenchR] = computegws(contactL,contactR,mu);
%% check whether a point is inside the convex hull. If its in, then its stable.
disturbance = [0.6;0.1;0.1;0.00;0.001;0];
ifInhull = inhull(disturbance',GWSUnion');
%% Visualization
% Mink will be very slow
% [~,hullVolumeMink] = convhulln(GWSMink');
% plot3dconvechullof6dwrenches((GWSMink'./max(GWSMink'))');

[~,hullVolumeUnion] = convhulln(GWSUnion');
plot3dconvechullof6dwrenches((GWSUnion'./max(abs(GWSUnion)'))');

plot6dwrenches((GWSMink'./max(GWSMink'))',(GWSUnion'./max(abs(GWSUnion)'))')



