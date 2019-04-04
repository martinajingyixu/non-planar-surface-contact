function surfaceInfo = definepressureparametricsurface(pathSurface,surfaceInfo,ifVisualize)
if nargin<3
    ifVisualize = 0;
end
syms u v r x y z ox oy oz real % u is Theta, v is Phi

lowerBoundU = surfaceInfo.lowerBoundU;
upperBoundU = surfaceInfo.upperBoundU;
lowerBoundV = surfaceInfo.lowerBoundV;
upperBoundV = surfaceInfo.upperBoundV;
r = surfaceInfo.r;
rToSample = surfaceInfo.rToSample;
coordCenterInt = surfaceInfo.coordCenterInt; % the pre-integrated com in mathematica
pressureType = surfaceInfo.pressureType;
surfaceType = surfaceInfo.surfaceType;

r1 = surfaceInfo.r1;
r2 = surfaceInfo.r2;
r3 = surfaceInfo.r3;

samplingRangeU =(upperBoundU-lowerBoundU)/30;
samplingRangeV = (upperBoundV-lowerBoundV)/30;

samplesU = lowerBoundU:samplingRangeU:upperBoundU;
samplesV = lowerBoundV:samplingRangeV:upperBoundV;

if lowerBoundU < 0 && upperBoundU > 0 
    samplesU = sort([lowerBoundU:samplingRangeU:upperBoundU,0]);
end
if lowerBoundV < 0 &&  upperBoundV > 0
    samplesV = sort([lowerBoundV:samplingRangeV:upperBoundV,0]);
end
rSampledCatesianCoord = parametric2cartesiansamples...
    (rToSample,samplesU,samplesV);
[rSampledCatesianCoord,downSampledIdx] = downsample3dpoints(rSampledCatesianCoord);

samlingCombi = combvec(samplesV,samplesU);
samlingCombi = samlingCombi(:,downSampledIdx); 

if isempty(coordCenterInt) % no pre-integrated coordCenter
    coordCenter = mean(rSampledCatesianCoord')';
else
    coordCenter = coordCenterInt;
end

if strcmp(pressureType,'uniform')
    pressureDist = 1;   
    pressureDistSyms = 1;
else
    pressureCenter = coordCenter;
    radiusPressure = max([r1,r2,r3])*1.8;
    
    if strcmp(surfaceType,'cylinder') ||  strcmp(surfaceType,'ellipticalCylinder')
        pressureCenter = [pressureCenter(1:2);v];
        pressureCenterSyms = [ox;oy;v];
    else
        pressureCenterSyms = [ox;oy;oz];
    end
    % radius of pressure is the longest distance from pressure center to
    % all points
    pressureDistSyms = hetzianpressure(r,radiusPressure,pressureCenterSyms);
    pressureDist = hetzianpressure(r,radiusPressure,pressureCenter);
end
surfaceInfo.pressureDist = pressureDist;
surfaceInfo.pressureDistSyms = pressureDistSyms;
surfaceInfo.rSampledCatesianCoord = rSampledCatesianCoord;
surfaceInfo.lengthArray = max(rSampledCatesianCoord')' - min(rSampledCatesianCoord')';
surfaceInfo.samlingCombi = samlingCombi;
save([pathSurface  'varsParametricForm.mat'],'surfaceInfo');

if ifVisualize
    figure
    syms p2plot(u,v)
    p2plot(u,v) = pressureDist;
    pressure = double(vpa(p2plot(samlingCombi(2,:),samlingCombi(1,:))));
    C = repmat([100],[1,numel(rSampledCatesianCoord(1,:))]);
    scatter3(rSampledCatesianCoord(1,:),rSampledCatesianCoord(2,:),...
        rSampledCatesianCoord(3,:),C,pressure,'filled')
    hold on
    axis equal
    colorbar
    hold on
    xlabel('x') % x-axis label
    ylabel('y') % y-axis label
    zlabel('z') % z-axis label
end
end

