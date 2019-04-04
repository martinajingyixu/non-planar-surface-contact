function computefrictionparfor(pathProcessingSurface,pathTwistSamples,idxStartSample,idxEndSample)
% function computefrictionparfor(pathProcessingSurface,dataTwist,dataSymbolics,idxStartSample,idxEndSample)

%Compute the integral of all samples by using parfor.
% the parametric form of the friction and the Twist samples are loaded from
% the .mat

dataTwist = load(pathTwistSamples);


allSamples = dataTwist.twistSampleVariables;

dataSymbolics = load([pathProcessingSurface '/FrictionParametricForm.mat']);
u = dataSymbolics.frictionParam.u;
v = dataSymbolics.frictionParam.v;

fx2intSyms = dataSymbolics.frictionParam.fx2intSyms;
fy2intSyms = dataSymbolics.frictionParam.fy2intSyms;
fz2intSyms = dataSymbolics.frictionParam.fz2intSyms;
taux2intSyms = dataSymbolics.frictionParam.taux2intSyms;
tauy2intSyms = dataSymbolics.frictionParam.tauy2intSyms;
tauz2intSyms = dataSymbolics.frictionParam.tauz2intSyms;


% vpx2int = dataSymbolics.velocityP.vpx2int;
% vpy2int = dataSymbolics.velocityP.vpy2int;
% vpz2int = dataSymbolics.velocityP.vpz2int;
% 
% vpNormalizedx2int = dataSymbolics.velocityP.vpNormalizedx2int;
% vpNormalizedy2int = dataSymbolics.velocityP.vpNormalizedy2int;
% vpNormalizedz2int = dataSymbolics.velocityP.vpNormalizedz2int;
% 
% vpProjectedx2int = dataSymbolics.frictionParam.velocityP.vpProjectedNotNormalizedVpx2int;
% vpProjectedy2int = dataSymbolics.frictionParam.velocityP.vpProjectedNotNormalizedVpy2int;
% vpProjectedz2int = dataSymbolics.frictionParam.velocityP.vpProjectedNotNormalizedVpz2int;
% 
% vpProjectedNormalizedx2int = dataSymbolics.frictionParam.velocityP.vpProjectedNormalizedVpx2int;
% vpProjectedNormalizedy2int = dataSymbolics.frictionParam.velocityP.vpProjectedNormalizedVpy2int;
% vpProjectedNormalizedz2int = dataSymbolics.frictionParam.velocityP.vpProjectedNormalizedVpz2int;
% 
% % area = dataSymbolics.surfaceArea;
% mass = dataSymbolics.surfaceInfo.mass;

%%
lowerBoundU = dataSymbolics.surfaceInfo.lowerBoundU;
upperBoundU = dataSymbolics.surfaceInfo.upperBoundU;

lowerBoundV = dataSymbolics.surfaceInfo.lowerBoundV;
upperBoundV = dataSymbolics.surfaceInfo.upperBoundV;

numCases = length(allSamples);
[fxArray,fyArray,fzArray,tauxArray,tauyArray,tauzArray] = deal(zeros(1,numCases));
% [vpx,vpy,vpz] = deal(zeros(numCases,1));
% [vpProjectedx,vpProjectedy,vpProjectedz] = deal(zeros(1,numCases));
% [vpx,vpy,vpz,vpNormalizedx,vpNormalizedy,vpNormalizedz] = deal(zeros(numCases,1));
% [vpProjectedx,vpProjectedy,vpProjectedz,...
%     vpProjectedNormalizedx,vpProjectedNormalizedy,vpProjectedNormalizedz] = deal(zeros(numCases,1));
% [vpProjectedx,vpProjectedy,vpProjectedz] = deal(zeros(numCases,1));
% [vpProjectedNormalizedx,vpProjectedNormalizedy,vpProjectedNormalizedz] = deal(zeros(numCases,1));

%[fxMinusArray,fyMinusArray,fzMinusArray,tauxMinusArray,tauyMinusArray,tauzMinusArray]...
%    = deal(zeros(numCases,1));

% for iSample = 1
parfor iSample = idxStartSample:idxEndSample
    iSample

    curSample = allSamples(:,iSample);

    aV= curSample(1);
    bV = curSample(2);
    cV = curSample(3);
    x0V= curSample(4);
    y0V= curSample(5);
    z0V= curSample(6);
    hV= curSample(7);

    fxArray(iSample) = integralcustimized(...
        (fx2intSyms(aV,bV,cV,x0V,y0V,z0V,hV,u,v))...
        ,lowerBoundU,upperBoundU,lowerBoundV,upperBoundV);
    fyArray(iSample) = integralcustimized(...
        (fy2intSyms(aV,bV,cV,x0V,y0V,z0V,hV,u,v))...
        ,lowerBoundU,upperBoundU,lowerBoundV,upperBoundV);
    fzArray(iSample) = integralcustimized(...
        (fz2intSyms(aV,bV,cV,x0V,y0V,z0V,hV,u,v))...
        ,lowerBoundU,upperBoundU,lowerBoundV,upperBoundV);
    tauxArray(iSample) = integralcustimized(...
        (taux2intSyms(aV,bV,cV,x0V,y0V,z0V,hV,u,v))...
        ,lowerBoundU,upperBoundU,lowerBoundV,upperBoundV);
    tauyArray(iSample) = integralcustimized(...
        (tauy2intSyms(aV,bV,cV,x0V,y0V,z0V,hV,u,v))...
        ,lowerBoundU,upperBoundU,lowerBoundV,upperBoundV);
    tauzArray(iSample) = integralcustimized(...
        (tauz2intSyms(aV,bV,cV,x0V,y0V,z0V,hV,u,v))...
        ,lowerBoundU,upperBoundU,lowerBoundV,upperBoundV);


%%
%     vpx(iSample) = integralcustimized(...
%         (vpx2int(aV,bV,cV,x0V,y0V,z0V,u,v))...
%         ,lowerBoundU,upperBoundU,lowerBoundV,upperBoundV)./mass;
% 
%     vpNormalizedx(iSample,1) = integralcustimized(...
%         (vpNormalizedx2int(aV,bV,cV,x0V,y0V,z0V,u,v))...
%         ,lowerBoundU,upperBoundU,lowerBoundV,upperBoundV)./mass;
% 
%     vpProjectedx(iSample) = integralcustimized(...
%         (vpProjectedx2int(aV,bV,cV,x0V,y0V,z0V,u,v))...
%         ,lowerBoundU,upperBoundU,lowerBoundV,upperBoundV)./mass;
% 
%     vpProjectedNormalizedx(iSample,1) = integralcustimized(...
%         (vpProjectedNormalizedx2int(aV,bV,cV,x0V,y0V,z0V,u,v))...
%         ,lowerBoundU,upperBoundU,lowerBoundV,upperBoundV)./mass;
% 
% 
% 
%     vpy(iSample) = integralcustimized(...
%         (vpy2int(aV,bV,cV,x0V,y0V,z0V,u,v))...
%         ,lowerBoundU,upperBoundU,lowerBoundV,upperBoundV)./mass;

%         vpNormalizedy(iSample,1) = integralcustimized(...
%             (vpNormalizedy2int(aV,bV,cV,x0V,y0V,z0V,u,v))...
%             ,lowerBoundU,upperBoundU,lowerBoundV,upperBoundV)./mass;
%         
%         vpProjectedy(iSample,1) = integralcustimized(...
%             (vpProjectedy2int(aV,bV,cV,x0V,y0V,z0V,u,v))...
%             ,lowerBoundU,upperBoundU,lowerBoundV,upperBoundV)./mass;
% %         
%         vpProjectedNormalizedy(iSample,1) = integralcustimized(...
%             (vpProjectedNormalizedy2int(aV,bV,cV,x0V,y0V,z0V,u,v))...
%             ,lowerBoundU,upperBoundU,lowerBoundV,upperBoundV)./mass;

% 
% 
%         vpz(iSample,1) = integralcustimized(...
%             (vpz2int(aV,bV,cV,x0V,y0V,z0V,u,v))...
%             ,lowerBoundU,upperBoundU,lowerBoundV,upperBoundV)./mass;
% 
%         vpNormalizedz(iSample,1) = integralcustimized(...
%             (vpNormalizedz2int(aV,bV,cV,x0V,y0V,z0V,u,v))...
%             ,lowerBoundU,upperBoundU,lowerBoundV,upperBoundV)./mass;
%         
%         vpProjectedz(iSample,1) = integralcustimized(...
%             (vpProjectedz2int(aV,bV,cV,x0V,y0V,z0V,u,v))...
%             ,lowerBoundU,upperBoundU,lowerBoundV,upperBoundV)./mass;
%         
%         vpProjectedNormalizedz(iSample,1) = integralcustimized(...
%             (vpProjectedNormalizedz2int(aV,bV,cV,x0V,y0V,z0V,u,v))...
%             ,lowerBoundU,upperBoundU,lowerBoundV,upperBoundV)./mass;

    
end
mat_name = [pathProcessingSurface '/friction_sphere_Twist_center_' num2str(idxStartSample) '_' num2str(idxEndSample) '.mat'];

save(mat_name,...
    'fxArray','fyArray','fzArray',...
    'tauxArray','tauyArray','tauzArray',...
    'idxStartSample','idxEndSample');

% mat_name_vel = [pathProcessingSurface '/velocity_sphere_Twist_center_' num2str(idxStartSample) '_' num2str(idxEndSample) '.mat'];
% save(mat_name_vel,'vpProjectedx','vpProjectedy','vpProjectedz');

% mat_name_vel = [pathProcessingSurface '/velocity_sphere_Twist_center_' num2str(idxStartSample) '_' num2str(idxEndSample) '.mat'];
% save(mat_name_vel,...
%     'vpx','vpy','vpz',...
%     'vpNormalizedx','vpNormalizedy','vpNormalizedz',...
%     'vpProjectedx','vpProjectedy','vpProjectedz',...
%     'vpProjectedNormalizedx','vpProjectedNormalizedy','vpProjectedNormalizedz');
clearvars -global
clearvars
clear functions
clear all
end

