%% Find idx and point that are interesting
fxP = abs(fx);
fyP = abs(fy);
fzP = abs(fz);
tauxP = abs(taux);
tauyP = abs(tauy);
tauzP = abs(tauz);

%%
idx = 1:numel(fx)*2;
varC1 = 2;
varC2 = 3;
%%
CORidx = find(CORSamples(:,1) == 0.5*pi & CORSamples(:,2) == 0);
idx = CORidx;
CORtest= CORSamples4Analysis(idx,:);
varC1 = 4;
varC2 = 5;
frictionWrench = [fx(idx),fy(idx),fz(idx),taux(idx),tauy(idx),tauz(idx)];
xAxisValue = CORSamples4Analysis(idx,varC1);
yAxisValue = CORSamples4Analysis(idx,varC2);
plot3Dfrictioncomponents(xAxisValue,yAxisValue,frictionWrench)
%%
i = 1;
idx = 468*(i-1)+1:468*i;
xAxisValue = 1:numel(idx);
%%
CORidx = find(CORSamples(:,1) == 0 & CORSamples(:,2) == 0);
idx = CORidx;
idx = idx(79:91);
xAxisValue = CORSamples4Analysis(idx,5)
frictionWrenchP = [fxP(idx),fyP(idx),fzP(idx),tauxP(idx),tauyP(idx),tauzP(idx)];

%%
idx = find(abs(COR(:,1))<1e-5 &COR(:,3)==0&COR(:,4)==0);
xAxisValue = COR(idx,6);
% plot3Dfrictioncomponents(xAxisValue,yAxisValue,frictionWrench)
frictionWrench = [fxM,fyM,fzM,tauxM,tauyM,tauzM];
% frictionWrenchLBABS = [fxMLBABS(idx),fyMLBABS(idx),fzMLBABS(idx),...
%     tauxMLBABS(idx),tauyMLBABS(idx),tauzMLBABS(idx)];

plot2Dfrictioncomponents(1:numel(fxM),frictionWrench,'b')
%%
plot(tauxM3D,tauyM3D,'.')
axis equal
hold on
plot(tauxM,tauyM,'r.')


%%
% COR, xaxis, z0 = 0,y0 = something
idx = find(CORSamples4Integral(:,2) ==0 &CORSamples4Integral(:,3) ==0&...
CORSamples4Integral(:,4) ==0 & abs(CORSamples4Integral(:,6))<1e-1);
frictionWrench = [fx(idx),fy(idx),fz(idx),...
    taux(idx),tauy(idx),tauz(idx)];
plot2Dfrictioncomponents(CORSamples4Integral(idx,5),frictionWrench,'b.')
plot2Dfrictioncomponents(CORSamples4Integral(idx,5),frictionWrench,'b')


%%
idx = 1:numel(fx)*2;
frictionWrench = [fxM(idx),fyM(idx),fzM(idx),tauxM(idx),tauyM(idx),tauzM(idx)];
xAxisValue = 1:numel(idx);
plot2Dfrictioncomponents(xAxisValue,frictionWrench)
