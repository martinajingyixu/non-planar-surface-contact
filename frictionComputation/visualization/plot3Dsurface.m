function plot3Dsurface(fx,fy,fz,taux,tauy,tauz)

figure
subplot(2,3,1)
plot3(fx(:,1),fx(:,2),fx(:,3),'o');
% fxFit = fit([fx(:,1),fx(:,2)],fx(:,3),'rat23');
% plot(fxFit,[fx(:,1),fx(:,2)],fx(:,3))
hold on
title('fx')


subplot(2,3,2)
plot3(fy(:,1),fy(:,2),fy(:,3),'o');
hold on
title('fy')


subplot(2,3,3)
plot3(fz(:,1),fz(:,2),fz(:,3),'o');
hold on
title('fz')

subplot(2,3,4)
plot3(taux(:,1),taux(:,2),taux(:,3),'o');
hold on
title('taux')

subplot(2,3,5)
plot3(tauy(:,1),tauy(:,2),tauy(:,3),'o');
hold on
title('tauy')

subplot(2,3,6)
plot3(tauz(:,1),tauz(:,2),tauz(:,3),'o');
hold on
title('tauz')


end

