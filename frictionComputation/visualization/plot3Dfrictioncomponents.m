function plot3Dfrictioncomponents(xAxisValue,yAxisValue,frictionWrench)

fx = frictionWrench(:,1);
fy = frictionWrench(:,2);
fz = frictionWrench(:,3);
taux = frictionWrench(:,4);
tauy = frictionWrench(:,5);
tauz = frictionWrench(:,6);

figure
subplot(2,3,1)
plot3(xAxisValue,yAxisValue,fx,'.');
hold on
title('fx')


subplot(2,3,2)
plot3(xAxisValue,yAxisValue,fy,'.');
hold on
title('fy')


subplot(2,3,3)
plot3(xAxisValue,yAxisValue,fz,'.');
hold on
title('fz')

subplot(2,3,4)
plot3(xAxisValue,yAxisValue,taux,'.');
hold on
title('taux')

subplot(2,3,5)
plot3(xAxisValue,yAxisValue,tauy,'.');
hold on
title('tauy')

subplot(2,3,6)
plot3(xAxisValue,yAxisValue,tauz,'.');
hold on
title('tauz')


end

