function plot2Dfrictioncomponents(xAxisValue,frictionWrench,color)

fx = frictionWrench(:,1);
fy = frictionWrench(:,2);
fz = frictionWrench(:,3);
taux = frictionWrench(:,4);
tauy = frictionWrench(:,5);
tauz = frictionWrench(:,6);

%%
figure
subplot(2,3,1)
hold on
plot(xAxisValue,fx,color);
title('fx')

subplot(2,3,2)
hold on
plot(xAxisValue,fy,color);
title('fy')


subplot(2,3,3)
plot(xAxisValue,fz,color);
title('fz')

subplot(2,3,4)
plot(xAxisValue,taux,color);
hold on
title('taux')

subplot(2,3,5)
plot(xAxisValue,tauy,color);
title('tauy')

subplot(2,3,6)
plot(xAxisValue,tauz,color);
hold on
title('tauz')

end

