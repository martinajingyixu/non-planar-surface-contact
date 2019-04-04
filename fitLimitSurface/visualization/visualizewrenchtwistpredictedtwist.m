function visualizewrenchtwistpredictedtwist(wrench,figHandle, twist1,twist2,colorWrench,color1,color2,mkSize,linewidth)

if (nargin < 9)
    linewidth = 2;
end
if (nargin < 8)
    mkSize = 10;
end
if (nargin < 7)
    color1 = '[0.2078,0.9098,0]';
    color2 = '[0.5569,0.2118,0.9098]';
    colorWrench = '[1,0.7176,0.2706]';
end
if (nargin < 3)
    figure;
else
    hold on;
    figure(figHandle);
end


numF = size(wrench,2);
numV = size(twist1,2);

assert(numF == numV);
% plot3(wrench(1,:), wrench(2,:), wrench(3,:), 'Marker', 'o', 'MarkerFaceColor', 'r', 'Markersize', mkSize, 'LineStyle', 'none');
% plot3(wrench(1,:), wrench(2,:), wrench(3,:), 'Marker', 'o', 'MarkerFaceColor', '[0.9290, 0.6940, 0.1250]', 'Markersize', mkSize, 'LineStyle', 'none');
% plot3(wrench(1,:), wrench(2,:), wrench(3,:), 'Marker', 'o', 'MarkerFaceColor', '[1, 0.5216, 0.2824]', 'Markersize', mkSize, 'LineStyle', 'none');
plot3(wrench(1,:), wrench(2,:), wrench(3,:), 'Marker', 'o', 'MarkerFaceColor', colorWrench, 'Markersize', mkSize, 'LineStyle', 'none');


for i = 1:numF
    p0 = [wrench(1,i), wrench(2,i), wrench(3,i)];
    vnorm = norm([twist1(1,i), twist1(2,i), twist1(3,i)]);
    p1 = p0 + [twist1(1,i), twist1(2,i), twist1(3,i)]./vnorm.*0.5;
    vectarrow(p0,p1,color1,linewidth);
    hold on
    vnorm = norm([twist2(1,i), twist2(2,i), twist2(3,i)]);
    p1 = p0 + [twist2(1,i), twist2(2,i), twist2(3,i)]./vnorm.*0.45;
    vectarrow(p0,p1,color2,linewidth*1.1);
    hold on
end

end


