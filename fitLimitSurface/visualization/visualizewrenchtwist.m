function visualizewrenchtwist(wrench, twist,figHandle, color,mkSize,linewidth)

if (nargin < 6)
    linewidth = 2;
end
if (nargin < 5)
    mkSize = 10;
end
if (nargin < 4)
    color = 'b';
end
if (nargin < 3)
    figure;
else
    hold on;
    figure(figHandle);
end


numF = size(wrench,2);
numV = size(twist,2);

assert(numF == numV);
% plot3(wrench(1,:), wrench(2,:), wrench(3,:), 'Marker', 'o', 'MarkerFaceColor', 'r', 'Markersize', mkSize, 'LineStyle', 'none');
% plot3(wrench(1,:), wrench(2,:), wrench(3,:), 'Marker', 'o', 'MarkerFaceColor', '[0.9290, 0.6940, 0.1250]', 'Markersize', mkSize, 'LineStyle', 'none');
% plot3(wrench(1,:), wrench(2,:), wrench(3,:), 'Marker', 'o', 'MarkerFaceColor', '[1, 0.5216, 0.2824]', 'Markersize', mkSize, 'LineStyle', 'none');
plot3(wrench(1,:), wrench(2,:), wrench(3,:), 'Marker', 'o', 'MarkerFaceColor', '[1, 0.4, 0.2824]', 'Markersize', mkSize, 'LineStyle', 'none');

hold on;
color2 = '[0.6627,0.2118,0.9098]';
vectarrow([0;0;0],[0;0;0.00001],color2,linewidth);
hold on;
for i = 1:numF
    p0 = [wrench(1,i), wrench(2,i), wrench(3,i)];
    vnorm = norm([twist(1,i), twist(2,i), twist(3,i)]).*5;
    p1 = p0 + [twist(1,i), twist(2,i), twist(3,i)]./vnorm;
    vectarrow(p0,p1,color,linewidth);
    hold on
end

end


