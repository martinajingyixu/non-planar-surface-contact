function [angle,coord] = samplesphere(numPoints)
%evenly sample a unit sphere. Input: number of point
% output, angle and xyz coordinate

% numPointsNotLimit = numPoints -6;
numPointsNotLimit = numPoints ;

indices = 0:(numPointsNotLimit-1) ;
indices = indices+ 0.5;

phi = asin(1 - 2.*indices./numPointsNotLimit);
theta = pi * (1 + sqrt(5)) .* indices;

% phi = [phi,0,0,0,0,pi/2,-pi/2];
% theta = [theta,0,pi,pi/2,-pi/2,0,0];
x = cos(theta) .* cos(phi);
y = sin(theta) .* cos(phi);
z = sin(phi);
% plot3(x,y,z,'.');
% 
phiLimit = [0,0,0,0,pi/2,-pi/2];
thetaLimit = [0,pi,pi/2,-pi/2,0,0];

xl = cos(thetaLimit) .* cos(phiLimit);
yl = sin(thetaLimit) .* cos(phiLimit);
zl = sin(phiLimit);

for iS = 1:length(xl)
    distances = sqrt(sum(bsxfun(@minus, [x;y;z], repmat([xl(iS);yl(iS);zl(iS)],[1,length(phi)])).^2,1));
    [~,minDistIdx] = min(distances);
    phi(minDistIdx) = [];
    theta(minDistIdx) = [];
    x(minDistIdx) = [];
    y(minDistIdx) = [];
    z(minDistIdx) = [];
end
phi = [phi,phiLimit];
theta = [theta,thetaLimit];
x = [x,xl];
y = [y,yl];
z = [z,zl];

angle = [phi;theta]';
coord = [x;y;z]';
end

