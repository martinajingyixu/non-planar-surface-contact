function [vp] = computevelocity(r,a,b,c,x0,y0,z0,h)

% h is the pitch. The body twist consists of the parallel velocity to the omega 
% and the perpendicular to omega

% r: a point on the surface. 
% [x0;y0;z0]: a point on the COR axis
omega = [a;b;c];

% parallel velocity v0
v0 = h .* omega;
d = r - [x0;y0;z0];
% vp= v0 + cross(d,omega);
vp= v0 + cross(omega,d);

end

