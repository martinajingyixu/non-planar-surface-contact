function [d] = perpendicular2vec(r,a,b,c,x0,y0,z0)
% compute d: the vector from the point r (pr called P) to another vector(e.g. the COR axis), d should be perpendicular to COR
% a,b,c,x0, y0,z0 define COR direction and location.
% sysmtem
xp = r(1);
yp = r(2);
zp = r(3);

temp_term = a*(xp - x0) + b*(yp - y0) + c*(zp - z0);
dx = (xp - x0 - a*temp_term);
dy = (yp - y0 - b*temp_term);
dz = (zp - z0 - c*temp_term);

d = [dx;dy;dz];

end

