function [vecCOR2Origin] = computevecCOR2Origin(Origin,a,b,c,x0,y0,z0)
% Compute the distance from origin to COR
% the orientation of the COR is defined by [a;b;c]
% COR goes through point x0,y0,z0
% the vector COR2Origin is perpendicular to COR. So
% dot(vecCOR2Origin.[a;b;c]) = 0

% origin is the friction center
Ox = Origin(1);
Oy = Origin(2);
Oz = Origin(3);

temp_term = a*(Ox - x0) + b*(Oy - y0) + c*(Oz - z0);
vecx = simplify(Ox - x0 - a*temp_term);
vecy = simplify(Oy - y0 - b*temp_term);
vecz = simplify(Oz - z0 - c*temp_term);

vecCOR2Origin = -[vecx;vecy;vecz];
% assert(double(simplify(dot(vecCOR2Origin,[a;b;c]))) == 0.0,...
%     'vecCOR2Origin should be perpendicular to COR axis');
end

