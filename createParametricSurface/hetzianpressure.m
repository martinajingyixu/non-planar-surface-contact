function [p] = hetzianpressure(r,radius,pressureCenter)
% compute Hetzian pressure distribtuion between cylinder and elastical half
% space. 
%https://en.wikipedia.org/wiki/Contact_mechanics#Contact_between_a_rigid_cylinder_with_flat-ended_and_an_elastic_half-space
% r: the parametric surface
% pressureCenter is a line, [0,0,v]

syms u v p real
% make p0 = 1, then scale the pressure based on p0. its linear.
p0 = 1;
distQuad = sum((r-pressureCenter).^2);
% p = p0./sqrt((1 - distQuad./radius^2));
p = p0.*sqrt(1 - distQuad/(radius^2));

end

