function [NtNormalized] = computesurfacenormal(r,var1,var2)
%Compute the normalized surface normal of the surface 
% the Nt normalized is the norm of the tangent plane

% compute the normal of the tangent plane of P: Nt
% r_theta: differential of r on theta
% r_phi: differential of r on phi

% substitue variables to char to solve diff.
var1ToString = char(var1);
var2ToString = char(var2);

r_subs_var1 = subs(r,var1,var1ToString);
r_subs_var1_var2 = subs(r_subs_var1,var2,var2ToString);

r_var1_diff = diff(r_subs_var1_var2,var1ToString);
r_var2_diff =  diff(r_subs_var1_var2,var2ToString);

% substitue back: from char to var
r_var1_diff_subs = subs(subs(r_var1_diff,var1ToString,var1),var2ToString,var2);
r_var2_diff_subs = subs(subs(r_var2_diff,var1ToString,var1),var2ToString,var2);


Nt = cross(r_var1_diff_subs,r_var2_diff_subs);% Nt_norm = simplify(Nt/norm(Nt));
NtNorm = simplify(norm(Nt));
NtNormalized = simplify(Nt./NtNorm);

end

