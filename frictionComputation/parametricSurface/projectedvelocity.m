function [projectedVp] = projectedvelocity(vp,NtNormalized)
% compute the velocity projected to the tangent plane of P

%% use vp directly for projected vp computation
    dot_product = (dot(vp,NtNormalized));
    vpp = vp - dot_product / ((norm(NtNormalized))^2) * NtNormalized;
%     vpp_norm = (norm(vpp));
%     vppx = vpp(1);
%     vppy = vpp(2);
%     vppz = vpp(3);
%% First normalize vp, then use it for projected vp, and not normalize again
    vpNorm = vp ./ norm(vp);
    dot_product = (dot(vpNorm,NtNormalized));
    vppNormalizedVp = vpNorm - dot_product / ((norm(NtNormalized))^2) * NtNormalized;

%%
%     projectedVp.normalized = (vpp./vpp_norm);
%     projectedVp.notNormallized = vpp;

    projectedVp.normalizedVp = vppNormalizedVp;
    projectedVp.notNormallizedVp = vpp;
    
%     projectedVp.lowerBoundAbs = ...
%        (vpp./(abs(vppx)+abs(vppy)+abs(vppz)));

end

