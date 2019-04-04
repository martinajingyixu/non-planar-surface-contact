function [fxSyms,fySyms,fzSyms,tauxSyms,tauySyms,tauzSyms]...
    = computefrictionVppNormbased(pressureDist,projectedVp,r,origin,dA)

% the direction of velocity is same as friction. The velocity is the one
% from the object, relative to the velocity to the table. The friction is
% the force that the object gives to the table. So the have the same
% direction.
    projectedVpNormalized = projectedVp ./norm(projectedVp);

    fxSyms = projectedVpNormalized(1)*pressureDist*dA;
    fySyms = projectedVpNormalized(2)*pressureDist*dA;
    fzSyms = projectedVpNormalized(3)*pressureDist*dA;
    fSyms = [fxSyms;fySyms;fzSyms];

%     fxMinusSyms = -projectedVpNormalized(1)*pressureDist*dA;
%     fyMinusSyms = -projectedVpNormalized(2)*pressureDist*dA;
%     fzMinusSyms = -projectedVpNormalized(3)*pressureDist*dA;

    l = r - origin;
    %tauSyms = cross(fSyms,l);
    tauSyms = cross(l,fSyms);

    tauxSyms = tauSyms(1);
    tauySyms = tauSyms(2);
    tauzSyms = tauSyms(3);

%     tauxMinusSyms = -tauxSyms;
%     tauyMinusSyms = -tauySyms;
%     tauzMinusSyms = -tauzSyms;
end

