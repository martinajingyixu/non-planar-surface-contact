function [FIntegrand,tauIntegrand] = computenormalforcetorque(pressureDist,NtNormalized,dA,r,origin)
%Compute the 6d force and torque of the pressure


FIntegrand = NtNormalized * dA * pressureDist;
l = r - origin;
tauIntegrand = cross(l,FIntegrand);

end

