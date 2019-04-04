function [velocityP] = velocityIntegrand(vp,projectedVp,symbols,dA,rho,thick)
%Prepare the velocity to be integrated

a = symbols.a;
b = symbols.b;
c = symbols.c;
x0 = symbols.x0;
y0 = symbols.y0;
z0 = symbols.z0;
u = symbols.u;
v = symbols.v;
h = symbols.h;

vpNormalized = vp./norm(vp);
vpProjectedVpNormalized = projectedVp.normalizedVp;
vpProjectedVpNotNormalized= projectedVp.notNormallizedVp;
%%

dMass = dA*rho*thick;

velocityP.vpx2int(a,b,c,x0,y0,z0,h,u,v) = simplify(vp(1))*dMass;
velocityP.vpy2int(a,b,c,x0,y0,z0,h,u,v) = simplify(vp(2))*dMass;
velocityP.vpz2int(a,b,c,x0,y0,z0,h,u,v) = simplify(vp(3))*dMass;

velocityP.vpNormalizedx2int(a,b,c,x0,y0,z0,h,u,v) = simplify(vpNormalized(1))*dMass;
velocityP.vpNormalizedy2int(a,b,c,x0,y0,z0,h,u,v) = simplify(vpNormalized(2))*dMass;
velocityP.vpNormalizedz2int(a,b,c,x0,y0,z0,h,u,v) = simplify(vpNormalized(3))*dMass;

velocityP.vpProjectedNormalizedVpx2int(a,b,c,x0,y0,z0,h,u,v) = simplify(vpProjectedVpNormalized(1))*dMass;
velocityP.vpProjectedNormalizedVpy2int(a,b,c,x0,y0,z0,h,u,v) = simplify(vpProjectedVpNormalized(2))*dMass;
velocityP.vpProjectedNormalizedVpz2int(a,b,c,x0,y0,z0,h,u,v) = simplify(vpProjectedVpNormalized(3))*dMass;

velocityP.vpProjectedNotNormalizedVpx2int(a,b,c,x0,y0,z0,h,u,v) = simplify(vpProjectedVpNotNormalized(1))*dMass;
velocityP.vpProjectedNotNormalizedVpy2int(a,b,c,x0,y0,z0,h,u,v) = simplify(vpProjectedVpNotNormalized(2))*dMass;
velocityP.vpProjectedNotNormalizedVpz2int(a,b,c,x0,y0,z0,h,u,v) = simplify(vpProjectedVpNotNormalized(3))*dMass;
end

