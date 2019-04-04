function [CORpoint,distOrigin2COR] = computePointCORGoesThrough(origin,a,b,c,dist1,dist2)
% Compute the point that COR goes through, given:
%1. the orientation of the COR, defined by [a;b;c]
%2. the vector from the origin to the COR (dist) should be perpendicular to COR 
%3. we sample the distance from origin to COR, based on dist1, dist2,
%compute dist3 since the dot product is 0.
% sampling dist1 and dist2 to get distX, distY, and distZ.


% divide to three cases. Compare value of a,b,c to avoid divide by 0.
%% this part is suitable is max(abs(a),abs(b),abs(c)) = a; 
% Means the COR is near x axis.
distYLA = dist1;
distZLA = dist2;
distXLA = -(distYLA*b+distZLA*c)/a;
PintLA = simplify(origin + [distXLA;distYLA;distZLA]);  


CORpoint.LargeA.x0 = PintLA(1);
CORpoint.LargeA.y0 = PintLA(2);
CORpoint.LargeA.z0 = PintLA(3);

distOrigin2COR.LargeA.distX = distXLA;
distOrigin2COR.LargeA.distY = distYLA;
distOrigin2COR.LargeA.distZ = distZLA;

%% this part is suitable is max(abs(a),abs(b),abs(c)) = b; 
% Means the COR is near y axis.
distXLB = dist1;
distZLB = dist2;
distYLB = -(distXLB*a+distZLB*c)/b;

PintLB = simplify(origin + [distXLB;distYLB;distZLB]);  

CORpoint.LargeB.x0 = PintLB(1);
CORpoint.LargeB.y0 = PintLB(2);
CORpoint.LargeB.z0 = PintLB(3);

distOrigin2COR.LargeB.distX = distXLB;
distOrigin2COR.LargeB.distY = distYLB;
distOrigin2COR.LargeB.distZ = distZLB;

%% this part is suitable is max(abs(a),abs(b),abs(c)) = c; 
% Means the COR is near z axis.
distXLC = dist1;
distYLC = dist2;
distZLC = -(distXLC*a+distYLC*b)/c;

PintLC = simplify(origin + [distXLC;distYLC;distZLC]);  

CORpoint.LargeC.x0 = PintLC(1);
CORpoint.LargeC.y0 = PintLC(2);
CORpoint.LargeC.z0 = PintLC(3);


distOrigin2COR.LargeC.distX = distXLC;
distOrigin2COR.LargeC.distY = distYLC;
distOrigin2COR.LargeC.distZ = distZLC;


end

