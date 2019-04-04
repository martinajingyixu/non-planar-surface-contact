function [ frictionalForce,frictionalTorque,vnp,vn] ...
    = frictiontriangularelement(elementInfo, pwc, CORAxis, pitch, mu)

%Compute max torque and max force of one element. 
%   !!Important: the vertics are all in world coordinate system!! Need to
%   transform them first to new local coordinate system with the Transformation matrix : TransMat!!
%   COR = [rx,ry,0]. Its the rotation axis parallel to z axis. Its in local
%   Coord system already

% pitch is the magnitude of the velocity, parallel to the COR/w.  

%% vertics of one element
% plane 1
% vertics of one element

elementCenter = elementInfo.center;
elementNormal = elementInfo.normal;
pressure = elementInfo.pressure;
area = elementInfo.area;

%% arm calculation: distance between element center to the contact normal, which goes through the projectedPWC
% torque arm should be the "Ortsvector", im Bezug auf O.
torqueArm =  elementCenter - pwc;
%% compute hand written equations
forceMagnitude = mu*pressure*area;
%% Compute COR related vectors to compute the direction of the friction force.
% 2.3.1: dn
% Pd0 = CORAxis.location + 1000 * CORAxis.direction;
% Pd1 = CORAxis.location - 1000 * CORAxis.direction;
% [intersectionPointCOR,~]=planelineintersect(CORAxis.direction,elementCenter,Pd0,Pd1);
% 
% dn = elementCenter - intersectionPointCOR;
% dnNorm = dn./norm(dn);
dn = elementCenter - CORAxis.location;
% rn = CORAxis.location - elementCenter;

CORAxisDirection = CORAxis.direction./norm(CORAxis.direction);
% 2.3.2
% vn = cross(CORAxisDirection,dnNorm); % no parallel velocity
% vn = vn./norm(vn);
% vn = pitch .* CORAxisDirection + cross(dn,CORAxisDirection);
vn = pitch .* CORAxisDirection + cross(CORAxisDirection,dn);
vn = vn./norm(vn);
% 2.3.3
vnp = vn - dot(vn,elementNormal) /(norm(elementNormal))^2 *elementNormal;
% forceDirction = -vnp./norm(vnp); % vnp is projected vn onto the element plane
forceDirction = vnp./norm(vnp); % vnp is projected vn onto the element plane

frictionalForce = forceDirction * forceMagnitude;
frictionalTorque = cross(torqueArm,frictionalForce);

end

