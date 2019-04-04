function [O,surfaceArea,forceSum] = frictioncenter(r,pressureDist,dA,var1,var2,rangeVar1,rangeVar2)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

lowerBoundVar1 = rangeVar1(1);
upperBoundVar1 = rangeVar1(2);

lowerBoundVar2 = rangeVar2(1);
upperBoundVar2 = rangeVar2(2);

surfaceArea = int(int(dA, var1,lowerBoundVar1,upperBoundVar1)...
    ,var2,lowerBoundVar2,upperBoundVar2);

forceSum = int(int(dA*pressureDist, var1,lowerBoundVar1,upperBoundVar1)...
    ,var2,lowerBoundVar2,upperBoundVar2);

Ox = int(int(r(1)*dA*pressureDist, var1,lowerBoundVar1,upperBoundVar1)...
    ,var2,lowerBoundVar2,upperBoundVar2)./forceSum;


Oy = int(int(r(2)*dA*pressureDist, var1,lowerBoundVar1,upperBoundVar1)...
    ,var2,lowerBoundVar2,upperBoundVar2)./forceSum;

Oz = int(int(r(3)*dA*pressureDist, var1,lowerBoundVar1,upperBoundVar1)...
    ,var2,lowerBoundVar2,upperBoundVar2)./forceSum;
O = vpa([Ox;Oy;Oz]);
% O = changepercision(vpa([Ox;Oy;Oz]));
% surfaceArea = changepercision(vpa(surfaceArea));
% forceSum = changepercision(vpa(forceSum));

end

