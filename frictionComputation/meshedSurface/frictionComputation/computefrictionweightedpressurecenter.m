function [pwc] = computefrictionweightedpressurecenter(centerArray,pressureArray,areaArray)
%computefrictionweightedpressurecenter Computes the pressure weighted
%friction center
% May not be on the contact area

    sumForce = sum(pressureArray.*areaArray);
     
    pwcx = sum(centerArray(1,:).*...
        pressureArray.*areaArray)./sumForce;
    pwcy = sum(centerArray(2,:).*...
        pressureArray.*areaArray)./sumForce;
    pwcz = sum(centerArray(3,:).*...
        pressureArray.*areaArray)./sumForce;

    pwc = [pwcx;pwcy;pwcz];
end

