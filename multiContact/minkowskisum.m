function [sum] = minkowskisum(setA,setB)
%every item in setA plus every item in setB and add it if unique
sum = zeros(6,length(setA) * length(setB));
iTotal = 1;
for iA = 1:length(setA)
    for iB = 1:length(setB)
        sum(:,iTotal) = setA(:,iA) + setB(:,iB);
        iTotal = iTotal+1;
    end
end
sum = unique(sum','rows')';
end

