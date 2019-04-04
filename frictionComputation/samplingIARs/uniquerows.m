function [B,idx] = uniquerows(A)
% get unique rows. Sometimes it has +0.00 and -0.00
% put + - 0.000 to 0
AProcessZeros = not(abs(A)<1e-10) .*A;
[B,idx] = uniquetol(AProcessZeros,'ByRows',0.00001);

end

