function [roundedA] = changepercision(A,numDigits)
%Change the percision of the matrix A to numDigits digits
if nargin == 1
    numDigits = 5;
end
roundedA = A;
roundedA(A<10^-numDigits) = 0;
roundedA = double(vpa(roundedA,numDigits));

end

