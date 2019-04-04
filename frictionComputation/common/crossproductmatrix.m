function [a_hat] = crossproductmatrix(a)

a_hat = [0,-a(3),a(2);a(3),0,-a(1);-a(2),a(1),0];

end

