function [ rotation_matrix ] = getRotationMatrix( anglex, angley, anglez)
%get rotation matrix given anglex, angley and anglez
%   Notice: angles are counter-clock wise!!
%https://en.wikipedia.org/wiki/Rotation_matrix

Rx = [1,0,0;0 cos(anglex), -sin(anglex),;0, sin(anglex), cos(anglex)];
Ry = [cos(angley) 0 sin(angley); 0 1 0; -sin(angley) 0 cos(angley)];
Rz = [cos(anglez) -sin(anglez) 0; sin(anglez) cos(anglez) 0;0,0,1];
rotation_matrix = Rz*Ry*Rx;
% rotation_matrix = (coord_old)^(-1) * coord_new;
end

