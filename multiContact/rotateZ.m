function [ Rz ] = rotateZ( anglez)
%get rotation matrix given anglex, angley and anglez
%   Notice: angles are counter-clock wise!!
%https://en.wikipedia.org/wiki/Rotation_matrix

Rz = [cos(anglez) -sin(anglez) 0; sin(anglez) cos(anglez) 0;0,0,1];
end

