function [ Ry ] = rotateY( angley)
%get rotation matrix given anglex, angley and anglez
%   Notice: angles are counter-clock wise!!
%https://en.wikipedia.org/wiki/Rotation_matrix
Ry = [cos(angley) 0 sin(angley); 0 1 0; -sin(angley) 0 cos(angley)];
end

