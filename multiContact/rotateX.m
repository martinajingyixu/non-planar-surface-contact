function [ Rx ] = rotateX( anglex)
%get rotation matrix given anglex, angley and anglez
%   Notice: angles are counter-clock wise!!
%https://en.wikipedia.org/wiki/Rotation_matrix

Rx = [1,0,0;0 cos(anglex), -sin(anglex),;0, sin(anglex), cos(anglex)];
end

