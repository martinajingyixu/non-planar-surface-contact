function rangeU = computerange(radius1, radius2,  halfLongAxis, shortAxis)

% compute range U from measured distances of the contact 

%% for elliptical cylinder

% http://www.petercollingridge.co.uk/tutorials/computational-geometry/finding-angle-around-ellipse/
% a = radius1;
% b = radius2;
% py = halfLongAxis;
% px = b - shortAxis;
% theta1 = atan(b / a * py / px);    

theta1 = atan(radius1 / radius2 * (radius2-shortAxis) / halfLongAxis);

% rangeU = [pi/2 - theta1, pi/2+theta1];
rangeU = [ theta1, pi-theta1];



end

