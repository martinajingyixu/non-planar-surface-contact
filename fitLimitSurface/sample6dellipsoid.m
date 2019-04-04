function [pointsRotBeforeRejection,pointsRot] = ...
    sample6dellipsoid(A,numPoints,iffilterpoints,ifshowfig)
% https://en.wikipedia.org/wiki/Ellipsoid#Generalised_equations
% The general form of an arbitrarily oriented ellipsoid, cnetered at v:
% (x-v)' * A * (x-v) = 1
% Here assume v = [0;0;0];
% A is a 3x3 positive definite matrix. The eigenvalues of A are the
% reciprocals of the squares of the semi-axis. The eigenvectors are the
% directions of the axes of the ellipsoid.
if nargin<3

    iffilterpoints = true; 
end

if nargin<4
    ifshowfig = false;
end

[axesDirections,radii2] = eig(A);
radii = sqrt(1./eig(A));
a1 = sqrt(1/radii2(1,1)); 
a2 = sqrt(1/radii2(2,2));
a3 = sqrt(1/radii2(3,3));
a4 = sqrt(1/radii2(4,4));
a5 = sqrt(1/radii2(5,5));
a6 = sqrt(1/radii2(6,6));

rotationMatrix = axesDirections;

%% sample the ellipsoid
% 1. samples the points with normal distrubtion, use the respective radii
% as the variances
% 2. scale them using the distance to the center (d)
% https://math.stackexchange.com/questions/973101/how-to-generate-points-uniformly-distributed-on-the-surface-of-an-ellipsoid
% fxx = [-a1,unifrnd(-a1,a1,[1,numPoints]),a1];
% fyy = [-a2,unifrnd(-a2,a2,[1,numPoints]),a2];
% fzz = [-a3,unifrnd(-a3,a3,[1,numPoints]),a3];
% txx = [-a4,unifrnd(-a4,a5,[1,numPoints]),a4];
% tyy = [-a5,unifrnd(-a5,a5,[1,numPoints]),a5];
% tzz = [-a6,unifrnd(-a6,a6,[1,numPoints]),a6];
% fxx = [unifrnd(-a1^2,a1^2,[1,numPoints])];
% fyy = [unifrnd(-a2^2,a2^2,[1,numPoints])];
% fzz = [unifrnd(-a3^2,a3^2,[1,numPoints])];
% txx = [unifrnd(-a4^2,a5^2,[1,numPoints])];
% tyy = [unifrnd(-a5^2,a5^2,[1,numPoints])];
% tzz = [unifrnd(-a6^2,a6^2,[1,numPoints])];
% 
a1S = 1.2*a1;if radii2(1,1)<1e-4 a1S = 0;end % if relative radii is 0, then not sampling this dimension
a2S = 1.2*a2;if radii2(2,2)<1e-4 a2S = 0;end
a3S = 1.2*a3;if radii2(3,3)<1e-4 a3S = 0;end
a4S = 1.2*a4;if radii2(4,4)<1e-4 a4S = 0;end
a5S = 1.2*a5;if radii2(5,5)<1e-4 a5S = 0;end
a6S = 1.2*a6;if radii2(6,6)<1e-4 a6S = 0;end


fxx = unique([-a1S:a1S/numPoints*2:a1S,0]);
fyy = unique([-a2S:a2S/numPoints*2:a2S,0]);
fzz = unique([-a3S:a3S/numPoints*2:a3S,0]);
txx = unique([-a4S:a4S/numPoints*2:a4S,0]);
tyy = unique([-a5S:a5S/numPoints*2:a5S,0]);
tzz = unique([-a6S:a6S/numPoints*2:a6S,0]);


combi = combvec(fxx,fyy,fzz,txx,tyy,tzz);
fxx = combi(1,:);
fyy = combi(2,:);
fzz = combi(3,:);
txx = combi(4,:);
tyy = combi(5,:);
tzz = combi(6,:);

lambda = sqrt(fxx.^2./a1^2+fyy.^2./a2^2+fzz.^2./a3^2 +...
    txx.^2./a4^2+tyy.^2./a5^2+tzz.^2./a6^2);

points = [fxx./lambda;fyy./lambda;fzz./lambda;txx./lambda;tyy./lambda;tzz./lambda];
pointsRot = rotationMatrix * points;
pointsRotBeforeRejection = pointsRot;
% if ifshowfig
%     plot6dwrenches(pointsRot')
% end

%% remove the points that are too close
% get the distance of each two points. 
if iffilterpoints
    pairwisedist = pdist(pointsRot','minkowski');
    % % convert the idx of distance to a pair of two points with pair_p2
    % % than pair_p1. Then remove the pair_p2
    distNum= numPoints-1:-1:1;
    pintsIdx = 1:1:numPoints-1;
    % pair_p1 = [];
    pair_p2 = [];
    for iPoints = 1:numPoints-1
    %     pair_p1 = [pair_p1,repmat(pintsIdx(iPoints),[1,distNum(iPoints)])];
        pair_p2 = [pair_p2,pintsIdx(iPoints)+1:pintsIdx(iPoints)+distNum(iPoints)];

    end
    % approximate the mean distance of two points by computing the surface area
    % then compute the number of rectangles. 

    % approximate surface area with 6d sphere
%     a = mean(radii);
%     surfaceArea =  a^5 * pi^3 *2;
%     meanDist = sqrt(surfaceArea./numPoints);
    meanDist = min(radii)/numPoints*10;
    
    %  pairs = [pair_p1(idxSmallPoints);pair_p2(idxSmallPoints)]';
    % remove the points which has small distance
    idxDel = unique(pair_p2(find(pairwisedist<meanDist)));
    pointsRot(:,idxDel) = [];
end
if ifshowfig
    plot6dwrenches(pointsRot')
end

end

