function [uniqueIdx] = getidxconvexhull(point3d)
%return unique index of the points which are the convex hull of the input vector
if size(point3d,1) == 3
    x = point3d(1,:)';
    y = point3d(2,:)';
    z = point3d(3,:)';
elseif size(point3d,2) == 3
    x = point3d(:,1);
    y = point3d(:,2);
    z = point3d(:,3);
end
    uniqueIdx = unique(convhull(x,y,z));

end

