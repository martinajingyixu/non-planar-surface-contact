function [twist] = computetwist(origin,twistSamples)
% origin can be friction center or center or mass (COM)
% compute twist from twist samples
omega = twistSamples(1:3,:);
dist = twistSamples(4:6,:);
pitch = twistSamples(7,:);

rn = repmat(origin,[1,length(dist)])-dist;

v = pitch .* omega + cross(omega,rn);


twist = [v;omega];
end



