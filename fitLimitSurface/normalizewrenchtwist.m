function [wrenchNormalized,twistNormalized] = normalizewrenchtwist(wrench,twist)
%     wrenchNormalized = wrench;
%     twistNormalized = twist;
%     wrenchNormalized(4:6,:) = wrench(4:6,:)./pho;
%     twistNormalized(4:6,:) = twist(4:6,:).* pho;
% 
    wrenchMaxOri = max(abs(wrench)');
    wrenchMax = wrenchMaxOri;
    wrenchMax(find(wrenchMax<1e-10)) = Inf;
    wrenchNormalized = (wrench'./wrenchMax)';
    twistNormalized = (twist'.* wrenchMaxOri)';

%     wrenchMax = max(abs(wrench)');
%     wrenchNormalized = (wrench'./wrenchMax)';
%     twistNormalized = (twist'.* wrenchMax)';
    twistNormalized = bsxfun(@rdivide, twistNormalized, sqrt(sum(twistNormalized.^2)));

end
