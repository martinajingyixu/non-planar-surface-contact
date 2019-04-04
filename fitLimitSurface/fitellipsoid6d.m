function [fittingRes] = ...
    fitellipsoid6d(wrench,twist,lambda,wrenchWeight,twistWeight)

%     disp('start fitting ellipsoidal LS');
    if norm(twist) == 0
        ifFitTwist = false;
    else
        ifFitTwist = true;
    end  
    [nDim, nSamples] = size(wrench);
    scale_min = eps;
    
    % convex optimization
    cvx_begin quiet
        variable fittedCoeff(nDim,nDim) semidefinite
        variables wrenchError(nSamples) s(nSamples) twistErrorNotNormalized(nSamples)
        minimize( lambda * norm(fittedCoeff, 'fro') + wrenchWeight * sum(wrenchError) + twistWeight * sum(twistErrorNotNormalized))
%         minimize( lambda * norm(fittedCoeff, 'fro') + wrenchWeight * norm(wrenchError)...
%             + twistWeight * norm(twistErrorNotNormalized))

        subject to
            for i = 1:nSamples
               norm(wrench(:,i)' * fittedCoeff * wrench(:,i) - 1) <= wrenchError(i)
            end
            if ifFitTwist
                for i = 1:nSamples
                   norm(fittedCoeff * wrench(:,i)- s(i) * twist(:,i)) <= twistErrorNotNormalized(i)
                   s(i) >= scale_min
                end
            end
    cvx_end

    % the twistErrorNotNormalized is different than the twistPreditionError, because the
    % magnitude of the predicted twist is smaller than 1
%     predTwist = fittedCoeff*wrench;
    predTwist = fittedCoeff*wrench;

    predTwistDir = bsxfun(@rdivide, predTwist, sqrt(sum(predTwist.^2)));
%     error2 = s.*twist' - predTwist';
%     twistError2 = sum(sqrt(sum(error2.^2,2)));

    twist = bsxfun(@rdivide, twist, sqrt(sum(twist.^2)));
    twistPreditionError = sqrt(sum((predTwistDir - twist).^2));
    
%     [~,pointsSampledEllipsoid] = sample6dellipsoid(fittedCoeff,1000);
%     [hullVolume,sampledWrench]=samplefitted6dellipsoid(wrench,fittedCoeff);
%     [sampledWrench,shapeArray]=samplefittedlimitsurface(wrench,fittedCoeff,samplingStepSize,2);
    fittingRes.fittedCoeff = fittedCoeff;
%     fittedLS.shapeArray = shapeArray;
%     fittedLS.sampledWrench = sampledWrench;
    fittingRes.wrenchFitError = wrenchError;
    fittingRes.twistErrorNotNormalized = twistErrorNotNormalized;
    fittingRes.predictedTwist = predTwistDir;
    fittingRes.twistPreditionError = twistPreditionError;    
    fittingRes.s = s;

        
    
%     if (ifVisualize)
%         if nDim == 3
%             r = [fittedCoeff(1,1), fittedCoeff(2,2), fittedCoeff(3,3), ...
%                 fittedCoeff(1,2)*2, fittedCoeff(1,3)*2, fittedCoeff(2,3)*2]; % *2 because of symmetry
%             figHandle = DrawEllipsoid(r, wrench);
%             visualizewrenchtwist(wrench,twist,figHandle);
%         elseif nDim == 6
%             visualizedfitted6dgeometry(wrench,sampledWrench,twist) ;           
%         end
% 
%     end
    %% the code to visualize ellipsoid
    %% 3d
%     syms fx fy fz tx ty tz real
%     wrench = [fx;fy;fz];
%     v = sym('v', [3 3]);
%     terms3d = wrench'*v*wrench; 
%     % fx^2*v1_1 + fy^2*v2_2 + fz^2*v3_3 + fx*fy*v1_2 + fx*fy*v2_1 + fx*fz*v1_3 + fx*fz*v3_1 + fy*fz*v2_3 + fy*fz*v3_2
%     %% 6d
%     wrench = [fx;fy;fz;tx;ty;tz];
%     v = sym('v', [6 6]);
%     terms6d = wrench'*v*wrench; 
%     %fx^2*v1_1 + fy^2*v2_2 + fz^2*v3_3 + tx^2*v4_4 + ty^2*v5_5 + tz^2*v6_6 + fx*fy*v1_2 + fx*fy*v2_1 + fx*fz*v1_3 + fx*fz*v3_1 + fy*fz*v2_3 + fy*fz*v3_2 + fx*tx*v1_4 + fx*tx*v4_1 + fx*ty*v1_5 + fy*tx*v2_4 + fy*tx*v4_2 + fx*ty*v5_1 + fx*tz*v1_6 + fy*ty*v2_5 + fz*tx*v3_4 + fz*tx*v4_3 + fy*ty*v5_2 + fx*tz*v6_1 + fy*tz*v2_6 + fz*ty*v3_5 + fz*ty*v5_3 + fy*tz*v6_2 + fz*tz*v3_6 + fz*tz*v6_3 + tx*ty*v4_5 + tx*ty*v5_4 + tx*tz*v4_6 + tx*tz*v6_4 + ty*tz*v5_6 + ty*tz*v6_5

end

