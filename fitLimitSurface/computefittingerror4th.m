function [meanWrenchError,meanTwistError] = computefittingerror4th(fittedCoeff,wrenchSamples,twistSamples)

    x1 = wrenchSamples(1,:);
    x2 = wrenchSamples(2,:);
    x3 = wrenchSamples(3,:);
    x4 = wrenchSamples(4,:);
    x5 = wrenchSamples(5,:);
    x6 = wrenchSamples(6,:);

    monomials4thOrder  = monomials4thorderpoly(x1,x2,x3,x4,x5,x6);
    gradient = gradient4thorderpoly(x1,x2,x3,x4,x5,x6);
    
    meanWrenchError = mean(sqrt((monomials4thOrder * fittedCoeff -1).^2));
    for i = 1:length(twistSamples)
        gradientOneSample(:,:) = gradient(i,:,:);
        predictedTwist(:,i) = gradientOneSample' * fittedCoeff;
    end
    
     predictedTwist = predictedTwist./vecnorm(predictedTwist);
    twist = twistSamples./vecnorm(twistSamples);
    meanTwistError = mean(acosd(dot(predictedTwist,twist)));

end

