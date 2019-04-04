function [meanWrenchError,meanTwistError] = computefittingerrorellip(fittedCoeff,wrenchSamples,twistSamples)
    
    wrenchError = zeros(length(wrenchSamples),1);
    for i=1:length(wrenchSamples)
        wrenchError(i) =  norm(wrenchSamples(:,i)' * fittedCoeff * wrenchSamples(:,i) - 1);
    end
    meanWrenchError = mean(wrenchError);
    predictedTwist= fittedCoeff * wrenchSamples;
    predictedTwist = predictedTwist./vecnorm(predictedTwist);
    twist = twistSamples./vecnorm(twistSamples);  
    meanTwistError = mean(acosd(dot(predictedTwist,twist)));

end

