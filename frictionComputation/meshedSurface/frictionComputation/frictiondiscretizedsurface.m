function [forceArray,torqueArray] = frictiondiscretizedsurface(TwistSample,elementInfoArray,pwc,numElements,mu)
    CORAxis.direction = TwistSample(1:3);
    CORAxis.location = TwistSample(4:6);
    pitch = TwistSample(end);
    [forceArray,torqueArray] = deal(zeros(3,numElements));

    for iElement = 1:numElements
        [ frictionalForce,frictionalTorque,projectedVel,Vel] ...
        = frictiontriangularelement(elementInfoArray{iElement}, pwc, CORAxis , pitch,mu);
        forceArray(:,iElement) = frictionalForce;
        torqueArray(:,iElement) = frictionalTorque;

    end
end

