function [rSampledCatesianCoord] = parametric2cartesiansamples...
    (r,sampledVar1,sampledVar2)
%Visualize the surface with the parametric form

    rSampledCatesianCoord = zeros(3,numel(sampledVar1)*numel(sampledVar2));
    iSample = 1;
    for curSampleVar1 = sampledVar1
        for curSampleVar2 = sampledVar2
            rSampledCatesianCoord(:,iSample) = ...
                r(curSampleVar1,curSampleVar2);
            iSample = iSample + 1;
        end
    end
   
end

