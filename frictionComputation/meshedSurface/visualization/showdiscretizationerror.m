function [errorArray] = showdiscretizationerror(wrenchErrorArray,surfaceTypeArray)

    numSurfaceType = numel(wrenchErrorArray);
    errorArray = cell(numSurfaceType,1);
    for iSurface = 1:numSurfaceType
        wrenchError = wrenchErrorArray{iSurface};
        wrenchDim = size(wrenchError{1},1);
        error = zeros(numel(wrenchError),wrenchDim);
        for iElementSize = 1:numel(wrenchError)
            error(iElementSize,:) = mean(wrenchError{iElementSize},2);    
        end
        errorArray{iSurface,1} = error;
    end
    for iWrench2Plot = 1:wrenchDim
        
        figure
        for iSurface = 1:numSurfaceType
%           plot(1:numel(wrenchError),errorArray{iSurface}(:,iWrench2Plot),...
%                 'DisplayName',surfaceTypeArray{iSurface});
            plot(1:numel(wrenchError),errorArray{iSurface}(:,iWrench2Plot));
            hold on
        end
        legend(surfaceTypeArray)

        hold off
    end

end

