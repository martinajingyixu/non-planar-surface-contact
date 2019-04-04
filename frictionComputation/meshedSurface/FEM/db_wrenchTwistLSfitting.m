function [contactInfo,fittingResults,wrenchTwist,objFittingResults] = ...
    db_wrenchTwistLSfitting(pathResults,FEMResults,isSqueezingX,config)
    [numLocation,numPressure] = size(FEMResults);
    elementsLeft = cell(numLocation,numPressure);
    elementsRight = cell(numLocation,numPressure);   
    fittingResultsLeft = cell(numLocation,numPressure);
    fittingResultsRight= cell(numLocation,numPressure);
    wrenchTwistLeft = cell(numLocation,numPressure);
    wrenchTwistRight = cell(numLocation,numPressure);
    objFittingResultsContact1 = NaN(numLocation,numPressure,4);
    objFittingResultsContact2 = NaN(numLocation,numPressure,4);
   
    for iPressure = 1:numPressure
        disp(['iPressure ' num2str(iPressure)]);
       parfor iLocation = 1:numLocation
            try                 
                contact = getFEMresults(FEMResults{iLocation,iPressure},isSqueezingX,config.ifVisualize);                
                hold on
                title(num2str(iLocation))
                [fitRes1,wrenchTwist1] = ...
                    computewrenchtwistfitls(contact.left,pathResults,config);
                
                [fitRes2,wrenchTwist2] = ...
                    computewrenchtwistfitlsdiscretizedsurfaces(contact.right,pathResults,config);   
                objFittingResultsContact1(iLocation,iPressure,:) = ...
                    [fitRes1.fittingErr.meanWrenchError4th,fitRes1.fittingErr.meanWrenchErrorEllip,...
                    fitRes1.fittingErr.meanTwistError4th,fitRes1.fittingErr.meanTwistErrorEllip];
                
                objFittingResultsContact2(iLocation,iPressure,:) = ...
                    [fitRes2.fittingErr.meanWrenchError4th,fitRes2.fittingErr.meanWrenchErrorEllip,...
                    fitRes2.fittingErr.meanTwistError4th,fitRes2.fittingErr.meanTwistErrorEllip];                

            if config.ifVerbose == 1
                displayfittingerror(fitRes1.fittingErr);
            end
            
            catch
                contact.left = [];
                contact.right = [];
                fitRes1 = [];
                fitRes2 = [];
                wrenchTwist1 = [];
                wrenchTwist2 = [];                
            end

            elementsLeft{iLocation,iPressure} = contact.left;
            elementsRight{iLocation,iPressure} = contact.right;
            fittingResultsLeft{iLocation,iPressure} = fitRes1;
            fittingResultsRight{iLocation,iPressure} = fitRes2;
            wrenchTwistLeft{iLocation,iPressure} = wrenchTwist1;
            wrenchTwistRight{iLocation,iPressure} = wrenchTwist2;
        end
    end
        contactInfo.elementsLeft = elementsLeft;
        contactInfo.elementsRight = elementsRight;
        fittingResults.left = fittingResultsLeft;
        fittingResults.right = fittingResultsRight;
        wrenchTwist.left = wrenchTwistLeft;
        wrenchTwist.right = wrenchTwistRight;
        objFittingResults.contact1 = reshape(objFittingResultsContact1,[numLocation*numPressure,4]) ;
        objFittingResults.contact2 = reshape(objFittingResultsContact2,[numLocation*numPressure,4]) ;
    
end

