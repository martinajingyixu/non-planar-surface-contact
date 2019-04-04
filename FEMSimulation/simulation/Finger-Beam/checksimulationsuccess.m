function [dbSimulationSuccess] = checksimulationsuccess(dbDir,objectParameters,objectRange)
%Check for each object, each grasp location and pressure, if the simulation
%is succesfull

    global CONSTANTS
    if ~exist('objectRange', 'var') || isempty(objectRange)
    objectRange = 1:height(objectParameters); end

%     numObject = height(objectParameters);
    numObject = numel(objectRange);
    dbSimulationSuccess = cell(numObject,1);
    
    for idxObject = 1:numel(objectRange)
        iObject = objectRange(idxObject);
        iObject
        objectInfo = table2struct(objectParameters(iObject, :));
        load([[CONSTANTS.DIR.SIMULATION '/' objectInfo.filename '/'] 'squeezeLocationArrayNoPenetration.mat'],'trueSqueezeLocationArray');
        simulationObjectRoot = [CONSTANTS.DIR.SIMULATION '/' objectInfo.filename '/'];

        dbSimulationSuccess{iObject}.xaxis = zeros(length(trueSqueezeLocationArray.xaxis),length(CONSTANTS.SIMULATION.PRESSURE_ARRAY));
        simulationObjectRootXAxis = [simulationObjectRoot 'xaxis/'];
        
        if ~objectInfo.symmetric
            dbSimulationSuccess{iObject}.yaxis = zeros(length(trueSqueezeLocationArray.yaxis),length(CONSTANTS.SIMULATION.PRESSURE_ARRAY));
            simulationObjectRootYAxis = [simulationObjectRoot 'yaxis/'];
        end
        
        for iPressure = 1:length(CONSTANTS.SIMULATION.PRESSURE_ARRAY)
            
            for iSqueezeLocation = 1:length(trueSqueezeLocationArray.xaxis)
                simulationResultPath = simulationresultpath(simulationObjectRootXAxis,...
            iSqueezeLocation,trueSqueezeLocationArray.xaxis,CONSTANTS.SIMULATION.PRESSURE_ARRAY(iPressure));
                ifSuccuess = issimulationsuccess(simulationResultPath);
                dbSimulationSuccess{iObject}.xaxis(iSqueezeLocation,iPressure) = ifSuccuess;

            end
            
            if ~objectInfo.symmetric
                for iSqueezeLocation = 1:length(trueSqueezeLocationArray.yaxis)
                    simulationResultPath = simulationresultpath(simulationObjectRootYAxis,...
                iSqueezeLocation,trueSqueezeLocationArray.yaxis,CONSTANTS.SIMULATION.PRESSURE_ARRAY(iPressure));
                    ifSuccuess = issimulationsuccess(simulationResultPath);

                    dbSimulationSuccess{iObject}.yaxis(iSqueezeLocation,iPressure) = ifSuccuess;

                end
            end
            
        end

    save([dbDir 'dbSimulationSuccess.mat'],'dbSimulationSuccess');
        
    end   
end

function [ pathSqueezeLocation ]= simulationresultpath(root,iSqueezeLocation,squeezeLocationArray,pressure)

    pathSqueezeLocation = [root '/pressure' num2str(pressure) '/' 'loc' num2str(iSqueezeLocation) '_x' num2str(squeezeLocationArray(1,iSqueezeLocation)) '_y' num2str(squeezeLocationArray(2,iSqueezeLocation)) '_z' num2str(squeezeLocationArray(3,iSqueezeLocation))];
end

