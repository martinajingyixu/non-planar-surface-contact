function [] = db_simulation(dbDir,objectParameters,objectRange)
global CONSTANTS
%db_simulation FEM simulation for db
% create directory for simulation dir if not exist
    createdir(CONSTANTS.DIR.SIMULATION);
    if ~exist('objectRange', 'var') || isempty(objectRange)
        objectRange = 1:height(objectParameters); end
    
    dbSimulationError = cell(size(objectRange));
%     for iObject = objectRange
    for idxObject = 1:numel(objectRange)
        iObject = objectRange(idxObject);
        % simulation for each object
        objectInfo = table2struct(objectParameters(iObject, :));
        fprintf('Model %5d/%5d: %s\n', iObject, height(objectParameters), objectInfo.filename);
        squeezePressureArray = CONSTANTS.SIMULATION.PRESSURE_ARRAY;
        savemat([CONSTANTS.DIR.SIMULATION '/' objectInfo.filename '/'], 'squeezePressureArray.mat',{'squeezePressureArray'},{squeezePressureArray});
        simulationError = FEMsimulation(dbDir,objectInfo);
        dbSimulationError{iObject}.xaxis = simulationError.xaxis;
        dbSimulationError{iObject}.yaxis = simulationError.yaxis;

        savemat([CONSTANTS.DIR.SIMULATION '/' objectInfo.filename '/'], 'dbSimulationError.mat',{'dbSimulationError'},{dbSimulationError});

    end
end

function [simulationSuccess] = FEMsimulation(dbDir, objectInfo)
global CONSTANTS

% load squeeze locations
   load([dbDir objectInfo.filename '.mat'],'squeezeLocationArray');
   load([dbDir objectInfo.filename '.mat'],'sall');
%    load([dbDir objectInfo.filename '.mat'],'rshift'); % for cien, assymetric objects, the shift of center on x axis
   load([dbDir objectInfo.filename '.mat'],'sshift'); % for cien, assymetric objects, the shift of center on x axis

   % sall represents the object geometry
   [squeezeLocationArrayNoPenetration.xaxis,squeezeLocationArrayNoPenetration.xaxisNegative] = contactlocationpreventpenetration...
   (objectInfo.r1, objectInfo.h, squeezeLocationArray.xaxis(3,:), sall,sshift);

   [squeezeLocationArrayNoPenetration.yaxis,squeezeLocationArrayNoPenetration.yaxisNegative]  = contactlocationpreventpenetration...
   (objectInfo.r2, objectInfo.h, squeezeLocationArray.yaxis(3,:), sall,0);
   
    trueSqueezeLocationArray.xaxis = [squeezeLocationArrayNoPenetration.xaxis';squeezeLocationArray.xaxis(2:3,:);];
    trueSqueezeLocationArray.yaxis = [squeezeLocationArray.yaxis(1,:);squeezeLocationArrayNoPenetration.yaxis';squeezeLocationArray.yaxis(3,:);];
    trueSqueezeLocationArray.xaxisNegative = [squeezeLocationArrayNoPenetration.xaxisNegative';squeezeLocationArray.xaxis(2:3,:);];
    trueSqueezeLocationArray.yaxisNegative = [squeezeLocationArray.yaxis(1,:);squeezeLocationArrayNoPenetration.yaxisNegative';squeezeLocationArray.yaxis(3,:);];

    savemat([CONSTANTS.DIR.SIMULATION '/' objectInfo.filename '/'], 'squeezeLocationArrayNoPenetration.mat',{'squeezeLocationArrayNoPenetration'},{squeezeLocationArrayNoPenetration});
    savemat([CONSTANTS.DIR.SIMULATION '/' objectInfo.filename '/'], 'squeezeLocationArrayNoPenetration.mat',{'trueSqueezeLocationArray'},{trueSqueezeLocationArray});

    simulationObjectRoot = [CONSTANTS.DIR.SIMULATION '/' objectInfo.filename '/'];

    if exist([simulationObjectRoot 'simulationSuccess.mat'],'file')
        load([simulationObjectRoot 'simulationSuccess.mat'],'simulationSuccess');
    else
        simulationSuccess.xaxis = zeros(length(squeezeLocationArray.xaxis),length(CONSTANTS.SIMULATION.PRESSURE_ARRAY));
        simulationSuccess.yaxis = zeros(length(squeezeLocationArray.yaxis),length(CONSTANTS.SIMULATION.PRESSURE_ARRAY));

    end


    % squeeze location is later used to specity the lowest point of the 
    % vh sensor, not the center
   
    % objectBounday only used to later select the obect for meshing
%     objectBoundayX = objectInfo.r1; 
%     objectBoundayY = objectInfo.r2;
    
    
    simulationObjectRootXAxis = [simulationObjectRoot 'xaxis/'];
    simulationObjectRootYAxis = [simulationObjectRoot 'yaxis/'];

    for iPressure = 1:length(CONSTANTS.SIMULATION.PRESSURE_ARRAY)
%     for iPressure = 5

        squeezePressure = CONSTANTS.SIMULATION.PRESSURE_ARRAY(iPressure);
%%         
        % run simulation with several squeeze locations along x axis
        isSuccessfulArray = runsimulation(dbDir,simulationObjectRootXAxis,objectInfo, ...
            CONSTANTS.MODEL.XAXIS,squeezePressure, ...
            trueSqueezeLocationArray.xaxis,trueSqueezeLocationArray.xaxisNegative,CONSTANTS.SIMULATION.SQUEEZE_XAXIS);
        simulationSuccess.xaxis(:,iPressure) = isSuccessfulArray;
        save([simulationObjectRoot 'simulationSuccess.mat'],'simulationSuccess');
%%        
        % run simulation with several squeeze locations along y axis
        if objectInfo.symmetric == false
            isSuccessfulArray = runsimulation(dbDir,simulationObjectRootYAxis,objectInfo, ...
                 CONSTANTS.MODEL.YAXIS,squeezePressure, ...
                trueSqueezeLocationArray.yaxis,trueSqueezeLocationArray.yaxisNegative,CONSTANTS.SIMULATION.SQUEEZE_YAXIS); % no center shift on y axis
            simulationSuccess.yaxis(:,iPressure) = isSuccessfulArray;
            save([simulationObjectRoot 'simulationSuccess.mat'],'simulationSuccess');
        end        
    end  

end

function [isSuccessfulArray] = runsimulation(dbDir,simulationResultsRoot,objectInfo, ...
    SIMULATION_PARAMETERS,squeezePressure, simulationLocationArray,simulationLocationArrayNegative,shouldRunSimulation)
    isSuccessfulArray = zeros(length(simulationLocationArray),1);
    dirIGESModels = [dbDir '/iges_models/'];
    
    for iSqueezeLocation = 1:length(simulationLocationArray)-2 % for the last 2 squeeze location, the vh sensor and the object does not have contact at all
%     for iSqueezeLocation = 1:3 % for the last 2 squeeze location, the vh sensor and the object does not have contact at all

        squeezeLocationPos = simulationLocationArray(:,iSqueezeLocation);
        squeezeLocationNeg = simulationLocationArrayNegative(:,iSqueezeLocation);

        simulationResultPath = simulationresultpath(simulationResultsRoot,...
            iSqueezeLocation,simulationLocationArray,squeezePressure);
        createdir(simulationResultPath);

        isSimulationed = issimulated(simulationResultPath);
        if ~isSimulationed && shouldRunSimulation
            
            isSuccessfulArray(iSqueezeLocation,1) = ...
                runansysfingerellipsoid(dirIGESModels,simulationResultPath,objectInfo,...
                SIMULATION_PARAMETERS,squeezeLocationPos,squeezeLocationNeg,squeezePressure);
            
            disp(['If Success:' num2str(isSuccessfulArray(iSqueezeLocation,1)) ' ' simulationResultPath ...
                ]);
        end

    end
        
end

function [squeezeLocationArrayNoPenetration,squeezeLocationArrayNoPenetrationNeg] = contactlocationpreventpenetration...
    (objectRadiusOneAxis, objectHeight, graspLocationHeights, ppObjectGeometry,centerShift)
% compute the contact location which deals with the penetration problem.
% Compute the object widest part within the sensor range of height. 
% This location is the sensor distance to object center axis. Because then
% there will be contact with center and object and no penetration.

    global CONSTANTS
    squeezeLocationArrayNoPenetration = zeros(length(graspLocationHeights),1);
    squeezeLocationArrayNoPenetrationNeg = zeros(length(graspLocationHeights),1);

    squeezeLocationMaxHeight = max(graspLocationHeights);
    stepSz = 0.0001;
     
    for iSqueezeLocation = 1:length(graspLocationHeights)
%         squeezeLocationRawDist = squeezeLocationDistToCenterAxis(iSqueezeLocation);
        squeezeLocationHeight = graspLocationHeights(iSqueezeLocation);
        heightsInRange = squeezeLocationHeight:stepSz:...
            min(squeezeLocationMaxHeight,squeezeLocationHeight+ 2*CONSTANTS.SIMULATION.SENSOR_HEIGHT);% the height is radius of the finger
%         heightsInRange = squeezeLocationHeight;
        
        try
            shiftedCenterArray = objectRadiusOneAxis * ppval(centerShift, heightsInRange/objectHeight);
        catch
            shiftedCenterArray = 0;
        end
        heightsInRangeFinger  = stepSz:stepSz:CONSTANTS.SIMULATION.SENSOR_HEIGHT;
%         heightArray = heightsInRange - squeezeLocationHeight - CONSTANTS.SIMULATION.SENSOR_HEIGHT;
%         heightsInRangeFinger = heightArray(heightArray>0);
        % (x/a)^2 + (z/c)^2 = 1; know a,c,z (z is the height), compute
        % possible x
        halfLength = sqrt((1- ((heightsInRangeFinger)./CONSTANTS.SIMULATION.SENSOR_HEIGHT).^2) ...
            .*(CONSTANTS.SIMULATION.FINGER_THICK^2)); % use ellipse function to compute the possible coordinate within the range
%         fingerLocationArray = [fliplr(halfLength),0,halfLength];
        fingerLocationArray = -[fliplr(halfLength),0,(halfLength)];
        fingerLocationArray = fingerLocationArray(1:numel(heightsInRange));
            
        
        % no penetration check for the finger because of the shape
        distancesInRangePos = objectRadiusOneAxis * ppval(ppObjectGeometry,heightsInRange/objectHeight)+shiftedCenterArray;
        distancesInRangeNeg = -objectRadiusOneAxis * ppval(ppObjectGeometry,heightsInRange/objectHeight)+shiftedCenterArray;
        fingerLocationArrayPos = fingerLocationArray + distancesInRangePos(1);
        fingerLocationArrayNeg = -fingerLocationArray + distancesInRangeNeg(1);
%         fingerLocationArrayPos = heightsInRangeFinger + distancesInRangePos(1);
%         fingerLocationArrayNeg = -heightsInRangeFinger+ distancesInRangeNeg(1); 
%         
        squeezeLocationArrayNoPenetration(iSqueezeLocation) = distancesInRangePos(1) - ...
            min(fingerLocationArrayPos - distancesInRangePos) - CONSTANTS.SIMULATION.FINGER_THICK + 0.5; % prevent penetration
        squeezeLocationArrayNoPenetrationNeg(iSqueezeLocation) = distancesInRangeNeg(1) -...
            max(fingerLocationArrayNeg -distancesInRangeNeg) + CONSTANTS.SIMULATION.FINGER_THICK + 0.5;

        %%
%         figure
%         objh = 0:5:objectHeight;
%         xx = objectRadiusOneAxis * ppval(ppObjectGeometry, objh./objectHeight);
%         plot(xx,objh,'.')
%         hold on
%         plot(distancesInRangePos,heightsInRange,'.')
%         hold on
%         plot(fingerLocationArray,heightsInRange,'.')
%          hold on
%         plot(fingerLocationArray + squeezeLocationArrayNoPenetration(iSqueezeLocation),heightsInRange ,'.')
%         axis equal
%         
        %         squeezeLocationArrayNoPenetration(iSqueezeLocation) = max(distancesInRangePos)-3;
%         squeezeLocationArrayNoPenetrationNeg(iSqueezeLocation) = min(distancesInRangeNeg)+3;
%         squeezeLocationArrayNoPenetration(iSqueezeLocation) = max(distancesInRangePos);
%         squeezeLocationArrayNoPenetrationNeg(iSqueezeLocation) = min(distancesInRangeNeg);
    end
    
end

function [simulationSuccess] = issimulated(simulationResultPath)

        if any(size(dir([simulationResultPath '/*.csv' ]),1))
            simulationSuccess = true;
        else
            simulationSuccess = false;
        end
end

function [ pathSqueezeLocation ]= simulationresultpath(root,iSqueezeLocation,squeezeLocationArray,pressure)

    pathSqueezeLocation = [root '/pressure' num2str(pressure) '/' 'loc' num2str(iSqueezeLocation) '_x' num2str(squeezeLocationArray(1,iSqueezeLocation)) '_y' num2str(squeezeLocationArray(2,iSqueezeLocation)) '_z' num2str(squeezeLocationArray(3,iSqueezeLocation))];
end