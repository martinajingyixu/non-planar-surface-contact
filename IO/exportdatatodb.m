function [FEMSimulationResults] = exportdatatodb(dbDir,objectParameters,objectRange)
%exportdatatodb Export the simulation results to .mat file
%read data from the directory of dirSimulationResult, then export the .mat
%file to dirDataExport

global CONSTANTS
if ~exist('objectRange', 'var') || isempty(objectRange)
    objectRange = 1:height(objectParameters); end

FEMSimulationResults = cell(numel(objectRange),1);
load([dbDir 'dbSimulationSuccess.mat'],'dbSimulationSuccess');


if ~exist(CONSTANTS.DIR.FEM_RESULTS, 'dir')
    mkdir(CONSTANTS.DIR.FEM_RESULTS); end
for idxObject = 1:numel(objectRange)
    iObject = objectRange(idxObject);
    iObject
    objectInfo = table2struct(objectParameters(iObject, :));
    load([[CONSTANTS.DIR.SIMULATION '/' objectInfo.filename '/'] 'squeezeLocationArrayNoPenetration.mat'],'trueSqueezeLocationArray');
    simulationObjectRoot = [CONSTANTS.DIR.SIMULATION '/' objectInfo.filename '/'];

    FEMSimulationResults{iObject}.xaxis = cell(length(trueSqueezeLocationArray.xaxis),length(CONSTANTS.SIMULATION.PRESSURE_ARRAY));
    simulationObjectRootXAxis = [simulationObjectRoot 'xaxis/'];

    if ~objectInfo.symmetric
        FEMSimulationResults{iObject}.yaxis = cell(length(trueSqueezeLocationArray.yaxis),length(CONSTANTS.SIMULATION.PRESSURE_ARRAY));
        simulationObjectRootYAxis = [simulationObjectRoot 'yaxis/'];
    end

    for iPressure = 1:length(CONSTANTS.SIMULATION.PRESSURE_ARRAY)
%     for iPressure = 1

        for iSqueezeLocation = 1:length(trueSqueezeLocationArray.xaxis)
%         for iSqueezeLocation = 1
            simulationResultPath = simulationresultpath(simulationObjectRootXAxis,...
        iSqueezeLocation,trueSqueezeLocationArray.xaxis,CONSTANTS.SIMULATION.PRESSURE_ARRAY(iPressure));
             if dbSimulationSuccess{iObject}.xaxis(iSqueezeLocation,iPressure) == 1
                 simulationData = loadCSVFiles(simulationResultPath);
                 try
                    simulationData.element_pressure = loadpressurefromtxtfile(simulationResultPath);
                 catch
                     warning('can not load text with importdata, use scan txt data instead');
   
                    simulationData.element_pressure = loadpressurefromtxtfilescantext(simulationResultPath,length(simulationData.element_info));
                 end
                 msg = ['Error in XAxis: Pressure data of Object ' int2str(iObject) ' Location ' int2str(iSqueezeLocation) ' Pressure ' int2str(iPressure) ' has problem!'];
                 assert(length(simulationData.element_info) == length(simulationData.element_pressure),msg);
                 FEMSimulationResults{iObject}.xaxis(iSqueezeLocation,iPressure) = {simulationData};

             end


        end


        if ~objectInfo.symmetric
            for iSqueezeLocation = 1:length(trueSqueezeLocationArray.yaxis)
                simulationResultPath = simulationresultpath(simulationObjectRootYAxis,...
            iSqueezeLocation,trueSqueezeLocationArray.yaxis,CONSTANTS.SIMULATION.PRESSURE_ARRAY(iPressure));
                 if dbSimulationSuccess{iObject}.yaxis(iSqueezeLocation,iPressure) == 1
                     simulationData = loadCSVFiles(simulationResultPath);
                 try
                    simulationData.element_pressure = loadpressurefromtxtfile(simulationResultPath);
                 catch
                     warning('can not load text with importdata, use scan txt data instead');
   
                    simulationData.element_pressure = loadpressurefromtxtfilescantext(simulationResultPath,length(simulationData.element_info));
                 end
                    msg = ['Error in YAxis: Pressure data of Object ' int2str(iObject) ' Location ' int2str(iSqueezeLocation) ' Pressure ' int2str(iPressure) ' has problem!'];
                     assert(length(simulationData.element_info) == length(simulationData.element_pressure),msg);
                     FEMSimulationResults{iObject}.yaxis(iSqueezeLocation,iPressure) = {simulationData};
                 end


            end
        end
        FEMResults = FEMSimulationResults{iObject};
        save([CONSTANTS.DIR.FEM_RESULTS '/FEMSimulationResultsObj' num2str(iObject) '.mat'],'FEMResults');
    end

%     save([CONSTANTS.DIR.FEM_RESULTS '/FEMSimulationResults.mat'],'FEMSimulationResults');

end 
% savefilesforeachobject(CONSTANTS.DIR.FEM_RESULTS, FEMSimulationResults)

end

function [simulationData] = loadCSVFiles(dirSimulationResult)

%% read available results

fileNamesCSV = {'element_area','element_info','element_kinematic_energy',...
    'element_mat','element_real','element_stiffness_energy',...
    'nodal_coord','nodal_displacement'};


simulationData = struct();
% element_area always can not import because the data is too small. So
% use import data to get all the non NaN data.
element_area_rough = importdata([dirSimulationResult '/element_area.csv']);
simulationData.element_area_rough = element_area_rough.data;

for iFileName = 1:length(fileNamesCSV)
    try
        varname = genvarname(fileNamesCSV{iFileName});
        tmpFileName = [dirSimulationResult '/' fileNamesCSV{iFileName} '.csv'];
        eval(['simulationData.' varname '=csvread( tmpFileName);']);

    catch
%         if strcmp(fileNamesCSV{iFileName},'nodal_coord') ||...
%                 strcmp(fileNamesCSV{iFileName},'nodal_displacement')
%             error(['Can not read the csv file of nodal info: ', fileNamesCSV{iFileName}]);
%         else
            warning(['Can not read the csv file: ', tmpFileName]);
%         end
    end

end

end

function [pressureArray] = loadpressurefromtxtfile(dirSimulationResult)

rawData = importdata([dirSimulationResult '/elementPressure.txt']);
try
rawTextData= rawData.textdata(:,1:2);
catch
    disp(dirSimulationResult)
    return
end
%     nonEmptyelements = find(~cellfun(@isempty,rawTextData(:,2))); % find non empty cell
tempNumbersArray = cell(length(rawTextData),2);

% search for lines which has the format '%d %.10f\n'
for iLine = 1:length(rawTextData)
    tmpString = [rawTextData{iLine,1} ' ' rawTextData{iLine,2}];
    tempNumbersArray(iLine,1:2) = textscan(tmpString,'%f %.10f\n');
end

idxNonEmptyCol1 = find(~cellfun(@isempty,tempNumbersArray(:,1)));
idxNonEmptyCol2 = find(~cellfun(@isempty,tempNumbersArray(:,2)));

idxBothNonEmpty = intersect(idxNonEmptyCol1,idxNonEmptyCol2);
pressureArray = cell2mat(tempNumbersArray(idxBothNonEmpty,:));
end

function [pressureArray] = loadpressurefromtxtfilescantext(dirSimulationResult,sizeElement)

fp = fopen([dirSimulationResult '/elementPressure.txt']);
pressureArray = zeros(sizeElement,2);
iElement = 1;
if (fp == -1)
    disp(['can not open ' dirSimulationResult '/elementPressure.txt']);
    return
end

 % find the text ELEPRES, read element until a break line
while ~feof(fp)
    %%
    str = fgets(fp);
    textStart = strfind(str,'ELEPRES');
    %%
    if ~isempty(textStart)
        str = fgets(fp);
        while true

            elementIdPressure= textscan(str,'%f %.10f \n');
            if ~isempty(elementIdPressure{1,1}) && ~isempty(elementIdPressure{1,2})
                pressureArray(iElement,:) = [elementIdPressure{1,1} elementIdPressure{1,2}];
                iElement = iElement + 1;
            end
            str = fgets(fp);
            if length(str)==81 || feof(fp) % read until a line break
                break;
            end 
     end

    end
end 
fclose(fp);
end

function [ pathSqueezeLocation ]= simulationresultpath(root,iSqueezeLocation,squeezeLocationArray,pressure)

pathSqueezeLocation = [root '/pressure' num2str(pressure) '/' 'loc' num2str(iSqueezeLocation) '_x' num2str(squeezeLocationArray(1,iSqueezeLocation)) '_y' num2str(squeezeLocationArray(2,iSqueezeLocation)) '_z' num2str(squeezeLocationArray(3,iSqueezeLocation))];

end

function [] = savefilesforeachobject(dbDir,FEMSimulationResults)
% save the FEM simulation restuls for each object separately such the file
% of each object can be loaded independently.
    numObjects = numel(FEMSimulationResults);
    save([dbDir '/FEMSimulationResults.mat'],'-append','numObjects')

    for iObject = 1:numObjects
        v = genvarname(['FEMSimulationResultsObj' int2str(iObject)]);
        eval([v '= FEMSimulationResults{iObject};']);
        save([dbDir '/FEMSimulationResults.mat'],'-append',v)

    end
end
