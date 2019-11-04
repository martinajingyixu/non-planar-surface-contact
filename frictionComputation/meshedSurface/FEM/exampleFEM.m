warning('off','all')
PACKAGE_ROOT = '~/3d-surface-contact';

addpath(PACKAGE_ROOT);
% loadconstantsFingerBeam; % type of finger
loadconstantsFingerEllipsoid;

addpackagespaths(PACKAGE_ROOT);

global CONSTANTS
FEMdir = CONSTANTS.DIR.FEM_RESULTS;

[dbDir,objectParameters] = db_parameters(PACKAGE_ROOT);
objectRange =1:height(objectParameters);
pathResults = ['frictionResults/FEM/' CONSTANTS.FINGER_TYPE '/'] ;
if ~exist(pathResults,'dir')
    mkdir(pathResults);
end
%% model generation

config.ifVisualize = 1;
config.ifSampling = 0;
config.ifConvex = 1;
config.lambda = 1; % 100
config.wrenchWeight = 10;
config.twistWeight = 1;
config.samplingStepSize = 10;
config.fineLevel = 0;
config.ifVerbose = 1; % display results
config.ifFitLS4th = 1;
config.ifFitLSEllip = 1;
%%

for idxObject = 1:numel(objectRange)
    iObject = objectRange(idxObject)
    objectInfo = table2struct(objectParameters(iObject, :));
    load([FEMdir '/FEMSimulationResultsObj' num2str(iObject) '.mat'],'FEMResults');
    % load FEM results for squeezing x axis
    isSqueezingX = true;
    disp('Sqeezing X axis');
    [contactInfox,fittingResultsx,wrenchTwistx,objFittingResultsx] = ...
        db_wrenchTwistLSfitting(pathResults,FEMResults.xaxis,isSqueezingX,config);

    contactInfo.xaxis = contactInfox;
    fittingResults.xaxis = fittingResultsx;
    wrenchTwist.xaxis = wrenchTwistx;
    objFittingResults.xaxis = objFittingResultsx;

    if ~objectInfo.symmetric
        disp('Sqeezing Y axis');
        isSqueezingX = false;
        [contactInfoy,fittingResultsy,wrenchTwisty,objFittingResultsy] = ...
            db_wrenchTwistLSfitting(pathResults,FEMResults.yaxis,isSqueezingX,config);
        contactInfo.yaxis = contactInfoy;
        fittingResults.yaxis = fittingResultsy;
        wrenchTwist.yaxis = wrenchTwisty;
        objFittingResults.yaxis = objFittingResultsy;
    end
    save([pathResults '/FEMContactObj' num2str(iObject) '.mat'],'contactInfo',...
        'fittingResults','objFittingResults','config','wrenchTwist');
end
