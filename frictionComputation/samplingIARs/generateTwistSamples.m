function [twistSampleVariables] = generateTwistSamples(frictionCenter,lengthArray,fineLevel,pathTwist)

if nargin<3
    fineLevel = 2;
end
l = mean(lengthArray);

if fineLevel == 6
    coeff = unique([0.01,0.05,0.1:0.1:0.3,1]);
    numPointsAngle = 500;
    pitchSamples = [-0.5,-0.25,0.25,0.5].*l;
    numAngles = 25;
elseif fineLevel == 5
    coeff = unique([0.1,0.2,0.3,1,0.5,0.1:0.025:0.5,0:0.05:0.1,0.3:0.01:0.4]); % 36
    numPointsAngle = 10000;
    pitchSamples = [0.25].*l;
    numAngles = 50;
elseif fineLevel == 4
    coeff = unique([0.01,0.05,0.1,0.2,0.3:0.02:0.4,0.5,1]); % 36
    numPointsAngle = 10000;
    pitchSamples = [0.25].*l;
    numAngles = 50;
elseif fineLevel == 3
    coeff = unique([0.01,0.05:0.05:0.1,0.1:0.1:0.3,0.3:0.025:0.4,0.5,1]); % 20
    numPointsAngle = 5000;
    pitchSamples = [0.25,0.5].*l;
    numAngles = 16;    
elseif fineLevel == 2
    coeff = unique([0.01,0.05,0.1,0.1:0.1:0.5,1]); 
    numPointsAngle = 5000;
    pitchSamples = [0.125,0.25,0.5,1,2].*l;
    numAngles = 16;
elseif fineLevel == 1
    coeff = unique([0.01,0.05,0.1,0.1:0.1:0.5,1]);
    numPointsAngle = 5000;
    pitchSamples = [0.25].*l;
    numAngles = 30;
elseif fineLevel == 0
    coeff = unique([0.01,0.05,0.1:0.1:0.3,1]);
    numPointsAngle = 500;
    pitchSamples = [0.25].*l;
    numAngles = 25;
elseif fineLevel == -1
    coeff = 0;
    numPointsAngle = 10;
    pitchSamples = 0;
    numAngles = 0;    
end
pitchSamples = unique([pitchSamples,0,-pitchSamples]);
%% test settings
distHalf = [coeff.*min(lengthArray),coeff.*max(lengthArray)];
dist1Samples = unique(floor([-distHalf,0,distHalf].*10000))./10000;
dist2Samples = dist1Samples;

%%
[angle,omega] = samplesphere(numAngles);

combi = combvec(1:length(angle),dist1Samples,dist2Samples)';

TwistSamples = [angle(combi(:,1),:),combi(:,end-1:end)];
TwistSamplesCoord = [omega(combi(:,1),:),combi(:,end-1:end)];


[TwistSamplesLimit,TwistSamplesCoordLimit]=getlimitTwistPoints(numPointsAngle);

TwistSamplesCoord = [TwistSamplesCoord;TwistSamplesCoordLimit];
TwistSamples = [TwistSamples;TwistSamplesLimit];

numSamples = size(TwistSamples,1);
TwistSamples4Integral = zeros(numSamples,6);
% TwistSamples4Integral = [];

%%
for iRow = 1:numSamples
    iRow;
    curSample = TwistSamples(iRow,:);
    dist1V = curSample(3);
    dist2V= curSample(4);
    a = TwistSamplesCoord(iRow,1);
    b = TwistSamplesCoord(iRow,2);
    c = TwistSamplesCoord(iRow,3);

    [~,idxField] = max(abs(TwistSamplesCoord(iRow,1:3)));
    if idxField ==1
        PintLA = frictionCenter + [-(dist1V*b+dist2V*c)/a;dist1V;dist2V];  

        x0V = PintLA(1);
        y0V = PintLA(2);
        z0V = PintLA(3);
    elseif idxField ==2
        PintLB = frictionCenter + [dist1V;-(dist1V*a+dist2V*c)/b;dist2V];  

        x0V = PintLB(1);
        y0V = PintLB(2);
        z0V = PintLB(3);
    elseif idxField==3

        PintLC = frictionCenter + [dist1V;dist2V;-(dist1V*a+dist2V*b)/c];  

        x0V = PintLC(1);
        y0V = PintLC(2);
        z0V = PintLC(3);
    end
    TwistSamples4Integral(iRow,:) = double([TwistSamplesCoord(iRow,1:3),x0V,y0V,z0V]);
end

twistSampleVariables = combvec((1:length(TwistSamples4Integral)),pitchSamples)';
twistSampleVariables = [TwistSamples4Integral(twistSampleVariables(:,1),:),twistSampleVariables(:,end)]';
if nargin == 4
    save([pathTwist '/twistSampleVariables.mat'],'coeff',...
        'numPointsAngle','pitchSamples','l','numAngles','twistSampleVariables');
end
end

function [TwistSamplesLimit,TwistSamplesCoordLimit]=getlimitTwistPoints(numPoints)
% get the limit axes. X Y Z-axes with infinte distance, and axes with zero
% distance

angleLimit = [0,0;0,pi;0,pi/2;0,-pi/2;pi/2,0;-pi/2,0];
omegaLimit = [1,0,0;-1,0,0;0,1,0;0,-1,0;0,0,1;0,0,-1];

combiLimit = combvec(1:length(angleLimit),[0,10000000,-10000000],[0,10000000,-10000000])';
combiLimit(find(combiLimit(:,2) == combiLimit(:,3)),:) = [];
TwistLimitSamples = [angleLimit(combiLimit(:,1),:),combiLimit(:,end-1:end)];
TwistSamplesCoordLimit = [omegaLimit(combiLimit(:,1),:),combiLimit(:,end-1:end)];

[angle0dist,omega0dist] = samplesphere(numPoints);
dist = [0];
combi0dist = combvec(1:length(angle0dist),dist,dist)';
TwistSamples0dist = [angle0dist(combi0dist(:,1),:),combi0dist(:,end-1:end)];
TwistSamplesCoord0dist = [omega0dist(combi0dist(:,1),:),combi0dist(:,end-1:end)];

TwistSamplesCoordLimit = [TwistSamplesCoord0dist;TwistSamplesCoordLimit];
TwistSamplesLimit = [TwistSamples0dist;TwistLimitSamples];

end