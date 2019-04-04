function [sampledWrench,shapeArray]=samplefittedlimitsurface(wrench,fittedCoeff,samplingStepSize,order)
 %%   
    maxd = max(abs(wrench'))*1.2;
    step = maxd / samplingStepSize;
%     maxd = max(max(wrench));
%     l1 = 1.1;
%     l2 = 0.7;
%     step1 = (l1*maxd-l2*maxd)/2;
%     step2 = l2*maxd/7;
%     disthalf = [0:step2:l1*maxd,l1*maxd:step1:l2*maxd];
%     valueArray =unique([-disthalf,disthalf]);

    
    shapeArray = [];
    sampledWrench = [];
    if order ==2
    fittedCoeffRearranged = [fittedCoeff(1,1),fittedCoeff(2,2),fittedCoeff(3,3),...
    fittedCoeff(4,4),fittedCoeff(5,5),fittedCoeff(6,6),...
    2*fittedCoeff(1,2),2*fittedCoeff(1,3),2*fittedCoeff(1,4),...
    2*fittedCoeff(1,5),2*fittedCoeff(1,6),2*fittedCoeff(2,3),...
    2*fittedCoeff(2,4),2*fittedCoeff(2,5),2*fittedCoeff(2,6),...
    2*fittedCoeff(3,4),2*fittedCoeff(3,5),2*fittedCoeff(3,6),...
    2*fittedCoeff(4,5),2*fittedCoeff(4,6),2*fittedCoeff(5,6)];

    elseif order ==4
        fittedCoeffRearranged = fittedCoeff';
    end
    % use partial combvec and partial for loop for memory and time
    % efficiency
    disp('sample LS')
    
    %%
%     tic
%     valueArray = -maxd:step:maxd;
%     combi3d = combvec(valueArray,valueArray);
%     combiLength = length(combi3d);
%     
%     for xx1 = valueArray
%         for xx2 = valueArray
%             for xx3  = valueArray
%                 for xx4  = valueArray
% 
%                 shapeArrayTemp = computeshape2(fittedCoeffRearranged,...
%                     repmat(xx1,[1,combiLength]),repmat(xx2,[1,combiLength]),...
%                     repmat(xx3,[1,combiLength]),...
%                     repmat(xx4,[1,combiLength]),...
%                     combi3d(1,:),combi3d(2,:));
%                 idx = find(shapeArrayTemp<1.01);
%                 shapeArray = [shapeArray,shapeArrayTemp(idx)];
%                 numIdx = numel(idx);
% %                 sampledWrench = [sampledWrench,[repmat([xx1;xx2;xx3],[1,numIdx]);...
% %                     combi3d(:,idx)]];
%                 sampledWrench = [sampledWrench,[repmat(xx1,[1,numIdx]);...
%                     repmat(xx2,[1,numIdx]);repmat(xx3,[1,numIdx]);...
%                     repmat(xx4,[1,numIdx]);...
%                     combi3d(:,idx)]];
%                 end
%             end
%         end
%     end
%     size(sampledWrench)
%     toc
 
    %%
    valueArray = -maxd:step:maxd;
    combi3d = combvec(valueArray,valueArray,valueArray);
    combiLength = length(combi3d);
    
    parfor iv = 1:numel(valueArray)
%         iv
        xx1 = valueArray(iv);
        for xx2 = valueArray
            for xx3  = valueArray
                shapeArrayTemp = computeshape(fittedCoeffRearranged,order,...
                    repmat(xx1,[1,combiLength]),repmat(xx2,[1,combiLength]),...
                    repmat(xx3,[1,combiLength]),...
                    combi3d(1,:),combi3d(2,:),combi3d(3,:));
                idx = find(abs(shapeArrayTemp-1)<0.01);
                shapeArray = [shapeArray,shapeArrayTemp(idx)];
                numIdx = numel(idx);
%                 sampledWrench = [sampledWrench,[repmat([xx1;xx2;xx3],[1,numIdx]);...
%                     combi3d(:,idx)]];
                sampledWrench = [sampledWrench,[repmat(xx1,[1,numIdx]);...
                    repmat(xx2,[1,numIdx]);repmat(xx3,[1,numIdx]);...
                    combi3d(:,idx)]];
            end
        end
    end
%     size(sampledWrench)

        
    %%
%     valueArray = -maxd:step:maxd;
%     combi4d = combvec(valueArray,valueArray,valueArray,valueArray);
%     combiLength = length(combi4d);
%     shapeArray = [];
%     sampledWrench = [];
%     for xx1 = valueArray
%         for xx2 = valueArray
%             shapeArrayTemp = computeshape2(fittedCoeffRearranged,...
%                 repmat(xx1,[1,combiLength]),repmat(xx2,[1,combiLength]),...
%                 combi4d(1,:),combi4d(2,:),combi4d(3,:),combi4d(4,:));
%             idx = find(abs(shapeArrayTemp-1)<1e-2);
%             shapeArray = [shapeArray,shapeArrayTemp(idx)];
%             numIdx = numel(idx);
% %             sampledWrench = [sampledWrench,[repmat([xx1;xx2],[1,numIdx]);...
% %                 combi4d(:,idx)]];
%             sampledWrench = [sampledWrench,[repmat(xx1,[1,numIdx]);...
%                 repmat(xx2,[1,numIdx]);...
%                 combi4d(:,idx)]];
%         end
%     end
%     toc
%     size(sampledWrench)

%     parfor iV = 1:numel(valueArray)
%         x1 = valueArray(iV);
%         for x2 = valueArray
%             for x3 =valueArray
%                 for x4 = valueArray
%                     for x5 = valueArray
%                         for x6 = valueArray
%                             shape = computeshape(fittedCoeffRearranged,order);
%                             if abs(shape-1)<1e-2
%                                 sampledWrench = [sampledWrench,[x1;x2;x3;x4;x5;x6]];
%                                 shapeArray = [shapeArray;shape];
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
%     if ifComputeHullValue
%         disp('compute convex hull')
%         [~,hullVolume] = convhulln(sampledWrench');
%     end

   
end

function shapeArray = computeshape2(fittedCoeffRearranged,x1,x2,x3,x4,x5,x6)

    monomials ={  x1.^2; x2.^2;x3.^2;x4.^2; x5.^2; x6.^2;...
        x1.*x2; x1.*x3; x1.*x4; x1.*x5; x1.*x6; ...
        x2.*x3; x2.*x4; x2.*x5; x2.*x6;  x3.*x4; x3.*x5; ...
        x3.*x6;  x4.*x5; x4.*x6;  x5.*x6; };

    shapeArray  = 0;
    for iTerm = 1:size(monomials,1)
        shapeArray = shapeArray + fittedCoeffRearranged(iTerm) * monomials{iTerm};
    end

%     monomials = [x1.^2; x2.^2;x3.^2;x4.^2; x5.^2; x6.^2;...
%         x1.*x2; x1.*x3; x1.*x4; x1.*x5; x1.*x6; ...
%         x2.*x3; x2.*x4; x2.*x5; x2.*x6;  x3.*x4; x3.*x5; ...
%         x3.*x6;  x4.*x5; x4.*x6;  x5.*x6; ];
%     shapeArray = fittedCoeffRearranged*monomials;
end

function shapeArray = computeshape(fittedCoeffRearranged,order,x1,x2,x3,x4,x5,x6)
    if order == 2
    monomials ={  x1.^2; x2.^2;x3.^2;x4.^2; x5.^2; x6.^2;...
        x1.*x2; x1.*x3; x1.*x4; x1.*x5; x1.*x6; ...
        x2.*x3; x2.*x4; x2.*x5; x2.*x6;  x3.*x4; x3.*x5; ...
        x3.*x6;  x4.*x5; x4.*x6;  x5.*x6; };



    elseif order == 4

          monomials ={ x1.^4; x1.^3.*x2; x1.^3.*x3; x1.^3.*x4; x1.^3.*x5;...
            x1.^3.*x6; x1.^2.*x2.^2; x1.^2.*x2.*x3; x1.^2.*x2.*x4; x1.^2.*x2.*x5;...
            x1.^2.*x2.*x6; x1.^2.*x3.^2; x1.^2.*x3.*x4; x1.^2.*x3.*x5; ...
            x1.^2.*x3.*x6; x1.^2.*x4.^2; x1.^2.*x4.*x5; x1.^2.*x4.*x6; ...
            x1.^2.*x5.^2; x1.^2.*x5.*x6; x1.^2.*x6.^2; x1.*x2.^3; ...
            x1.*x2.^2.*x3; x1.*x2.^2.*x4; x1.*x2.^2.*x5; x1.*x2.^2.*x6; ...
            x1.*x2.*x3.^2; x1.*x2.*x3.*x4; x1.*x2.*x3.*x5; x1.*x2.*x3.*x6; ...
            x1.*x2.*x4.^2; x1.*x2.*x4.*x5; x1.*x2.*x4.*x6; x1.*x2.*x5.^2; ...
            x1.*x2.*x5.*x6; x1.*x2.*x6.^2; x1.*x3.^3; x1.*x3.^2.*x4; ...
            x1.*x3.^2.*x5; x1.*x3.^2.*x6; x1.*x3.*x4.^2; x1.*x3.*x4.*x5;...
            x1.*x3.*x4.*x6; x1.*x3.*x5.^2; x1.*x3.*x5.*x6; x1.*x3.*x6.^2; ...
            x1.*x4.^3; x1.*x4.^2.*x5; x1.*x4.^2.*x6; x1.*x4.*x5.^2; ...
            x1.*x4.*x5.*x6; x1.*x4.*x6.^2; x1.*x5.^3; x1.*x5.^2.*x6; ...
            x1.*x5.*x6.^2; x1.*x6.^3; x2.^4; x2.^3.*x3; x2.^3.*x4; x2.^3.*x5; ...
            x2.^3.*x6; x2.^2.*x3.^2; x2.^2.*x3.*x4; x2.^2.*x3.*x5; x2.^2.*x3.*x6;...
            x2.^2.*x4.^2; x2.^2.*x4.*x5; x2.^2.*x4.*x6; x2.^2.*x5.^2; ...
            x2.^2.*x5.*x6; x2.^2.*x6.^2; x2.*x3.^3; x2.*x3.^2.*x4; ...
            x2.*x3.^2.*x5; x2.*x3.^2.*x6; x2.*x3.*x4.^2; x2.*x3.*x4.*x5; ...
            x2.*x3.*x4.*x6; x2.*x3.*x5.^2; x2.*x3.*x5.*x6; x2.*x3.*x6.^2; ...
            x2.*x4.^3; x2.*x4.^2.*x5; x2.*x4.^2.*x6; x2.*x4.*x5.^2; ...
            x2.*x4.*x5.*x6; x2.*x4.*x6.^2; x2.*x5.^3; x2.*x5.^2.*x6; ...
            x2.*x5.*x6.^2; x2.*x6.^3; x3.^4; x3.^3.*x4; x3.^3.*x5; x3.^3.*x6; ...
            x3.^2.*x4.^2; x3.^2.*x4.*x5; x3.^2.*x4.*x6; x3.^2.*x5.^2;...
            x3.^2.*x5.*x6; x3.^2.*x6.^2; x3.*x4.^3; x3.*x4.^2.*x5; ...
            x3.*x4.^2.*x6; x3.*x4.*x5.^2; x3.*x4.*x5.*x6; x3.*x4.*x6.^2; ...
            x3.*x5.^3; x3.*x5.^2.*x6; x3.*x5.*x6.^2; x3.*x6.^3; x4.^4; ...
            x4.^3.*x5; x4.^3.*x6; x4.^2.*x5.^2; x4.^2.*x5.*x6; x4.^2.*x6.^2; ...
            x4.*x5.^3; x4.*x5.^2.*x6; x4.*x5.*x6.^2; x4.*x6.^3; x5.^4; ...
            x5.^3.*x6; x5.^2.*x6.^2; x5.*x6.^3; x6.^4 };

    end
        shapeArray  = 0;
    for iTerm = 1:size(monomials,1)
        shapeArray = shapeArray + fittedCoeffRearranged(iTerm) * monomials{iTerm};
    end
end



