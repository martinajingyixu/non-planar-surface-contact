function [downSampledWrench,idxDownsampledWrench] = downsamplewrench(wrench,numSteps,pathResults)
if nargin == 1
%     numSteps = 3.4;
    numSteps = 2.5; %77 % 60 used for parametric surface
%       numSteps = 2.2; % 60 used for FEM
%       numSteps = 2.1; % 60 used for FEM

%     numSteps = 2.08; % 52 samples
%     numSteps = 1.25; % 52 samples


end
idxDownsampledWrench = [];
if size(wrench,2) == 6
    wrench = wrench';
end
rangeMin = min(wrench');
rangeMax = max(wrench');
stepSz = (rangeMax - rangeMin)./numSteps; % very fine downsample
% step = (rangeMax - rangeMin)./6; % default for showing the results of
% LS 
% step = (rangeMax - rangeMin)./5;

% step = (rangeMax - rangeMin)./3;

stepSz(isnan(stepSz)) = 0;

x1V = rangeMin(1):stepSz(1):rangeMax(1);if isempty(x1V) x1V = 0; end
x2V = rangeMin(2):stepSz(2):rangeMax(2);if isempty(x2V) x2V = 0; end
x3V = rangeMin(3):stepSz(3):rangeMax(3);if isempty(x3V) x3V = 0; end
x4V = rangeMin(4):stepSz(4):rangeMax(4);if isempty(x4V) x4V = 0; end
x5V = rangeMin(5):stepSz(5):rangeMax(5);if isempty(x5V) x5V = 0; end
x6V = rangeMin(6):stepSz(6):rangeMax(6);if isempty(x6V) x6V = 0; end

% x1V(isnan(x1V)) = 0;
% x2V(isnan(x2V)) = 0;
% x3V(isnan(x3V)) = 0;
% x4V(isnan(x4V)) = 0;
% x5V(isnan(x5V)) = 0;
% x6V(isnan(x6V)) = 0;

% x1V(isnan(x1V)) = 0;
% x2V(isnan(x2V)) = 0;
% x3V(isnan(x3V)) = 0;
% x4V(isnan(x4V)) = 0;
% x5V(isnan(x5V)) = 0;
% x6V(isnan(x6V)) = 0;
parfor ix1 = 1:numel(x1V)
    if numel(x1V) <=1 % means no points in this dimension needs to be downsampled
        idx1 = 1:length(wrench);
    else
        idx1 = find(wrench(1,:)>= x1V(ix1) & wrench(1,:)<=x1V(ix1)+stepSz(1));
        if numel(idx1) == 0
            continue
        end
    end
    

    for ix2 = 1:numel(x2V)
%    for ix2 = 1:numel(x2V)
        if numel(x2V) <= 1
            idx2 = idx1;
        else
            idx2 = intersect(idx1,find(wrench(2,:)>= x2V(ix2) & wrench(2,:)<=x2V(ix2)+stepSz(2)));
            if numel(idx2) == 0
                continue
            end
        end
    
        

        for ix3 = 1:numel(x3V)
            if numel(x3V) <= 1
                idx3 = idx2;
            else
                idx3 =  intersect(idx2,find(wrench(3,:)>= x3V(ix3) & wrench(3,:)<=x3V(ix3)+stepSz(3)));
                if numel(idx3) == 0
                    continue
                end
            end
            

            for ix4 = 1:numel(x4V)
                if numel(x4V)<= 1
                    idx4 = idx3;
                else
                   idx4 = intersect(idx3, find(wrench(4,:)>= x4V(ix4) & wrench(4,:)<=x4V(ix4)+stepSz(4)));
                   if numel(idx4) == 0
                       continue
                   end
                end
                

                for ix5 = 1:numel(x5V)
                    
                    if numel(x5V) <=1
                       idx5 = idx4;
                    else
                       idx5=  intersect(idx4,find(wrench(5,:)>= x5V(ix5) & wrench(5,:)<=x5V(ix5)+stepSz(5)));
                       
                        if numel(idx5) == 0
                            continue
                        end
                    end

                    for ix6 = 1:numel(x6V)
                       if numel(x6V) <= 1
                           idx = idx5;
                       else
                           idx=  intersect(idx5,find(wrench(6,:)>= x6V(ix6) & wrench(6,:)<=x6V(ix6)+stepSz(6)));
                       end
                        if ~isempty(idx)
                            % Find the sample which is nearst to the
                            % centerid of the samples in grid
                            [~,idxMinDist] = min(sqrt(sum((wrench(:,idx) -  mean(wrench(:,idx),2)).^2)));
                            idxDownsampledWrench = [idxDownsampledWrench,idx(idxMinDist)];
                        end             
                    end
                end
            end
        end
    end
end
% for iAxis1 = 1:6
%     for iAxis2 = 1:6
%         for iAxis3 = 1:6
% 
%             if iAxis1 < iAxis2 && iAxis2 < iAxis3 
%                 plotDim = [iAxis1,iAxis2,iAxis3];
%                 uniqueIdx = getidxconvexhull(wrench(plotDim,:));
%                 idxArray = [idxArray,uniqueIdx'];             
%             end
%         end
%     end
% end
% idxArray = unique(idxArray);
idxDownsampledWrench = unique(idxDownsampledWrench);
downSampledWrench = wrench(:,idxDownsampledWrench);

if nargin == 3
    save([pathResults '/wrenchTwist.mat'],'downSampledWrench','idxDownsampledWrench','-append');
end

end


