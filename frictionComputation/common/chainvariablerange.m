function [ allSamples ] = chainvariablerange( variableRangeArray )
%For the range of each varialbe, chreate a chain matrix of all possible
%variable combinations in it

% number of cells are the number of variables
%     numVariables = numel(variableRangeArray);
%     varLengthArray = zeros(numVariables,1);
%     for iVar = 1:numVariables
%         varLengthArray(iVar) = numel(variableRangeArray{iVar});
%     end
%     numAllSamples = prod(varLengthArray);
%     allSamples = zeros(numAllSamples,numVariables);
%     iSample = 1;
%     for iVar = 1:numVariables
%         for iRange = variableRangeArray{iVar}
%             allSamples(iSample,iVar) = iRange;
%             iSample = iSample + 1;
%         end
%     end

%     allSamples = combvec(variableRangeArray{1},variableRangeArray{2},...
%         variableRangeArray{3},...
%         variableRangeArray{4},variableRangeArray{5})';
    allSamples = combvec(variableRangeArray{1},variableRangeArray{2},...
        variableRangeArray{3},...
        variableRangeArray{4})';
end

