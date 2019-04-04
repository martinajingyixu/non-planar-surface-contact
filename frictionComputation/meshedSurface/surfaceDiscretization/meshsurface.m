function [elementArray] = meshsurface(nodesCoordMatrix,idxNodeRows, ifConnectFirstLastNode)
%mesh the parametric surface. given the idx of nodes of each row and the
%coordinate of the nodes
% nodesCoordMatrix: a 3xN matrix, each column is a 3d coord of a node.
% output: the elements of the surface, which consist of indices of nodes
% counter-clockwise. A 3xN matrix

if nargin<3
    ifConnectFirstLastNode = false;
end
    

elementArray = [];
% elementArrayU = [];
% elementArrayL = [];

numRows = numel(idxNodeRows);
numNodes = size(nodesCoordMatrix,2);
processedNodeL = zeros(1,numNodes);
widestRow = 1;
numNodesRow = 0;

for iRow = 1:numel(idxNodeRows)
    if numel(idxNodeRows{iRow})>numNodesRow
        widestRow = iRow;
        numNodesRow = numel(idxNodeRows{iRow});
    end
        
end

% if numel(idxNodeRows{1}) >= numel(idxNodeRows{end})
%     ifPyramid = true; % upper thin bottom wide. always build the mesh starting with wide part.
% else
%     ifPyramid = false;
%     idxNodeRows = flip(idxNodeRows);    
% end

% if numel(idxNodeRows{1}) >= numel(idxNodeRows{end})
%     ifPyramid = true; % upper thin bottom wide. always build the mesh starting with wide part.
% else
%     ifPyramid = false;
%     idxNodeRows = flip(idxNodeRows);    
% end



    if widestRow<numRows
    for iRow = widestRow:numRows-1
        if numel(idxNodeRows{iRow}) >= numel(idxNodeRows{iRow+1})
             ifPyramid = true;
        else
             ifPyramid = false;
        end
        % two neighbor nodes create one edge. 
        if numel(idxNodeRows{iRow}) == 1
            continue
        end
            edgeArray = [idxNodeRows{iRow} idxNodeRows{iRow}(1)];
            
        idxHigherRow = iRow + 1;

        for iEdge = 1:numel(edgeArray)-1

            curEdge = [edgeArray(iEdge) edgeArray(iEdge+1)];
            leftPoint = getmiddlepointcoord(nodesCoordMatrix,[edgeArray(iEdge) edgeArray(iEdge)]);
            rightPoint = getmiddlepointcoord(nodesCoordMatrix,[edgeArray(iEdge+1) edgeArray(iEdge+1)]);

            midPoint = getmiddlepointcoord(nodesCoordMatrix,curEdge);
            % find the nearst point in the upper row
            if idxHigherRow <= numRows
                nodeIdxArray = idxNodeRows{idxHigherRow};
                idxNearstNode = getidxnearstpoint(leftPoint,nodeIdxArray,nodesCoordMatrix);
                % counter-clockwise
                if iEdge == numel(edgeArray)-1 && ifConnectFirstLastNode == false

                else
                    if ifPyramid
                        elementArray = [elementArray;  [edgeArray(iEdge),edgeArray(iEdge+1), idxNearstNode]];
%                         elementArray = [elementArray;  [edgeArray(iEdge), idxNearstNode,edgeArray(iEdge+1)]];

                    else
                        elementArray = [elementArray;  [edgeArray(iEdge+1), edgeArray(iEdge), idxNearstNode]];
                    end
                end

                if numel(nodeIdxArray)>1
                    idx = find(nodeIdxArray == idxNearstNode);
                    if ~ifConnectFirstLastNode && idx == 1
                        continue
                    end
                    if idx  >1
                        leftNode = idxNearstNode-1;
                    else
                        leftNode = nodeIdxArray(end);

                    end
                    if processedNodeL(1,leftNode) == 0
                        if ifPyramid
                            elementArray = [elementArray;  [edgeArray(iEdge), idxNearstNode,leftNode]];
%                             elementArray = [elementArray;  [edgeArray(iEdge),leftNode, idxNearstNode]];

                        else
                            elementArray = [elementArray;  [idxNearstNode,edgeArray(iEdge), leftNode]];
                        end
                        processedNodeL(1,leftNode) = 1;
                    end

                end
            end

        end


    end
    end

    %%
    if widestRow>1
    if numel(idxNodeRows{widestRow}) <= numel(idxNodeRows{widestRow-1})
         ifPyramid = true;
    else
         ifPyramid = false;
         idxNodeRows = [flip(idxNodeRows(1:widestRow)),idxNodeRows(widestRow+1:end)];    

    end
    for iRow = 1:widestRow

        % two neighbor nodes create one edge. 
        if numel(idxNodeRows{iRow}) == 1
            continue
        end
            edgeArray = [idxNodeRows{iRow} idxNodeRows{iRow}(1)];
            
        idxHigherRow = iRow + 1;

        for iEdge = 1:numel(edgeArray)-1

            curEdge = [edgeArray(iEdge) edgeArray(iEdge+1)];
            leftPoint = getmiddlepointcoord(nodesCoordMatrix,[edgeArray(iEdge) edgeArray(iEdge)]);
            rightPoint = getmiddlepointcoord(nodesCoordMatrix,[edgeArray(iEdge+1) edgeArray(iEdge+1)]);

            midPoint = getmiddlepointcoord(nodesCoordMatrix,curEdge);
            % find the nearst point in the upper row
            if idxHigherRow <= numRows
                nodeIdxArray = idxNodeRows{idxHigherRow};
                idxNearstNode = getidxnearstpoint(leftPoint,nodeIdxArray,nodesCoordMatrix);
                % counter-clockwise
                if iEdge == numel(edgeArray)-1 && ifConnectFirstLastNode == false

                else
                    if ifPyramid
%                         elementArray = [elementArray;  [edgeArray(iEdge),edgeArray(iEdge+1), idxNearstNode]];
                        elementArray = [elementArray;  [edgeArray(iEdge), idxNearstNode,edgeArray(iEdge+1)]];

                    else
                        elementArray = [elementArray;  [edgeArray(iEdge+1), edgeArray(iEdge), idxNearstNode]];
                    end
                end

                if numel(nodeIdxArray)>1
                    idx = find(nodeIdxArray == idxNearstNode);
                    if ~ifConnectFirstLastNode && idx == 1
                        continue
                    end
                    if idx  >1
                        leftNode = idxNearstNode-1;
                    else
                        leftNode = nodeIdxArray(end);

                    end
                    if processedNodeL(1,leftNode) == 0
                        if ifPyramid
%                             elementArray = [elementArray;  [edgeArray(iEdge), idxNearstNode,leftNode]];
                            elementArray = [elementArray;  [edgeArray(iEdge),leftNode, idxNearstNode]];

                        else
                            elementArray = [elementArray;  [idxNearstNode,edgeArray(iEdge), leftNode]];
                        end
                        processedNodeL(1,leftNode) = 1;
                    end

                end
            end

        end


    end
    end
    
    
elementArray = elementArray';
end

function midPoint = getmiddlepointcoord(nodesCoordMatrix,curEdge)

    midPoint =  (nodesCoordMatrix(:,curEdge(1)) + ...
        nodesCoordMatrix(:,curEdge(2)))./2;

end

function area = trianglearea()

    area = 0.5*sqrt((x2*y3-x3*y2)^2 + (x3*y1-x1*y3)^2 + (x1*y2-x2*y1)^2);

end

function idxNearstNode = getidxnearstpoint(point,nodeIdxArray,nodesCoordMatrix)

% find the nearst point of the given point in the node array and return the idx of the
% nodeArray
    pointArray = repmat(point,[1,numel(nodeIdxArray)]);
    nodeCandidatesCoordArray = nodesCoordMatrix(:,nodeIdxArray);
    distArray = sum((pointArray - nodeCandidatesCoordArray).^2);
    [~,idxMinDist] = min(distArray);
    idxNearstNode = nodeIdxArray(idxMinDist);
end

