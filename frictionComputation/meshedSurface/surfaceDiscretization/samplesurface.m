function [nodesCoord,nodesCoordMatrix,idxNodeRows] = samplesurface...
        (parametricSurface,rangeU,rangeV,sizeElem,ifVisualize)

%     (parametricSurface,rangeU,rangeV,numEleU,numEleV,surfacePressure,ifVisualize)
% discretize the parametric surface with given number of Elements in
% U and V, which are variables of parametric surface.
% specify the element size to control the number of elements

% First sample the heights/rows, then for each row, evenly sample the
% points
%% determine the variable which controls the height (z axis)
syms u v varHeight varRow real

r = parametricSurface(u,v);
varHeight = symvar(r(3));
zExpr(varHeight) = r(3);
% for parabolloid, the height component is the first one
if numel(varHeight)>1
 varHeight = symvar(r(1));
 zExpr(varHeight) = r(1);
end
if varHeight == u
    varRow = v;
    rangeVarHeight = rangeU;
    rangeVarRow= rangeV;
else
    varRow = u;
    rangeVarHeight = rangeV;
    rangeVarRow= rangeU;
end

rNewParam(varHeight,varRow) = r;
step = (max(rangeVarHeight)-min(rangeVarHeight))/10;
samplesRangeVarHeight = min(rangeVarHeight):step:max(rangeVarHeight);

heightMin = min(double(zExpr(samplesRangeVarHeight)));
heightMax = max(double(zExpr(samplesRangeVarHeight)));

%%
if (sizeElem> heightMax-heightMin)
    warning(['for a better sampling, size of element has to be smaller than the surface height, within the range '...
        num2str(heightMin) ' and ' num2str(heightMax)]);
end
% sizeElemRows = sizeElem./2*sqrt(3);
sizeElemRows = sizeElem./sqrt(3);

sizeElemRows = min(sizeElemRows,heightMax-heightMin);
% numElemRows = ceil((sin(max(rangeU))- sin(min(rangeU)))/sizeElemRows);
numElemRows = max(2,ceil((heightMax- heightMin)/sizeElemRows));
stepHeight = (heightMax- heightMin)/(numElemRows-1);
% monomials = children(zExpr(varHeight));
% constant = 0;
% for m = monomials
%     if isempty(symvar(m))
%         constant = constant + m;
%     end
% end
zExprInverse(varHeight) = finverse(zExpr(varHeight));
% stepSinU = (sin(max(rangeU))- sin(min(rangeU)))/(numElemRows-1);
heightValArray = double(zExprInverse(heightMin:stepHeight:heightMax));
% sample the sphere circle by circle (row by row). the circle is determined
% by the u value. The number of points on each circle is determined by the
% element size. The larger the size, the more number of points on the
% circle. 
uvPairs = [];
idxNodeRows = cell(size(heightValArray));
numCreatedNodes = 0;
for iRow = 1:numel(heightValArray)
    
    curHeight = heightValArray(iRow);
%     circleRadius = abs(cos(curHeight));
    circleParam(varRow) = rNewParam(curHeight,varRow);
    curParam = rNewParam(curHeight,varRow);
    try
        coeffs1 = coeffs(curParam(1));
        coeffs2 = coeffs(curParam(2));
        if ~isempty(coeffs1) && ~isempty(coeffs2)
            circleRadius = sqrt(double(coeffs1(end))^2 + double(coeffs2(end))^2);
        else
            circleRadius = 0;
        end
    catch
        circleRadius = 0;
    end


    if circleRadius < 1e-3
        circleRadius = 0; % means there is only 1 point. put the point inside of the ndoe list.
        valPointOneRow = min(rangeVarRow);
    else
        pointStep = 2*asin(sizeElem./2./circleRadius);
        numElemInCirle = ceil((max(rangeVarRow) - min(rangeVarRow))./pointStep);
        realStep = (max(rangeVarRow) - min(rangeVarRow))./numElemInCirle;
        valPointOneRow = min(rangeVarRow) :realStep: max(rangeVarRow);
    end
    % if circleRadius = 0, means its a point and not a circle. Just add the
    % point into the node list
%     circle = parametricSurface(curHeight,varRow);

    
    

    % check distance between two points. If its too short, select partial
    % of them.
%     if numel(valPointOneRow)>1
%         diffFirstSecondPoint = sum((double(circleParam(valPointOneRow(1))).^2 ...
%         - double(circleParam(valPointOneRow(2)))).^2);
%         stepValCircle = ceil(sizeElem / diffFirstSecondPoint);
%         valPointOneRowLessPoints = valPointOneRow(1:stepValCircle:end);
%     end
%     valPointOneRow= valPointOneRowLessPoints;


    if numel(valPointOneRow)>1
        diffFirstEndPoint = sum((double(circleParam(valPointOneRow(1)))...
    - double(circleParam(valPointOneRow(end)))).^2);
        if diffFirstEndPoint<1e-6
            valPointOneRow = valPointOneRow(1:end-1);
        end
    end
    idxNodeRows{iRow} = numCreatedNodes+1:numCreatedNodes+numel(valPointOneRow);
    numCreatedNodes = numCreatedNodes+numel(valPointOneRow);
    pairs = [repmat(curHeight,size(valPointOneRow)); valPointOneRow];
    uvPairs = [uvPairs pairs];
%     vValArray = [vValArray vValCircle];
end

%%
nodeCartesian = rNewParam(uvPairs(1,:),uvPairs(2,:));

% nodeCartesian = parametricSurface(uvPairs(1,:),uvPairs(2,:));
nodesCoord.x = double(nodeCartesian{1});
nodesCoord.y = double(nodeCartesian{2});
nodesCoord.z = double(nodeCartesian{3});
nodesCoordMatrix = [nodesCoord.x;nodesCoord.y;nodesCoord.z];
center = mean(nodesCoordMatrix')';

nodesCoordMatrix =  nodesCoordMatrix - center;
nodesCoord.x  = nodesCoord.x - center(1);
nodesCoord.y  = nodesCoord.y - center(2);
nodesCoord.z  = nodesCoord.z - center(3);

%%
if ifVisualize
    figure;
    plot3(nodesCoord.x,nodesCoord.y,nodesCoord.z,'bo');
    axis equal;
end

end

