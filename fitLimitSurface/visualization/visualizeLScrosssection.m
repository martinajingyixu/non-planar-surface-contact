function visualizeLScrosssection(wrench,twist,fittedCoeff,order,sampledLS)

if order == 4
    figpath = '/home/jingyi/promotion/programming/softcontactanalysis/figures/LS4thOrder/descretizedEllipsoid/';
elseif order ==2
    figpath = '/home/jingyi/promotion/programming/softcontactanalysis/figures/LSEllpisoid/descretizedEllipsoid/';
end
if ~exist(figpath,'dir')
    mkdir(figpath);
end
    if nargin == 4
        sampledLS = [];
    end
    if size(wrench,2) == 6
        wrench = wrench';
    end
    if size(twist,2) == 6
        twist = twist';
    end

    maxd = max(max(abs(wrench')))*2;
    step = maxd / 30;

    [ fxx, fyy, fzz ] = meshgrid( -maxd:step:maxd, -maxd:step:maxd, -maxd:step:maxd);
    wrenchgrid = {zeros(size(fxx)),zeros(size(fxx)),zeros(size(fxx)),...
        zeros(size(fxx)),zeros(size(fxx)),zeros(size(fxx))};    
    
    p = zeros(1,4);
    
    for iAxis1 = 1:6
        for iAxis2 = 1:6
            for iAxis3 = 1:6
                if iAxis1 < iAxis2 && iAxis2 < iAxis3 
                    
                    plotDim = [iAxis1,iAxis2,iAxis3];
%                     plotDim = [1,4,5];

                    zeroDimArray = setdiff(1:6,plotDim);

                    curwrenchgrid = wrenchgrid;
                    curwrenchgrid{plotDim(1)} = fxx;
                    curwrenchgrid{plotDim(2)} = fyy;
                    curwrenchgrid{plotDim(3)} = fzz;
                    shape = computeshape(fittedCoeff,curwrenchgrid,order);
                                       
                    %% way1: compute isosurface with meshgrid
%                     figHandle = figure('units','normalized','outerposition',[0 0 1 1]);
                    figHandle = figure;
                    axis tight;
                    p(1) = patch(isosurface(fxx,fyy,fzz,shape,1));
                    set( p(1), 'FaceColor', '[0.2510, 0.8863, 1]','FaceAlpha', 0.8, 'EdgeColor', 'none' );

                    view(50, 40);

                    camlight
                    hold on
                    %% way2: compute isosurface with manuelly
                    % In sampled limit surface, find 3 axis which are 0.
                    % The rest are the point to plot
                    if ~isempty(sampledLS) 
                        thresLS = 0.1;
                        zeroDimArray = setdiff(1:6,plotDim); % plot the crosssection of 6D LS
                        idx1 = find(abs(sampledLS(zeroDimArray(1),:))<thresLS);
                        idx2 = find(abs(sampledLS(zeroDimArray(2),:))<thresLS);
                        idx3 = find(abs(sampledLS(zeroDimArray(3),:))<thresLS);
                        idx = intersect(intersect(idx1,idx2),idx3);
                        wrench2plot = sampledLS(plotDim,idx);
                        k = convhulln(wrench2plot');
                        p=patch('Faces',k,'Vertices',wrench2plot');
                        set( p, 'FaceColor', 'b','FaceAlpha', 0.1, 'EdgeColor', 'none' );
                        hold on
                        camlight
                        hold on
                    end
                    
                    %% Plot the samples on the cross section.
                    % Find the samples where other 3 dimensions are 0.
%                     threshold = 0.08;
                    threshold = 0.1;

                    idxD1 = find(abs(wrench(zeroDimArray(1),:))<threshold);
                    idxD2 = find(abs(wrench(zeroDimArray(2),:))<threshold);
                    idxD3 = find(abs(wrench(zeroDimArray(3),:))<threshold);
                    idx = intersect(intersect(idxD1,idxD2),idxD3);
                    % downsample the wrench to be plotted
                    wrench2plotdownsample = zeros(6,numel(idx));
                    wrench2plotdownsample(plotDim(1:3),:) = wrench(plotDim(1:3),idx);
                    [~,idxd] = downsamplewrench(wrench2plotdownsample,5);
                    idx = idx(idxd);
                    colorWrench = '[1,0.7176,0.2706]';
                    color1 = '[0.7843,0.9098,0]';
                    color2 = '[0.8706205,0.4118,0.9098]';
                    if norm(twist) == 0
                        plot3( wrench(plotDim(1),idx),wrench(plotDim(2),idx),wrench(plotDim(3),idx)...
                            , 'Marker', 'o', 'MarkerFaceColor', 'r', 'Markersize', 8, 'LineStyle', 'none')
                        hold on
                    else
                       
                        if order == 2
                            predictedTwist = fittedCoeff * wrench(:,idx);
                        elseif order ==4
                            predictedTwist = computegradient(fittedCoeff,wrench(:,idx))';
                        end

                        hold on

                        visualizewrenchtwistpredictedtwist(wrench(plotDim,idx),figHandle, ...
                            twist(plotDim,idx), predictedTwist(plotDim,:),colorWrench,color1,color2);


                    end

                    xlabelTitle = getlabel(iAxis1);
                    ylabelTitle = getlabel(iAxis2);
                    zlabelTitle = getlabel(iAxis3);
 
                    xlabel(xlabelTitle,'FontSize',34,'interpreter','latex');
                    ylabel(ylabelTitle,'FontSize',34,'interpreter','latex');
                    zlabel(zlabelTitle,'FontSize',34,'interpreter','latex');


%                     print(gcf, [figpath 'fig-with-Legend' num2str(figname) '.pdf'], '-dpdf','-bestfit');
%                     savefig([figpath 'fig-Legend' num2str(figname) '.fig'])
%                     figname = figname+1;
                    
                end
            end
        end
    end  
end

function shape = computeshape(fittedCoeff,curwrenchgrid,order)
x1 = curwrenchgrid{1};
x2 = curwrenchgrid{2};
x3 = curwrenchgrid{3};
x4 = curwrenchgrid{4};
x5 = curwrenchgrid{5};
x6 = curwrenchgrid{6};

if order == 2
    fittedCoeffRearranged = [fittedCoeff(1,1),fittedCoeff(2,2),fittedCoeff(3,3),...
    fittedCoeff(4,4),fittedCoeff(5,5),fittedCoeff(6,6),...
    2*fittedCoeff(1,2),2*fittedCoeff(1,3),2*fittedCoeff(1,4),...
    2*fittedCoeff(1,5),2*fittedCoeff(1,6),2*fittedCoeff(2,3),...
    2*fittedCoeff(2,4),2*fittedCoeff(2,5),2*fittedCoeff(2,6),...
    2*fittedCoeff(3,4),2*fittedCoeff(3,5),2*fittedCoeff(3,6),...
    2*fittedCoeff(4,5),2*fittedCoeff(4,6),2*fittedCoeff(5,6)];

    monomials ={  x1.^2; x2.^2;x3.^2;x4.^2; x5.^2; x6.^2;...
        x1.*x2; x1.*x3; x1.*x4; x1.*x5; x1.*x6; ...
        x2.*x3; x2.*x4; x2.*x5; x2.*x6;  x3.*x4; x3.*x5; ...
        x3.*x6;  x4.*x5; x4.*x6;  x5.*x6; };
    
elseif order == 4
    fittedCoeffRearranged = fittedCoeff;
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
    x5.^3.*x6; x5.^2.*x6.^2; x5.*x6.^3; x6.^4};
end

shape  = 0;
for iTerm = 1:size(monomials,1)
    shape = shape + fittedCoeffRearranged(iTerm) * monomials{iTerm};
end

end

function predictedTwist = computegradient(fittedCoeff,wrench)

    x1 = wrench(1,:);
    x2 = wrench(2,:);
    x3 = wrench(3,:);
    x4 = wrench(4,:);
    x5 = wrench(5,:);
    x6 = wrench(6,:);
    degree = 4;
    numDim = 6;
    numMonomials = nchoosek(degree + numDim-1,numDim-1);
    
    zeroArray = zeros(size(x1));
    diffX1 = [ 4.*x1.^3; 3.*x1.^2.*x2; 3.*x1.^2.*x3; 3.*x1.^2.*x4; 3.*x1.^2.*x5; 3.*x1.^2.*x6; 2.*x1.*x2.^2; 2.*x1.*x2.*x3; 2.*x1.*x2.*x4; 2.*x1.*x2.*x5; 2.*x1.*x2.*x6; 2.*x1.*x3.^2; 2.*x1.*x3.*x4; 2.*x1.*x3.*x5; 2.*x1.*x3.*x6; 2.*x1.*x4.^2; 2.*x1.*x4.*x5; 2.*x1.*x4.*x6; 2.*x1.*x5.^2; 2.*x1.*x5.*x6; 2.*x1.*x6.^2; x2.^3; x2.^2.*x3; x2.^2.*x4; x2.^2.*x5; x2.^2.*x6; x2.*x3.^2; x2.*x3.*x4; x2.*x3.*x5; x2.*x3.*x6; x2.*x4.^2; x2.*x4.*x5; x2.*x4.*x6; x2.*x5.^2; x2.*x5.*x6; x2.*x6.^2; x3.^3; x3.^2.*x4; x3.^2.*x5; x3.^2.*x6; x3.*x4.^2; x3.*x4.*x5; x3.*x4.*x6; x3.*x5.^2; x3.*x5.*x6; x3.*x6.^2; x4.^3; x4.^2.*x5; x4.^2.*x6; x4.*x5.^2; x4.*x5.*x6; x4.*x6.^2; x5.^3; x5.^2.*x6; x5.*x6.^2; x6.^3; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray]';
    diffX2 = [ zeroArray; x1.^3; zeroArray; zeroArray; zeroArray; zeroArray; 2.*x1.^2.*x2; x1.^2.*x3; x1.^2.*x4; x1.^2.*x5; x1.^2.*x6; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; 3.*x1.*x2.^2; 2.*x1.*x2.*x3; 2.*x1.*x2.*x4; 2.*x1.*x2.*x5; 2.*x1.*x2.*x6; x1.*x3.^2; x1.*x3.*x4; x1.*x3.*x5; x1.*x3.*x6; x1.*x4.^2; x1.*x4.*x5; x1.*x4.*x6; x1.*x5.^2; x1.*x5.*x6; x1.*x6.^2; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; 4.*x2.^3; 3.*x2.^2.*x3; 3.*x2.^2.*x4; 3.*x2.^2.*x5; 3.*x2.^2.*x6; 2.*x2.*x3.^2; 2.*x2.*x3.*x4; 2.*x2.*x3.*x5; 2.*x2.*x3.*x6; 2.*x2.*x4.^2; 2.*x2.*x4.*x5; 2.*x2.*x4.*x6; 2.*x2.*x5.^2; 2.*x2.*x5.*x6; 2.*x2.*x6.^2; x3.^3; x3.^2.*x4; x3.^2.*x5; x3.^2.*x6; x3.*x4.^2; x3.*x4.*x5; x3.*x4.*x6; x3.*x5.^2; x3.*x5.*x6; x3.*x6.^2; x4.^3; x4.^2.*x5; x4.^2.*x6; x4.*x5.^2; x4.*x5.*x6; x4.*x6.^2; x5.^3; x5.^2.*x6; x5.*x6.^2; x6.^3; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray]';
    diffX3 = [ zeroArray; zeroArray; x1.^3; zeroArray; zeroArray; zeroArray; zeroArray; x1.^2.*x2; zeroArray; zeroArray; zeroArray; 2.*x1.^2.*x3; x1.^2.*x4; x1.^2.*x5; x1.^2.*x6; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; x1.*x2.^2; zeroArray; zeroArray; zeroArray; 2.*x1.*x2.*x3; x1.*x2.*x4; x1.*x2.*x5; x1.*x2.*x6; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; 3.*x1.*x3.^2; 2.*x1.*x3.*x4; 2.*x1.*x3.*x5; 2.*x1.*x3.*x6; x1.*x4.^2; x1.*x4.*x5; x1.*x4.*x6; x1.*x5.^2; x1.*x5.*x6; x1.*x6.^2; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; x2.^3; zeroArray; zeroArray; zeroArray; 2.*x2.^2.*x3; x2.^2.*x4; x2.^2.*x5; x2.^2.*x6; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; 3.*x2.*x3.^2; 2.*x2.*x3.*x4; 2.*x2.*x3.*x5; 2.*x2.*x3.*x6; x2.*x4.^2; x2.*x4.*x5; x2.*x4.*x6; x2.*x5.^2; x2.*x5.*x6; x2.*x6.^2; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; 4.*x3.^3; 3.*x3.^2.*x4; 3.*x3.^2.*x5; 3.*x3.^2.*x6; 2.*x3.*x4.^2; 2.*x3.*x4.*x5; 2.*x3.*x4.*x6; 2.*x3.*x5.^2; 2.*x3.*x5.*x6; 2.*x3.*x6.^2; x4.^3; x4.^2.*x5; x4.^2.*x6; x4.*x5.^2; x4.*x5.*x6; x4.*x6.^2; x5.^3; x5.^2.*x6; x5.*x6.^2; x6.^3; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray]';
    diffX4 = [ zeroArray; zeroArray; zeroArray; x1.^3; zeroArray; zeroArray; zeroArray; zeroArray; x1.^2.*x2; zeroArray; zeroArray; zeroArray; x1.^2.*x3; zeroArray; zeroArray; 2.*x1.^2.*x4; x1.^2.*x5; x1.^2.*x6; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; x1.*x2.^2; zeroArray; zeroArray; zeroArray; x1.*x2.*x3; zeroArray; zeroArray; 2.*x1.*x2.*x4; x1.*x2.*x5; x1.*x2.*x6; zeroArray; zeroArray; zeroArray; zeroArray; x1.*x3.^2; zeroArray; zeroArray; 2.*x1.*x3.*x4; x1.*x3.*x5; x1.*x3.*x6; zeroArray; zeroArray; zeroArray; 3.*x1.*x4.^2; 2.*x1.*x4.*x5; 2.*x1.*x4.*x6; x1.*x5.^2; x1.*x5.*x6; x1.*x6.^2; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; x2.^3; zeroArray; zeroArray; zeroArray; x2.^2.*x3; zeroArray; zeroArray; 2.*x2.^2.*x4; x2.^2.*x5; x2.^2.*x6; zeroArray; zeroArray; zeroArray; zeroArray; x2.*x3.^2; zeroArray; zeroArray; 2.*x2.*x3.*x4; x2.*x3.*x5; x2.*x3.*x6; zeroArray; zeroArray; zeroArray; 3.*x2.*x4.^2; 2.*x2.*x4.*x5; 2.*x2.*x4.*x6; x2.*x5.^2; x2.*x5.*x6; x2.*x6.^2; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; x3.^3; zeroArray; zeroArray; 2.*x3.^2.*x4; x3.^2.*x5; x3.^2.*x6; zeroArray; zeroArray; zeroArray; 3.*x3.*x4.^2; 2.*x3.*x4.*x5; 2.*x3.*x4.*x6; x3.*x5.^2; x3.*x5.*x6; x3.*x6.^2; zeroArray; zeroArray; zeroArray; zeroArray; 4.*x4.^3; 3.*x4.^2.*x5; 3.*x4.^2.*x6; 2.*x4.*x5.^2; 2.*x4.*x5.*x6; 2.*x4.*x6.^2; x5.^3; x5.^2.*x6; x5.*x6.^2; x6.^3; zeroArray; zeroArray; zeroArray; zeroArray; zeroArray]';
    diffX5 = [ zeroArray; zeroArray; zeroArray; zeroArray; x1.^3; zeroArray; zeroArray; zeroArray; zeroArray; x1.^2.*x2; zeroArray; zeroArray; zeroArray; x1.^2.*x3; zeroArray; zeroArray; x1.^2.*x4; zeroArray; 2.*x1.^2.*x5; x1.^2.*x6; zeroArray; zeroArray; zeroArray; zeroArray; x1.*x2.^2; zeroArray; zeroArray; zeroArray; x1.*x2.*x3; zeroArray; zeroArray; x1.*x2.*x4; zeroArray; 2.*x1.*x2.*x5; x1.*x2.*x6; zeroArray; zeroArray; zeroArray; x1.*x3.^2; zeroArray; zeroArray; x1.*x3.*x4; zeroArray; 2.*x1.*x3.*x5; x1.*x3.*x6; zeroArray; zeroArray; x1.*x4.^2; zeroArray; 2.*x1.*x4.*x5; x1.*x4.*x6; zeroArray; 3.*x1.*x5.^2; 2.*x1.*x5.*x6; x1.*x6.^2; zeroArray; zeroArray; zeroArray; zeroArray; x2.^3; zeroArray; zeroArray; zeroArray; x2.^2.*x3; zeroArray; zeroArray; x2.^2.*x4; zeroArray; 2.*x2.^2.*x5; x2.^2.*x6; zeroArray; zeroArray; zeroArray; x2.*x3.^2; zeroArray; zeroArray; x2.*x3.*x4; zeroArray; 2.*x2.*x3.*x5; x2.*x3.*x6; zeroArray; zeroArray; x2.*x4.^2; zeroArray; 2.*x2.*x4.*x5; x2.*x4.*x6; zeroArray; 3.*x2.*x5.^2; 2.*x2.*x5.*x6; x2.*x6.^2; zeroArray; zeroArray; zeroArray; x3.^3; zeroArray; zeroArray; x3.^2.*x4; zeroArray; 2.*x3.^2.*x5; x3.^2.*x6; zeroArray; zeroArray; x3.*x4.^2; zeroArray; 2.*x3.*x4.*x5; x3.*x4.*x6; zeroArray; 3.*x3.*x5.^2; 2.*x3.*x5.*x6; x3.*x6.^2; zeroArray; zeroArray; x4.^3; zeroArray; 2.*x4.^2.*x5; x4.^2.*x6; zeroArray; 3.*x4.*x5.^2; 2.*x4.*x5.*x6; x4.*x6.^2; zeroArray; 4.*x5.^3; 3.*x5.^2.*x6; 2.*x5.*x6.^2; x6.^3; zeroArray]';
    diffX6 = [ zeroArray; zeroArray; zeroArray; zeroArray; zeroArray; x1.^3; zeroArray; zeroArray; zeroArray; zeroArray; x1.^2.*x2; zeroArray; zeroArray; zeroArray; x1.^2.*x3; zeroArray; zeroArray; x1.^2.*x4; zeroArray; x1.^2.*x5; 2.*x1.^2.*x6; zeroArray; zeroArray; zeroArray; zeroArray; x1.*x2.^2; zeroArray; zeroArray; zeroArray; x1.*x2.*x3; zeroArray; zeroArray; x1.*x2.*x4; zeroArray; x1.*x2.*x5; 2.*x1.*x2.*x6; zeroArray; zeroArray; zeroArray; x1.*x3.^2; zeroArray; zeroArray; x1.*x3.*x4; zeroArray; x1.*x3.*x5; 2.*x1.*x3.*x6; zeroArray; zeroArray; x1.*x4.^2; zeroArray; x1.*x4.*x5; 2.*x1.*x4.*x6; zeroArray; x1.*x5.^2; 2.*x1.*x5.*x6; 3.*x1.*x6.^2; zeroArray; zeroArray; zeroArray; zeroArray; x2.^3; zeroArray; zeroArray; zeroArray; x2.^2.*x3; zeroArray; zeroArray; x2.^2.*x4; zeroArray; x2.^2.*x5; 2.*x2.^2.*x6; zeroArray; zeroArray; zeroArray; x2.*x3.^2; zeroArray; zeroArray; x2.*x3.*x4; zeroArray; x2.*x3.*x5; 2.*x2.*x3.*x6; zeroArray; zeroArray; x2.*x4.^2; zeroArray; x2.*x4.*x5; 2.*x2.*x4.*x6; zeroArray; x2.*x5.^2; 2.*x2.*x5.*x6; 3.*x2.*x6.^2; zeroArray; zeroArray; zeroArray; x3.^3; zeroArray; zeroArray; x3.^2.*x4; zeroArray; x3.^2.*x5; 2.*x3.^2.*x6; zeroArray; zeroArray; x3.*x4.^2; zeroArray; x3.*x4.*x5; 2.*x3.*x4.*x6; zeroArray; x3.*x5.^2; 2.*x3.*x5.*x6; 3.*x3.*x6.^2; zeroArray; zeroArray; x4.^3; zeroArray; x4.^2.*x5; 2.*x4.^2.*x6; zeroArray; x4.*x5.^2; 2.*x4.*x5.*x6; 3.*x4.*x6.^2; zeroArray; x5.^3; 2.*x5.^2.*x6; 3.*x5.*x6.^2; 4.*x6.^3]';

    gradient(:,:,1) = diffX1;
    gradient(:,:,2) = diffX2;
    gradient(:,:,3) = diffX3;
    gradient(:,:,4) = diffX4;
    gradient(:,:,5) = diffX5;
    gradient(:,:,6) = diffX6;

    gradientOneSample = zeros(numMonomials,numDim);
    numSamples = length(x1);
    predictedTwist = zeros(numSamples,numDim);
     for iSample = 1:numSamples
        gradientOneSample(:,:) = gradient(iSample,:,:);
        predictedTwist(iSample,:) = [gradientOneSample' * fittedCoeff]';
     end
end

function labelTitle = getlabel(iAxis)
%     switch iAxis
%         case 1
%             labelTitle = 'Fx';
%         case 2
%             labelTitle = 'Fy';
%         case 3 
%             labelTitle = 'Fz';
%         case 4
%             labelTitle = 'Taux';
%         case 5 
%             labelTitle = 'Tauy';
%         case 6
%             labelTitle = 'Tauz';
%     end
    switch iAxis
        case 1
            labelTitle = '$$f_x$$';
        case 2
            labelTitle = '$$f_y$$';
        case 3 
            labelTitle = '$$f_z$$';
        case 4
            labelTitle = '$$\tau_x$$';
        case 5 
            labelTitle = '$$\tau_y$$';
        case 6
            labelTitle = '$$\tau_z$$';
    end
end





