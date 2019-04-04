function visualizedfitted6dgeometry(wrench,fittedWrench,twist,predictedTwsit)
if nargin<=3
    predictedTwsit = [];
end
if size(wrench,2) == 6
    wrench = wrench';
end
for iAxis1 = 1:6
    for iAxis2 = 1:6
        for iAxis3 = 1:6

            if iAxis1 < iAxis2 && iAxis2 < iAxis3 

                plotDim = [iAxis1,iAxis2,iAxis3];
                xlabelTitle = getlabel(iAxis1);
                ylabelTitle = getlabel(iAxis2);
                zlabelTitle = getlabel(iAxis3);
                figHandle = figure;

                wrench2plot = fittedWrench(plotDim,:);
                k = convhulln(wrench2plot');
                axis tight;

                p=patch('Faces',k,'Vertices',wrench2plot');
                set( p, 'FaceColor', 'g','FaceAlpha', 0.8, 'EdgeColor', 'none' );
                view(-10, 20);
                camlight
                hold on
%                 trimesh(k,wrench2plot(1,:),wrench2plot(2,:),wrench2plot(3,:),...
%                     'FaceAlpha', 0.25, 'EdgeColor', 'g')
                hold on
                if norm(twist) == 0
                    plot3( wrench(plotDim(1),:),wrench(plotDim(2),:),wrench(plotDim(3),:)...
                        , 'Marker', 'o', 'MarkerFaceColor', 'r', 'Markersize', 8, 'LineStyle', 'none')
                    hold on
                else
                    uniqueIdx = getidxconvexhull(wrench(plotDim,:));
%                     uniqueIdx = 1:length(wrench(plotDim,:));
                    visualizewrenchtwist(wrench(plotDim,uniqueIdx),twist(plotDim,uniqueIdx),figHandle);
                    hold on
                    if ~isempty(predictedTwsit)
                        visualizewrenchtwist(wrench(plotDim,uniqueIdx),...
                            predictedTwsit(plotDim,uniqueIdx),figHandle,'m');
                    end
                end

                xlabel(xlabelTitle);
                ylabel(ylabelTitle);
                zlabel(zlabelTitle);
                title([xlabelTitle ' ' ylabelTitle ' ' zlabelTitle]);

            else
                continue
            end
        end
    end
end
end

function labelTitle = getlabel(iAxis)
switch iAxis
    case 1
        labelTitle = 'Fx';
    case 2
        labelTitle = 'Fy';
    case 3 
        labelTitle = 'Fz';
    case 4
        labelTitle = 'Taux';
    case 5 
        labelTitle = 'Tauy';
    case 6
        labelTitle = 'Tauz';
end
end

