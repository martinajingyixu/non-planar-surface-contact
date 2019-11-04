function displayfittingerror(fittedRes)
% its possible when fitting failes then the fittedRes = []
try
    disp(['Wrench error (4th, ellip): ' num2str(fittedRes.meanWrenchError4th)...
        ', ' num2str(fittedRes.meanWrenchErrorEllip) '. ',...
    'Twist error (4th, ellip): ' num2str(fittedRes.meanTwistError4th)...
    ', ' num2str(fittedRes.meanTwistErrorEllip)])
catch
end

end

