function [] = savemat(dir,matName,varNameArray,varValue)
%savemat savemat first create directory is not exist, then save the
%variables to mat
% input varArray should be cell
    createdir(dir);
    for iVar = 1:length(varValue)
        varName = varNameArray{iVar};
        saveStruct.(varName) =  varValue{iVar};
        if exist([dir matName],'file')
            save([dir matName],'-struct','saveStruct',varName, '-append');
        else
            save([dir matName],'-struct','saveStruct',varName);
        end
        
    end
end
