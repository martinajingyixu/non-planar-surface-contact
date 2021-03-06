function [integratedFriction] = mergeintegraldata(pathProcessingSurface,idxStart,idxEnd,integratedFriction)
% Merge friction Array from different files
    filename = [pathProcessingSurface '/friction_sphere_Twist_center_' ...
        num2str(idxStart) '_' num2str(idxEnd) '.mat'];
    storedFriction = load(filename);
    
    integratedFriction.fxArray(idxStart:idxEnd) = ...
        storedFriction.fxArray(idxStart:idxEnd);
    
    integratedFriction.fyArray(idxStart:idxEnd)= ...
    storedFriction.fyArray(idxStart:idxEnd);

    integratedFriction.fzArray(idxStart:idxEnd) = ...
    storedFriction.fzArray(idxStart:idxEnd);

    integratedFriction.tauxArray(idxStart:idxEnd) = ...
    storedFriction.tauxArray(idxStart:idxEnd);

    integratedFriction.tauyArray(idxStart:idxEnd)= ...
    storedFriction.tauyArray(idxStart:idxEnd);

    integratedFriction.tauzArray(idxStart:idxEnd) = ...
    storedFriction.tauzArray(idxStart:idxEnd);

%%
%     integratedFriction.fxMinusArray(idxStart:idxEnd) = ...
%         storedFriction.fxMinusArray(idxStart:idxEnd);
%     
%     integratedFriction.fyMinusArray(idxStart:idxEnd)= ...
%     storedFriction.fyMinusArray(idxStart:idxEnd);
% 
%     integratedFriction.fzMinusArray(idxStart:idxEnd) = ...
%     storedFriction.fzMinusArray(idxStart:idxEnd);
% 
%     integratedFriction.tauxMinusArray(idxStart:idxEnd) = ...
%     storedFriction.tauxMinusArray(idxStart:idxEnd);
% 
%     integratedFriction.tauyMinusArray(idxStart:idxEnd)= ...
%     storedFriction.tauyMinusArray(idxStart:idxEnd);
% 
%     integratedFriction.tauzMinusArray(idxStart:idxEnd) = ...
%     storedFriction.tauzMinusArray(idxStart:idxEnd);
%%
%     integratedFriction.fxIntFailed (idxStart:idxEnd)= ...
%         storedFriction.fxIntFailed(idxStart:idxEnd);
%     
%     integratedFriction.fyIntFailed(idxStart:idxEnd) = ...
%     storedFriction.fyIntFailed (idxStart:idxEnd);
% 
%     integratedFriction.fzIntFailed(idxStart:idxEnd)  = ...
%     storedFriction.fzIntFailed (idxStart:idxEnd);
% 
%     integratedFriction.tauxIntFailed(idxStart:idxEnd)  = ...
%     storedFriction.tauxIntFailed (idxStart:idxEnd);
% 
%     integratedFriction.tauyIntFailed(idxStart:idxEnd) = ...
%     storedFriction.tauyIntFailed (idxStart:idxEnd);
% 
%     integratedFriction.tauzIntFailed(idxStart:idxEnd)  = ...
%     storedFriction.tauzIntFailed(idxStart:idxEnd);

end

