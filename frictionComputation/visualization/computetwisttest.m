function [twist,dotProdArray,twistCrossOmega] = computetwisttest(com,velProjArray,velArray,corDir,corLoc)
    numSamples = size(velProjArray,1);
    twist = zeros(numSamples,6);
    twistCrossOmega = zeros(numSamples,6);
    dotProdArray = zeros(numSamples,1);
    for iSample = 1:numSamples      
        %%
        omega = corDir(iSample,:)'; % angular velocity
        vnp = velProjArray(iSample,:)';
        vn = velArray(iSample,:)';
        vnn = vn - vnp;
        dotProdArray(iSample,1) = dot(vnn,vnp);
        omegaNewCross = cross(-vnn,vnp)./(norm(vnn)^2);
%         ra = rc + cross(omega,vc)./((norm(omega))^2);
%         omegaNew = omega ./norm(raNew).*norm(ra);
%         h = dot(omega,vc)./((norm(omega))^2);
%         omegaNew = h*omega;
        r = com - corLoc(iSample,:)';
        omegaNew = cross(r,vnp) ./(norm(r)^2);
        twist(iSample,:) = [vnp;omegaNew]';
        twistCrossOmega(iSample,:) = [vnp;omegaNewCross]';

    end
    
end

%%
% function [twist] = computetwisttest(com,velArrayProjected,velArray,corDir)
%     numSamples = size(velArray,1);
%     twist = zeros(6,numSamples);
%     
%     rc = com;
%     for iSample = 1:numSamples
% 
%         %%
%         omega = corDir(iSample,:)'; % angular velocity
%         vc = velArray(iSample,:)';
%         vcp = velArrayProjected(iSample,:)';
%         vo = vcp - vc;
%         assert(abs(dot(omega,vc))<1e-4);
% %         ra = rc + cross(omega,vc)./((norm(omega))^2);
% 
% %         h = dot(omega,vc)./((norm(omega))^2);
% %         twist(:,iSample) = [cross(ra,omega);h*omega];
%         twist(:,iSample) = [vc;vo];
% 
%     end
%     
% end
% 
