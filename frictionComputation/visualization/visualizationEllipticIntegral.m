%%



% idxm = 1:length(m);
plot3(phi_t0(idx),m(idx),fittedModel.ellipticFt0(idx),'b.');
hold on
plot3(phi_t0(idx),m(idx),tan(phi_t0(idx)),'b.');
hold on

plot3(phi_t1(idx),m(idx),ellipticResults.ellipticFt0(idx),'r.');

plot3(phi_t0(idx),m(idx),ellipticResults.ellipticEt0(idx),'r.');


%%
idxphi = find(phi_t0<0.2*pi & phi_t0 > 0.01*pi);
idxm = find(m<-2);
idxIe0 = find(abs(ellipticResults.ellipticFt0)<20);

idx = intersect(idxm,idxIe0);
idx = intersect(idx,idxIe0A);
idx = intersect(idx,idxphi);

plot3(phi_t0(idx),m(idx),tan(phi_t0(idx)),'b.');
hold on
% plot3(phi_t0(idx),m(idx),ellipticResults.ellipticFt0(idx),'r.');
plot3(phi_t0(idx),m(idx),ellipticResults.IMe10EllipticPit0(idx),'r.');




%%
idxm = find(abs(m)<10);
idxt0 = find(abs(ellipticResults.ellipticFt0)<10);
idx = intersect(idxm,idxt0);

phit0real = real(phi_t0);
mreal = real(m);
% reale0 = real(ellipticResults.ellipticFt0);
reale0 = real(ellipticResults.ellipticFt0);

f = fit([phit0real(idx),mreal(idx)],reale0(idx),'poly44');
plot(f,[phit0real(idx),mreal(idx)],reale0(idx))
%%
idxm = find(abs(m)<10);
idxt0 = find(abs(ellipticResults.IMe10EllipticPit1)<10);
idx = intersect(idxm,idxt0);

phit0real = real(phi_t0);
mreal = real(m);
% reale0 = real(ellipticResults.ellipticFt0);
reale0 = real(ellipticResults.IMe10EllipticPit1);

f = fit([phit0real(idx),mreal(idx)],reale0(idx),'poly44');
plot(f,[phit0real(idx),mreal(idx)],reale0(idx))




%%
idxm = find(abs(m)<20);
idxIe0 = find(abs(tauxAppr(:))<20);
idxIe0A = find(abs(taux(:))<20);
idx = intersect(idxm,idxIe0);
idx = intersect(idx,idxIe0A);
% plot3(phi_t0(idx),m(idx),IMinuse9Array(idx),'b.');
% plot3(phi_t0(idx),m(idx),termIMinuse2(idx,column),'b.');
hold on
plot3(phi_t0(idx),m(idx),tauxAppr(idx),'r.');
hold on
plot3(phi_t0(idx),m(idx),taux(idx),'g.');
%%
idxm = find(abs(m)<20);
idxIe0 = find(abs(tauzAppr)<20);
idxIe0A = find(abs(tauz)<20);
idx = intersect(idxm,idxIe0);
idx = intersect(idx,idxIe0A);
% plot3(phi_t0(idx),m(idx),IMinuse9Array(idx),'b.');
% plot3(phi_t0(idx),m(idx),termIMinuse2(idx,column),'b.');
hold on
plot3(phi_t0(idx),m(idx),tauzAppr(idx),'r.');
hold on
axis equal
plot3(phi_t0(idx),m(idx),tauz(idx),'g.');
axis equal


%%
idxm = find(abs(m)<20);
idxIe0 = find(abs(Ie0Array)<10);
idxIe0A = find(abs(Ie0approcArray)<10);
idx = intersect(idxm,idxIe0);
idx = intersect(idx,idxIe0A);
figure
plot3(phi_t0(idx),m(idx),Ie0Array(idx),'b.');
axis equal

figure
plot3(phi_t0(idx),m(idx),Ie0approcArray(idx),'r.');

figure
plot3(phi_t0(idx),m(idx),Ie0Array(idx),'b.');
hold on
plot3(phi_t0(idx2),m(idx2),IMinuse9Array(idx2),'r.');
hold on
plot3(phi_t0(idxe10),m(idxe10),IMinuse10Array(idxe10),'g.');
hold on
plot3(phi_t0(idxe2),m(idxe2),IMinuse2Array(idxe2),'b.');

I = [Ie0Array,IMinuse9Array,IMinuse10Array,IMinuse2Array];
%%
idxm = find(abs(m)<20);
idxIe0 = find(abs(IMinuse9Array)<10);
idxIe0A = find(abs(Ie0Array)<10);
idx = intersect(idxm,idxIe0);
idx = intersect(idx,idxIe0A);
figure
plot3(phi_t0(idx),m(idx),IMinuse9Array(idx),'b.');
hold on
plot3(phi_t0(idx),m(idx),0.4*Ie0Array(idx),'r.');
%%
idxm = find(abs(m)<20);
idxIe0 = find(abs(IMinuse10Array)<10);
idxIe0A = find(abs(IMe10approcArray)<10);
idx = intersect(idxm,idxIe0);
idx = intersect(idx,idxIe0A);
plot3(phi_t0(idx),m(idx),IMinuse10Array(idx),'b.');
hold on
plot3(phi_t0(idx),m(idx),IMe10approcArray(idx),'r.');
%%
idxm = find(abs(m)<20);
idxIe0 = find(abs(IMinuse2Array)<10);
idxIe0A = find(abs(IMe2approcArray)<10);
idxe2 = intersect(idxm,idxIe0);
idxe2 = intersect(idxe2,idxIe0A);
plot3(phi_t0(idx),m(idx),IMinuse2Array(idx),'b.');
hold on
plot3(phi_t0(idx),m(idx),IMe2approcArray(idx),'r.');
%%
scatter3(phi_t0(idx),m(idx),ne10(idx),ellipticResults.IMe10EllipticPit0(idx)+10,ellipticResults.IMe10EllipticPit0(idx)+10)
hold on
colorbar
%%
ground_truth_t0 = ellipticResults.ellipticFt0(idx);
ground_truth_t1 = ellipticResults.ellipticFt1(idx);
ground_truth_int = ground_truth_t1 - ground_truth_t0;
approc_t0 = phi_t0(idx);
approc_t1 = phi_t1(idx);
diff_int = approc_t1 - approc_t1;

idx1 = find(abs(Ie0approc)<3);
idx2 = find(abs(Ie0Array)<3);
idx = intersect(idx1,idx2);
plot(1:numel(real(Ie0approc(idx))),real(Ie0approc(idx)),'b');
hold on
plot(1:numel(real(Ie0Array(idx))),real(Ie0Array(idx)),'r');