close all; clear all; clc;
cd D:\ETH\Matlab\SpringMassEnergy\share_110319
%path(path,'..\scripts\')

N = 150; %120;
K = 1.2*N;
M = 1.2/N;
c = sqrt(K/M)

tmax = pi;
dt = 0.0005; %0.001;
t = 0:dt:tmax;
nt = length(t);

nmasses = N;

fc = 5; force_fn = zeros(nt,nmasses);
force_fn(:,1) = ricker_time(t,fc,1/fc)';
figure; plot(t,force_fn);

u_masses = manySprings2(K, M, dt, nt, 0, 0, nmasses, force_fn);
u_springs = springs_from_masses(nt, nmasses, u_masses);

figure; imagesc(u_masses); colorbar;
figure; imagesc(1:N, t, u_masses); colorbar;

KE = KineticE(M, u_masses, dt);
PE = PotentialE(K, u_springs);
PEavgx = (PE(:,1:nmasses)+PE(:,2:nmasses+1))/2;
PEavg = (PEavgx(1:end-1,:)+PEavgx(2:end,:))/2;

figure; 
sh1=subplot(1,3,1); imagesc(KE); ca=caxis; colorbar;
sh2=subplot(1,3,2); imagesc(PEavg); caxis(ca); colorbar;
sh3=subplot(1,3,3); imagesc(KE+PEavg); caxis(ca); colorbar;

trcno = 60;
figure; hold on; plot(KE(:,trcno),'b'); plot(PEavg(:,trcno),'r'); plot(KE(:,trcno)+PEavg(:,trcno),'g');

PEtrc = PE(:,trcno+1);
KEtrc = KE(:,trcno);
PEtrc_before = PE(:,trcno);
PEtrc_avg = PEavg(:,trcno);

%-------------------------------------------------------------------------
% Now truncate the medium with an IBC
%-------------------------------------------------------------------------
nmasses2 = 60;
nmute = 2500;

force2 = zeros(nt,nmasses2);
force2(:,1) = force_fn(:,1);
u_boundary = zeros(nt,1);
f_boundary = zeros(nt,1);
u_masses2_u = manySpringsFixU2(K, M, dt, nt, 0, 0, nmasses2, force2, u_boundary);
u_masses2_f = manySpringsFixF2(K, M, dt, nt, 0, 0, nmasses2, force2, f_boundary);

force2 = zeros(nt,nmasses2);
u_boundary = u_masses(:,nmasses2);
u_boundary(nmute:end,1) = 0;
f_boundary = K*u_springs(:,nmasses2+1);
f_boundary(nmute:end,1) = 0;
u_masses2_IBCu = manySpringsFixU2(K, M, dt, nt, 0, 0, nmasses2, force2, u_boundary);
u_masses2_IBCf = manySpringsFixF2(K, M, dt, nt, 0, 0, nmasses2, force2, f_boundary);

force2 = zeros(nt,nmasses2);
force2(:,1) = force_fn(:,1);
u_masses2_allu = manySpringsFixU2(K, M, dt, nt, 0, 0, nmasses2, force2, u_boundary);
u_masses2_allf = manySpringsFixF2(K, M, dt, nt, 0, 0, nmasses2, force2, f_boundary);

figure; sc = 1e4; %000; %e4;
sh1=subplot(1,3,1); imagesc(u_masses2_u); title('force'); ca = caxis;
sh2=subplot(1,3,2); imagesc(u_masses2_IBCu); title('IBC'); caxis(ca);
sh3=subplot(1,3,3); imagesc(u_masses2_allu); title(['all (x' num2str(sc) ')']); caxis(ca/sc);

figure; sc = 1e4; %000; %e4;
sh1=subplot(1,3,1); imagesc(u_masses2_f); title('force'); ca = caxis;
sh2=subplot(1,3,2); imagesc(u_masses2_IBCf); title('IBC'); caxis(ca);
sh3=subplot(1,3,3); imagesc(u_masses2_allf); title(['all (x' num2str(sc) ')']); caxis(ca/sc);

v_boundary = (u_boundary(2:end)-u_boundary(1:end-1))/dt;
w_force = v_boundary.*f_boundary(1:end-1);
figure; hold on; plot(c*(PEtrc(1:end-1)+KEtrc),'b'); plot(-w_force,'r--'); 

% verify last mass in truncated medium moves exactly like in extended medium...
u_boundary_trunc = u_masses2_allf(:,nmasses2);
figure; hold on; plot(u_boundary,'b'); plot(u_boundary_trunc,'r--');

en_x = sum(KE(1000,:))
en_t = sum(c*dt*KE(1:1500,50))
