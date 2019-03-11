function [u_springs] = springs_from_masses(nt, nmasses, u_masses)

u_springs = zeros(nt, nmasses+1);
u_springs(:,1) = u_masses(:,1);
u_springs(:,2:nmasses) = u_masses(:,2:nmasses)-u_masses(:,1:nmasses-1);
u_springs(:,nmasses+1) = u_masses(:,nmasses);
