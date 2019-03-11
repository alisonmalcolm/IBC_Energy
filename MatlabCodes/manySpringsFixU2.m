function [u] = manySpringsFixU2(k, m, dt, nt, u10, un0, nmasses, force_fn, u_boundary)

u = zeros(nt,nmasses);
u(1,1) = u10;
u(nt,1) = u10;

for ii=2:nt-1
    u(ii+1,1) = (k/m)*dt^2*(u(ii,2)-2*u(ii,1))+2*u(ii,1)-u(ii-1,1)+(1/m)*force_fn(ii+1,1)*dt^2;
    u(ii+1,2:nmasses-1) = (k/m)*dt^2*(u(ii,3:nmasses)-2*u(ii,2:nmasses-1)+u(ii,1:nmasses-2))+2*u(ii,2:nmasses-1)-u(ii-1,2:nmasses-1)+(1/m)*force_fn(ii+1,2:nmasses-1)*dt^2;
    u(ii+1,nmasses) = u_boundary(ii+1,1);
end

