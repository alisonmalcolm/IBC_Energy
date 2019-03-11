function [KE] = KineticE(m, u, dt)

KE = 0.5*m*((u(2:end,:)-u(1:end-1,:))/dt).^2;