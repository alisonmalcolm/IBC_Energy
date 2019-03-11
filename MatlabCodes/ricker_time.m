function [result] = ricker_time(ts,fc,reft)

  result=(1-2*pi^2*fc^2*(ts-reft).^2).*exp(-(pi*fc*(ts-reft)).^2);

return
