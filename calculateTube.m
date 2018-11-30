function [top_limit,bottom_limit] = calculate_tube(rho_0, time_window, A_on, A_off, on)
%CALCULATE_TUBE Summary of this function goes here
%   Detailed explanation goes here
r = size(A_on,1);

for i=1:time_window
    top_limit(i) = on'*(A_on^i) * rho_0;
end
for i=1:time_window
    bottom_limit(i) = on'*(A_off^i) * rho_0;
end

end

