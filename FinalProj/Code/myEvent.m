function [value, isterminal, direction] = myEvent(~, x)

global mu Re

%%% Define e from results of integration step
H = norm([x(1) x(2) x(3)]);     % Magnitude of angular momentum  
e = norm([x(4) x(5) x(6)]);     % Eccentricity

%%% Define other parameters in terms of state
a = H^2/(mu*(1-e^2));           % Semi-major axis [m]
rp = a*(1-e);                   % Perigee radius [m]
hp = rp-Re;                     % Perigee altitude [m]

value      = hp <= 130;
isterminal = 1;   % Stop the integration
direction  = 0;

end