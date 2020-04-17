function  [dxdt] = project_function_gurfil(~,x,mu,delta,~,~,Re,rho_p0,r_p0)

% disp('----------------------------------------------------------');

%%% Define H, Hhat, e, ehat from results of integration step
H = norm([x(1) x(2) x(3)]);     % Magnitude of angular momentum  
e = norm([x(4) x(5) x(6)]);     % Eccentricity

%%% Define other parameters in terms of state
a = H^2/(mu*(1-e^2));           % Keplerian semi-major axis
rp = a*(1-e);                   % Radius of periapsis [m]
hp = rp-Re;                     % Perigee height [m]
hp0 = r_p0-Re;
H
%%% Compute for atmospheric density (perigee altitude must be < 1000 km)
[~,H_p] = atmosphere_gurfil(hp*1e-3)      % Input perigee altitude in km
rho = rho_p0*exp((r_p0-rp)/(H_p*1e3));
% %%% Compute for atmospheric density (perigee altitude must be < 1000 km)
% rho = rho_p0*exp((r_p0-rp)/Hp);

%%% Terms needed for the differential equations
Hvec = x(1:3);
evec = x(4:6);
B = 2*delta;
z = a*e/(H_p*1e3);
K1 = (1+3*e^2)/(8*z*(1-e^2));
K2 = (3*e^2-4*e-3)/(8*z*(1-e^2));

%%% Differential equations for angular momentum (H) and eccentricity vector (e)
dHvec = -0.5*B*sqrt(mu*(1-e^2)/(2*a*pi*z))*rho*(1+K1)*Hvec;
devec = -B*((1+e)/(a*sqrt(2*pi*z)))*rho*(1+K2)*H*evec;

%%% Save differential equations to output variable for function
dxdt = [dHvec;devec];

end