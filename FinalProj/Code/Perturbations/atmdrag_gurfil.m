function [ dxdt ] = atmdrag_gurfil(x,delta,rho_p0,r_p0,H_p0)

%ATMDRAG_GURFIL Summary of this function goes here
%   Detailed explanation goes here

global Re mu_earth 

%
% Define H, Hhat, e, ehat from results of integration step
%
H = norm([x(1) x(2) x(3)]);     % Magnitude of angular momentum  
e = norm([x(4) x(5) x(6)]);     % Eccentricity

%
% Define other parameters in terms of state
%
a = H^2/(mu_earth*(1-e^2));     % Keplerian semi-major axis
rp = a*(1-e);                   % Radius of periapsis [m]
hp = rp-Re;                     % Perigee height [m]

%
% (A) The Gurfil Way
% Compute for atmospheric density (perigee altitude must be < 1000 km)
%
rho = rho_p0*exp((r_p0-rp)/(H_p0));
H_p = H_p0;

%
% Terms needed for the differential equations
%
Hvec = x(1:3);
evec = x(4:6);
B = 2*delta;
z = a*e/(H_p);
K1 = (1+3*e^2)/(8*z*(1-e^2));
K2 = (3*e^2-4*e-3)/(8*z*(1-e^2));

%
% Differential equations for angular momentum (H) and eccentricity vector (e)
%
dHvec = -0.5*B*sqrt(mu_earth*(1-e^2)/(2*a*pi*z))*rho*(1+K1)*Hvec;
devec = -B*((1+e)/(a*sqrt(2*pi*z)))*rho*(1+K2)*H*evec./e;

%
% Save differential equations to output variable for function
%
dxdt = [dHvec;devec];

end

