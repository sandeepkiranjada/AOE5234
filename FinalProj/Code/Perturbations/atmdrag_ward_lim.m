function  [dxdt] = atmdrag_ward_lim(x,delta,wa,rho_p0,r_p0,H_p0)

%ATMDRAG_WARD Summary of this function goes here
%   Detailed explanation goes here

global Re mu_earth zhat

%
% Define H, Hhat, e, ehat from results of integration step
%
H_vec = [x(1); x(2); x(3)];     % Angular momentum vector
e_vec = [x(4); x(5); x(6)];     % Eccentricity vector
%
H = norm(H_vec);                % Magnitude of angular momentum  
e = norm(e_vec);                % Eccentricity
%
Hhat = H_vec/H;                 % Angular momentum unit vector
ehat = e_vec/e;                 % Eccentricity unit vector

%    
% Define other parameters in terms of the states
%
a = H^2/(mu_earth*(1-e^2));     % Semi-major axis [m]
rp = a*(1-e);                   % Perigee radius [m]
hp = rp-Re;                     % Perigee altitude [m]
inc = acos(dot(zhat,Hhat));     % Inclination [rad]

% %
% % (A) Simply use atmosphere function from Dr. Rosengren
% % Compute for atmospheric density (perigee altitude must be < 1000 km)
% %
% [rho,H_p] = atmosphere_og(hp);          % Input perigee altitude in km
% beta = 1/(H_p*1e3);                     % Inverse of scale height [1/m]

%
% (B) H_p0 is constant
% Compute for atmospheric density (perigee altitude must be < 1000 km)
%
rho = rho_p0*exp((r_p0-rp)/(H_p0));
% [~,rho_p_all] = atmosnrlmsise00(hp,0,0,2010,1,0,'None');     % (B) Input perigee altitude in m
% rho_p = rho_p_all(6);
% rho = rho_p *  exp((-a*e)/(H_p0));
% rho = rho_p0*exp((r_p0-a)/(H_p0)) ;% *  exp((-a*e)/(H_p0));
beta = 1/(H_p0);                        % Inverse of scale height [1/m]

% %
% % (C) H_p0 changes
% % Compute for atmospheric density (perigee altitude must be < 1000 km)
% %
% [~,H_p] = atmosphere_og(hp);            % Input perigee altitude in km and compute new scale height
% rho = rho_p0*exp((r_p0-a)/(H_p));        %
% beta = 1/(H_p);                         % Inverse of scale height [1/m]

%
% Compute for Bessel coefficients
%
K = beta*a*e

exp_I_nu = @(alpha,z) (1 - (alpha-1)/8/z + (alpha-1)*(alpha-3)/(2*(8*z)^2) - (alpha-1)*(alpha-9)*(alpha-25)/(6*(8*z)^3)) / sqrt(2*pi*z);

if K<300
    I0 = besseli(0,K)*exp(-K);
    I1 = besseli(1,K)*exp(-K);
    I2 = besseli(2,K)*exp(-K);
else
    I0 = exp_I_nu(0,K);
    I1 = exp_I_nu(4,K);
    I2 = exp_I_nu(4*4,K);
end

%
% Compute for terms in cross product needed for differential equations
%
ehatperp = cross(Hhat,ehat);
c1 = ((1+e^2)*I0-2*e*I1)*dot(ehatperp,zhat)*ehat-0.5*(1-e^2)*(I0-I2)*dot(ehat,zhat)*ehatperp;
c2 = 2*beta*a*(1-e^2);
c3 = 0.5*(1-e^2)*(I0-I2)*dot(ehat,zhat)*ehatperp;

%
% Terms needed for the differential equations
%
dH = -((delta*rho*H^2)/a)*(I0 + e/(2*beta*a*(1-e^2))*I1 - (2*wa*a^2*cos(inc)/H)*((1+e^2)*I0-2*e*I1));
dHhat = delta*rho*wa*a*cross(c1,Hhat);
de = -2*(delta*rho*H/a)*((1-(2-e^2)/c2)*I1 + (1-1/c2)*e*I0 - ((2*wa*a^2*(1-e^2)*cos(inc))/H)*(I1-e*I0));
dehat = -delta*rho*wa*a*cross(c3,Hhat);

%
% Differential equations for angular momentum (H) and eccentricity vector (e)
%
dHvec = dH*Hhat + H*dHhat;
devec = de*ehat + e*dehat;

%
% Save differential equations to output variable for function
%
dxdt = [dHvec;devec];

end