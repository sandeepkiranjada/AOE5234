function  [dxdt] = project_function_ward_withJ2lunisolar(~,x,delta,wa,zhat,rho_p0,rp0,H_p0)

global Re mu_earth mu_sun mu_moon

% disp('----------------------------------------------------------');

%%% Define H, Hhat, e, ehat from results of integration step
H = norm([x(1) x(2) x(3)]);     % Magnitude of angular momentum  
Hhat = [x(1) x(2) x(3)]/H;      % Angular momentum unit vector
e = norm([x(4) x(5) x(6)]);     % Eccentricity
ehat = [x(4) x(5) x(6)]/e;      % Eccentricity unit vector

%%% Define other parameters in terms of state
a = H^2/(mu_earth*(1-e^2));     % Semi-major axis [m]
rp = a*(1-e);                   % Perigee radius [m]
hp = rp-Re;                     % Perigee altitude [m]
inc = acos(dot(zhat,Hhat));     % Inclination [rad]

%%% (A) H_p0 is constant
%%% Compute for atmospheric density (perigee altitude must be < 1000 km)
rho = rho_p0*exp((rp0-a)/(H_p0));
beta = 1/(H_p0);                       % Inverse of scale height [1/m]

%%% Compute for Bessel coefficients
K = beta*a*e;
I0 = besseli(0,K);
I1 = besseli(1,K);
I2 = besseli(2,K);

%%% Compute for terms in cross product needed for differential equations
ehatperp = cross(Hhat,ehat);
c1 = ((1+e^2)*I0-2*e*I1)*dot(ehatperp,zhat)*ehat-0.5*(1-e^2)*(I0-I2)*dot(ehat,zhat)*ehatperp;
c2 = 2*beta*a*(1-e^2);
c3 = 0.5*(1-e^2)*(I0-I2)*dot(ehat,zhat)*ehatperp;

%%% Terms needed for the differential equations
dH = -((delta*rho*H^2)/a)*(I0 + e/(2*beta*a*(1-e^2))*I1 - (2*wa*a^2*cos(inc)/H)*((1+e^2)*I0-2*e*I1));
dHhat = delta*rho*wa*a*cross(c1,Hhat);
de = -2*(delta*rho*H/a)*((1-(2-e^2)/c2)*I1 + (1-1/c2)*e*I0 - ((2*wa*a^2*(1-e^2)*cos(inc))/H)*(I1-e*I0));
dehat = -delta*rho*wa*a*cross(c3,Hhat);

%%% Differential equations for angular momentum (H) and eccentricity vector (e)
dHvec = dH*Hhat + H*dHhat;
devec = de*ehat + e*dehat;

%%% Save differential equations to output variable for function
dxdt = [dHvec';devec'];

end