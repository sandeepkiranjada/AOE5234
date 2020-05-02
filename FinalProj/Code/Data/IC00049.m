%
% Start date and time of simulation
%
Mjd_UTC_Epoch = Mjday(1960,8,12,22,52,32.6369199156761);
% [a,i,W,w,e,M] a  inc   RAAN   argp   ecc   M
%
% Spacecraft properties
%
Cd = 2;                                   % Drag coefficient of the spacecraft
AMR = 11;                               % S/m: Area to mass ratio of the spacecraft [m^2/kg]
delta = 0.5*AMR*Cd;                         % Ballistic coefficient;
%
% Orbit properties
%
a0 = 7984.93475761389*1e3;
i0 = 0.825144415385987;
raan0 = 4.44347154605910;
argp0 = 0.162045275772745;
e0 = 0.00829310414818048;
M0 = -0.136563170407848;
nu0 = truanamoly(M0,e0);
%
% Orbit properties
r_p0 = a0*(1-e0);                              % Radius of perigee
hp0 = r_p0-Re;                                 % Perigee altitude
r_a0 = (2*a0-r_p0);                            % Radius of apogee
ha0 = r_a0-Re;                                 % Apogee altitude
%
H0 = sqrt(a0*mu_earth*(1-e0^2));               % Initial angular momentum
%
% Atmospheric properties
%
we = 7.2921159e-5;                          % Angular velocity of the earth [rad/s]
wa = we;                                    % Angular velocity of the atmosphere in z-direction [rad/s]


