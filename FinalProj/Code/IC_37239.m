%
% Start date and time of simulation
%
Mjd_UTC_Epoch = Mjday(2010,11,28,9,8,13.2351440191269);
% [a,i,W,w,e,M]
%
% Spacecraft properties
%
Cd = 2.2;                                   % Drag coefficient of the spacecraft
AMR = 0.0153;0.02;                          % S/m: Area to mass ratio of the spacecraft [m^2/kg]
delta = 0.5*AMR*Cd;                         % Ballistic coefficient;
%
% Orbit properties
%
a0 = 24361.7482201275*1e3;
i0 = 0.0302181724765577;
raan0 = 3.10964407613065;
argp0 = 3.02952610730407;
e0 = 0.728072678688366;
M0 = -2.64967023526191;
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