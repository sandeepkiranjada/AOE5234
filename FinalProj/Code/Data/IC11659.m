[0.700223341795105,0.0712517649170450,2765629.24680352,6407067.76520317,542.102080371019,-7488.04530340572,5586.78254628662,2958.01979559507]

%
% Start date and time of simulation
%
Mjd_UTC_Epoch = Mjday(1980,6,8,20,41,48.3446457982063);
% [a,i,W,w,e,M] a  inc   RAAN   argp   ecc   M
%
% Spacecraft properties
%
Cd = 2.2;                                   % Drag coefficient of the spacecraft
AMR = 0.0153;0.02;                          % S/m: Area to mass ratio of the spacecraft [m^2/kg]
delta = 0.5*AMR*Cd;                         % Ballistic coefficient;
%
% Orbit properties
%
a0 = 21902.5901287491*1e3;
i0 = 0.314506803384171;
raan0 = 1.16298457218633;
argp0 = 5.74073951853275;
e0 = 0.700223341795105;
M0 = 0.0712517649170450;
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