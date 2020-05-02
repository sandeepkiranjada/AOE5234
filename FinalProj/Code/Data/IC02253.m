load 02253_data

line_no = 1; % This is the one that broke up
%
% Start date and time of simulation
%
Mjd_UTC_Epoch = Mjday(realdata(line_no,1),realdata(line_no,2),...
                realdata(line_no,3),realdata(line_no,4),...
                realdata(line_no,5),realdata(line_no,6));
% [a,i,W,w,e,M] a  inc   RAAN   argp   ecc   M
%
% Spacecraft properties
%
Cd = 2;                                   % Drag coefficient of the spacecraft
AMR = ((30.48)/2)^2*pi/56.7;              % S/m: Area to mass ratio of the spacecraft [m^2/kg]
                                          % https://nssdc.gsfc.nasa.gov/nmc/spacecraft/display.action?id=1966-056A
delta = 0.5*AMR*Cd;                       % Ballistic coefficient;
%
% Orbit properties
%
a0 = realdata(line_no,7)*1e3;
i0 = realdata(line_no,8);
raan0 = realdata(line_no,9);
argp0 = realdata(line_no,10);
e0 = realdata(line_no,11);
M0 = realdata(line_no,12);
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

