%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Initial conditions for PasComSat/OV1-8 (LEO)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Load real data: realdata = [YYYY,MM,DD,hh,mm,ss,a,inc,RAAN,argp,ecc,M,rx,ry,rz,vx,vy,vz]
%
load 02324_data

%
% Define which row of real data to use as initial condition
%
line_no = 10;

%
% Start date and time of simulation
%
Mjd_UTC_Epoch = Mjday(realdata(line_no,1),realdata(line_no,2),...
                realdata(line_no,3),realdata(line_no,4),...
                realdata(line_no,5),realdata(line_no,6));

%           
% Spacecraft properties
%
Cd = 2;                                   % Drag coefficient of the spacecraft
AMR = ((9.1)/2)^2*pi/10.1;                % S/m: Area to mass ratio of the spacecraft [m^2/kg]
                                          % https://nssdc.gsfc.nasa.gov/nmc/spacecraft/display.action?id=1966-063A
                                          % No ref for Diameter
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
r0 = realdata(line_no,13:15)';
v0 = realdata(line_no,16:18)';
r_p0 = a0*(1-e0);                              % Radius of perigee
hp0 = r_p0-Re;                                 % Perigee altitude
r_a0 = (2*a0-r_p0);                            % Radius of apogee
ha0 = r_a0-Re;                                 % Apogee altitude
H0 = sqrt(a0*mu_earth*(1-e0^2));               % Initial angular momentum

%
% Atmospheric properties
%
we = 7.2921159e-5;                          % Angular velocity of the earth [rad/s]
wa = we;                                    % Angular velocity of the atmosphere in z-direction [rad/s]


