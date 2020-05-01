%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                     AME 5234 Project
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
% close all
flag_save = 0;
addpath('./Perturbations')
addpath('./../No-Averaged/Matlab codes')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  Define Global Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% my_constants
project_constants
global PC eopdata

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Numerical integration parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

no_yrs = 4;
tf = no_yrs*(365.25*(24*(60*60)));
tspan = [0 tf];
options = odeset('RelTol',1e-12,'AbsTol',1e-12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         Read Ephemeris
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen('eop19620101.txt','r');
eopdata = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 inf]);
fclose(fid);
load DE430Coeff.mat
PC = DE430Coeff;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                Spacecraft Initial Conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% noradID = 'Gurfil';     % Gurfil
% noradID = '00049';      % Echo 1A (LEO)
% noradID = '02253';      % PAGEOS-A (Polar)
% noradID = '02324';      % PasComSat/OV1-8 (LEO)
% noradID = '11659';      % Ariane 1
noradID = '37239';      % Ariane 5 R/B
% noradID = '16657';      % Ariane 3 R/B
% noradID = '19218';      % Ariane 44LP R/B

switch noradID
    case 'Gurfil'
        %
        % Start date and time of simulation
        %
        Mjd_UTC_Epoch = Mjday(2015, 01, 01, 00, 00, 00);
        %
        % Spacecraft properties
        %
        Cd = 2.2;                                   % Drag coefficient of the spacecraft
        AMR = 0.02;                                 % S/m: Area to mass ratio of the spacecraft [m^2/kg]
        delta = 0.5*AMR*Cd;                         % Ballistic coefficient;
        %
        % Orbit properties
        %
        hp0 = 250*1e3;                              % Initial perigee height [m]
        ha0 = 35943*1e3;                            % Initial apogee hegiht [m]
        r_p0 = Re+hp0;                              % Initial perigee radius [m]
        r_a0 = Re+ha0;                              % Initial apogee radius [m]
        %
        a0 = (r_a0+r_p0)/2;                         % Initial semi-major axis [m]
        e0 = 1-r_p0/a0;                             % Initial eccentricity
        i0 = deg2rad(6);                            % Initial inclination [rad]
        argp0 = deg2rad(178);                       % Initial argument of perigee [rad]
        raan0 = deg2rad(60);                        % Initial RAAN [rad]
        M0 = 0;                                     % Initial mean anomaly [rad]
        %
        H0 = sqrt(a0*mu_earth*(1-e0^2));            % Initial angular momentum
        %
        % Atmospheric properties
        %
        we = 7.2921159e-5;                          % Angular velocity of the earth [rad/s]
        wa = we;                                    % Angular velocity of the atmosphere in z-direction [rad/s]
        %
    case '37239' % Ariane 5 R/B'
        %
        % Start date and time of simulation
        %
        Mjd_UTC_Epoch = Mjday(2010,11,28,9,8,13.2351440191269);
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
    case '19218'
%         [1988,6,23,6,47,34.0918216109276,24655.7171757147,0.174811240538349,4.21788744534367,3.22086962168996,0.734115683001826,2.78425661127548]
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       Initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Compute for atmospheric density (perigee altitude must be < 1000 km)
%
[rho_p0,H_p0] = atmosphere_og(hp0);         % (A) Input perigee altitude
% [rho_p0,H_p0] = atmosphere_gurfil(hp0);     % (B) Input perigee altitude in m

%
% Define initial conditions for integrator based on spacecraft initial conditions
%
Hvec0 = H0*[sin(raan0)*sin(i0) -cos(raan0)*sin(i0) cos(i0)]';
evec0 = e0*[(cos(argp0)*cos(raan0)-cos(i0)*sin(argp0)*sin(raan0)) ...
    (cos(argp0)*sin(raan0)+cos(i0)*sin(argp0)*cos(raan0)) ...
    (sin(argp0)*sin(i0))]';
x0 = [Hvec0;evec0];

%
% Define flag for averaging
%
% (1) for singly averaged
% (2) for doubly averaged
% (3) for triply averaged
avg_flag = 1;

%
% Define atmospheric drag model
%
% (1) for Gurfil
% (2) for Ward
drag_model = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Numerically integrate equations of motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Numerically integrate equations of motion
[t_integrator,xx] = ode113(@(t,x) project_function(t,x,delta,wa,rho_p0,r_p0,H_p0,avg_flag,drag_model,Mjd_UTC_Epoch),tspan,x0,options);
%     [t_intW,xx_W] = ode113(@(t,x) project_function(t,x,delta,wa,zhat,Re,rho_p0,r_p0,H_p0,avg_flag,2,Mjd_UTC_Epoch),tspan,x0,options); flag = 0;
%     [t_intG,xx_G] = ode113(@(t,x) project_function(t,x,delta,wa,zhat,Re,rho_p0,r_p0,H_p0,avg_flag,1,Mjd_UTC_Epoch),tspan,x0,options); flag = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Convert results of integration to Keplerian elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Define H, Hhat, e, ehat from results of integration
H = Mag([xx(:,1)' ; xx(:,2)' ; xx(:,3)']);  % Magnitude of angular momentum
Hhat = [xx(:,1) xx(:,2) xx(:,3)]./H';       % Angular momentum unit vector
evec = [xx(:,4) xx(:,5) xx(:,6)];           % Eccentricity vector
e = Mag([xx(:,4)' ; xx(:,5)' ; xx(:,6)']);  % Eccentricity
ehat = [xx(:,4) xx(:,5) xx(:,6)]./e';       % Eccentricity unit vector

%%% Define some constants for convenience
zhatvec = repmat(zhat',length(t_integrator),1)';       % Replicate zhat into 3xn
xhatvec = repmat(xhat',length(t_integrator),1)';       % Replicate xhat into 3xn
c1 = cross(zhatvec,Hhat');                             % zhat x Hhat
c2 = Mag(c1);                                          % |zhat x Hhat|

%%% Define other Keplerian elements
a = H.^2./(mu_earth*(1-e.^2));
inc = rad2deg(acos(dot(zhatvec,Hhat')));
raan = rad2deg(asin(dot(xhatvec,Hhat')./c2));
argp = rad2deg(acos(dot(evec',c1)./(e.*c2)));

rp = a.*(1-e);                              % Radius of perigee
hp = rp-Re;                                 % Perigee altitude
ra = (2*a-rp);                              % Radius of apogee
ha = ra-Re;                                 % Apogee altitude

idx = 1
project_plotting_basic

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                                   Save Mat files
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if drag_model == 1
%     if avg_flag == 1
%         str1 = 'Gurfil_1';
%     elseif avg_flag == 2
%         str1 = 'Gurfil_2';
%     elseif avg_flag == 3
%         str1 = 'Gurfil_3';
%     end
% elseif drag_model == 2
%     if wa == 0
%         str1 = 'Ward_0';
%     else
%         str1 = 'Ward_wa';
%     end
% end
% dataname = sprintf([noradID ' ' str1 ' data.mat']);
% save(fullfile(pwd,'Data',dataname),'realdata');