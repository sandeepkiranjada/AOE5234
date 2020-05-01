%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                     AME 559 Project
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all
flag_save = 0;

global Re J2 mu_earth mu_sun mu_moon

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Numerical integration parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

no_yrs = 10;
tf = no_yrs*(365*(24*(60*60)));
tspan = [0 tf];
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
% options = odeset('RelTol',1e-12,'AbsTol',1e-12,'Events',@myEvent);

AMRvec = [0.02  0.02  0.02  0.02  0.01  0.005];     % S/m: Area to mass ratio of the spacecraft [m^2/kg]
ha0vec = [35943 20000 10000 35943 35943 35943]*1e3; % Initial apogee height [m]
% AMRvec = [0.05  0.02  0.02  0.02  0.01  0.005];     % S/m: Area to mass ratio of the spacecraft [m^2/kg]
% ha0vec = [25000 20000 10000 35943 35943 35943]*1e3; % Initial apogee height [m]
format longg

for idx = 1%:length(AMRvec)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                        System parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    J2 = 0.0010826267;                          % Geopotential perturbation
    mu_earth = 3.986004418e14;                  % Gravitational parameter of the Earth [m^3/s^2]
    mu_sun = 132712440041.939400e9;             % Gravitational parameter of the Sun [m^3/s^2]
    mu_moon = 4.902800144440050e+12;            % Gravitational parameter of the Moon [m^3/s^2]
    Cd = 2.2;                                   % Drag coefficient of the spacecraft
    AMR = AMRvec(idx);                          % S/m: Area to mass ratio of the spacecraft [m^2/kg]
    delta = 0.5*AMR*Cd;                         % Ballistic coefficient;
    we = 7.2921159e-5;                          % Angular velocity of the earth [rad/s]
    wa = we;0;0.2*we;                           % Angular velocity of the atmosphere in z-direction [rad/s]
    xhat = [1;0;0];                             % Inertial X-direction
    yhat = [0;1;0];                             % Inertial Y-direction
    zhat = [0;0;1];                             % Inertial Z-direction
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                       Initial conditions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Re = 6378*1e3;                              % Radius of the Earth [m]
    
    hp0 = 250*1e3;                              % Initial perigee height [m]
    ha0 = ha0vec(idx);                          % Initial apogee hegiht [m]
    rp0 = Re+hp0;                               % Initial perigee radius [m]
    ra0 = Re+ha0;                               % Initial apogee radius [m]
    
    a0 = (ra0+rp0)/2;                           % Initial semi-major axis [m]
    e0 = 1-rp0/a0;                              % Initial eccentricity
    i0 = deg2rad(6);                            % Initial inclination [rad]
    argp0 = deg2rad(178);                       % Initial argument of perigee [rad]
    raan0 = deg2rad(60);                        % Initial RAAN [rad]
    M0 = 0;                                     % Initial mean anomaly [rad]
    
    H0 = sqrt(a0*mu_earth*(1-e0^2));                  % Initial angular momentum
    
    jd0 = juliandate(datetime([2015 01 01 00 00 00]));
    
    %% Compute for atmospheric density (perigee altitude must be < 1000 km)
	[rho_p0,H_p0] = atmosphere_og(hp0);          % Input perigee altitude
    % [rho_p0,H_p0] = atmosphere_gurfil(hp0);          % Input perigee altitude
    % [rho_0,H_p0,h_p0,rho_p0] = atmosphere(hp0*1e-3);  % Input perigee altitude in km
    % r_p0 = (h_p0*1e3)+Re;                             % Initial perigee radius [m]
    % H_p0 = H_p0*1e3;                                  % Initial scale height [m]
    
    %%% Define initial conditions
    Hvec0 = H0*[sin(raan0)*sin(i0) -cos(raan0)*sin(i0) cos(i0)]';
    evec0 = e0*[(cos(argp0)*cos(raan0)-cos(i0)*sin(argp0)*sin(raan0)) ...
                (cos(argp0)*sin(raan0)+cos(i0)*sin(argp0)*cos(raan0)) ...
                (sin(argp0)*sin(i0))]';
    x0 = [Hvec0;evec0];
    
%     check = [rp0;a0;e0;rho_p0;H_p0];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                           Numerically integrate equations of motion
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Numerically integrate equations of motion
    [t_integrator,xx] = ode113(@(t,x) project_function_ward(t,x,delta,wa,zhat,rho_p0,rp0,H_p0,1,1,jd0),tspan,x0,options); flag = 0;
    % [t_integrator,xx] = ode113(@(t,x) project_function_gurfil(t,x,mu_earth,delta,wa,zhat,Re,rp0,rho_p0,H_p0),tspan,x0,options); flag = 1;
    
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                         Plot results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    project_plotting
    
end

% %{

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         Define Legend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q = get(0,'children');
for f = 1:length(q)
    figure(f);
    grid on
    
    %%% Modify properties of figures for aesthetics  
    set(findall(gcf,'Type','line'),'LineWidth',1)
    set(gca,'FontSize',12);
    if f <= 8
        legend('h_a(0) = 35,943 km','h_a(0) = 20,000 km','h_a(0) = 10,000 km','Location','Best');
    else
        legend('AMR = 0.02 m^2/kg','AMR = 0.01 m^2/kg','AMR = 0.005 m^2/kg','Location','Best')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                        Save Images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_save == 1
    for f = 1:length(q)
        figure(f);
        if flag == 0
            print(q(f),['Figure Ward ' num2str(length(q)+1-f,'%01.0f')],'-dpng','-r300');
        elseif flag == 1
            print(q(f),['Figure Gurfil' num2str(length(q)+1-f,'%01.0f')],'-dpng','-r300');
        end
    end
    
else
end
% %}