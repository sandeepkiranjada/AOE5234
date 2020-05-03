%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                     AME 5234 Project
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all
flag_save = 0;
addpath('./Formulation');
addpath('./Perturbations v1')
addpath('./Data')
addpath('./Post')
addpath('./../No-Averaged/Matlab codes')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  Define Global Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% my_constants
project_constants
global eopdata

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
% noradID = '16657';      % Ariane 3 R/B
% noradID = '19218';      % Ariane 44LP R/B
noradID = '37239';      % Ariane 5 R/B

switch noradID
    case 'Gurfil'
        ICGURFIL
    case '00049'  % Echo 1A (LEO)
        IC00049
    case '02253'  % PAGEOS-A (Polar)
        IC02253
    case '02324'  % PasComSat/OV1-8 (LEO)
        IC02324   
    case '11659'  % Ariane 1
        IC11659  
%     case '16657'  % Ariane 3 R/B
%         IC16657   
    case '19218'  % Ariane 44LP R/B
        IC19218
    case '37239'  % Ariane 5 R/B
        IC37239    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       Initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Compute for atmospheric density (perigee altitude must be < 1000 km)
%
[rho_p0,H_p0] = atmosphere_og(hp0);         % (A) Input perigee altitude
% [~,H_p0] = atmosphere_gurfil(hp0);     % (B) Input perigee altitude in m
% [~,rho_p0_all] = atmosnrlmsise00(hp0,0,0,2010,1,0,'None');     % (B) Input perigee altitude in m
% rho_p0 = rho_p0_all(6);
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
% drag_model = 2;

%
% Define which perturbations to consider
%   atm drag    lunisolar   J2
% [ 1/0        1/0        1/0 ]
pert_fac = [1 1 1]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Numerically integrate equations of motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Numerically integrate equations of motion
%

% Gurfil 
drag_model = 1;
[t_integrator,xx] = ode113(@(t,x) project_function(t,x,delta,wa,rho_p0,r_p0,H_p0,avg_flag,drag_model,...
                                                        Mjd_UTC_Epoch,pert_fac,PC),tspan,x0,options);
gurfil = milankovitch2coe(t_integrator,xx);                                                    

% Ward wa = we
drag_model = 2;
[t_integrator,xx] = ode113(@(t,x) project_function(t,x,delta,wa,rho_p0,r_p0,H_p0,avg_flag,drag_model,...
                                                        Mjd_UTC_Epoch,pert_fac,PC),tspan,x0,options);
wardw = milankovitch2coe(t_integrator,xx);   

% Ward wa = 0
wa = 0;
[t_integrator,xx] = ode113(@(t,x) project_function(t,x,delta,wa,rho_p0,r_p0,H_p0,avg_flag,drag_model,...
                                                        Mjd_UTC_Epoch,pert_fac,PC),tspan,x0,options);
ward0 = milankovitch2coe(t_integrator,xx); 

project_plotting_GurfWardReal
