%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                     AME 5234 Project
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear
close all
flag_save_matfile = 0;
flag_save_figs = 0;
addpath('./Initial Conditions')
addpath('./Formulation');
addpath('./Perturbations')
addpath('./Data')
addpath('./Post')
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Set up flags for which simulations to run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flag_gurfil = 1;    % Gurfil
flag_wardw  = 1;    % Ward wa = we
flag_ward0  = 1;    % Ward wa = 0
flag_naw    = 0;    % Non-averaged wa = we
flag_na0    = 0;    % Non-averaged wa = 0
flags = {'flag_gurfil','flag_wardw','flag_ward0','flag_naw','flag_na0'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  Define Global Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
project_constants
global eopdata

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

% noradID = 'Gurfil';     no_yrs = 1/365;10;   % Gurfil
% noradID = '00049';      no_yrs = 05;   % Echo 1A (LEO)
% noradID = '02253';      no_yrs = 51;   % PAGEOS-A (Polar)
% noradID = '02324';      no_yrs = 11;   % PasComSat/OV1-8 (LEO)
% noradID = '11659';      no_yrs = 03;   % Ariane 1
% noradID = '16657';      no_yrs = 05;   % Ariane 3 R/B
% noradID = '19218';      no_yrs = 04;   % Ariane 44LP R/B
% noradID = '37239';      no_yrs = 04;   % Ariane 5 R/B
noradID = '37239 mod';  no_yrs = 15;   % Ariane 5 R/B Mod
% noradID = 'ICLEO';      no_yrs = 2;    %LEO

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
    case '16657'  % Ariane 3 R/B
        IC16657   
    case '19218'  % Ariane 44LP R/B
        IC19218
    case '37239'  % Ariane 5 R/B
        IC37239
    case '37239 mod'
        IC37239mod
    case 'ICLEO'
        ICLEO
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Numerical integration parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tf = no_yrs*(365.25*(24*(60*60)));
tspan = [0 tf];
% options = odeset('RelTol',1e-12,'AbsTol',1e-12);
options = odeset('RelTol',1e-12,'AbsTol',1e-12,'Events',@myEvent);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Optimizing eopdata for Non-Averaged Simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
JD = Mjd_UTC_Epoch + 2400000.5;
i_PC = find(PC(:,1)<=JD & JD<=PC(:,2),1,'first');
PC = PC(i_PC:end,:);
mjd = (floor(Mjd_UTC_Epoch));
i_epo = find(mjd==eopdata(4,:),1,'first');
if isempty(i_epo)~=1
    eopdata = eopdata(:,i_epo:end);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       Initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Compute for initial atmospheric density (perigee altitude must be < 1000 km)
%
[rho_p0,H_p0] = atmosphere_og(hp0);         % (A) Input perigee altitude

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
% Define which perturbations to consider
%
%   atm drag   lunisolar   J2
% [   1/0         1/0      1/0 ]
pert_fac = [1 1 1]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Numerically integrate equations of motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---------------------------------------------------------------------------------------------- %
%                                             Gurfil                                             %
% -------------------------------------------------------------------------------------------- %
if flag_gurfil
    %
    % Define drag model: (1) for Gurfil, (2) for Ward
    drag_model = 1;
    %
    % Run numerical integrator
    [t_integrator,xx] = ode113(@(t,x) project_function(t,x,delta,wa,rho_p0,r_p0,H_p0,avg_flag,drag_model,...
                                                        Mjd_UTC_Epoch,pert_fac,PC),tspan,x0,options);
    %
    % Convert integrator output 'xx' from Milankovitch elements to classical orbital elements (coe)
    gurfil = milankovitch2coe(t_integrator,xx);
    clear t_integrator xx
end

% ---------------------------------------------------------------------------------------------- %
%                                        Ward: v_atm = we                                        %
% ---------------------------------------------------------------------------------------------- %
if flag_wardw
    %
    % Define drag model: (1) for Gurfil, (2) for Ward
    drag_model = 2;
    %
    % Run numerical integrator
    [t_integrator,xx] = ode113(@(t,x) project_function(t,x,delta,wa,rho_p0,r_p0,H_p0,avg_flag,drag_model,...
                                                        Mjd_UTC_Epoch,pert_fac,PC),tspan,x0,options);
    %
    % Convert integrator output 'xx' from Milankovitch elements to classical orbital elements (coe)
    wardw = milankovitch2coe(t_integrator,xx);
    clear t_integrator xx
end

% ---------------------------------------------------------------------------------------------- %
%                                         Ward: v_atm = 0                                        %
% ---------------------------------------------------------------------------------------------- %
if flag_ward0
    %
    % Define drag model: (1) for Gurfil, (2) for Ward
    drag_model = 2;
    %
    % Run numerical integrator                                   v--- wa = v_atm = 0
    [t_integrator,xx] = ode113(@(t,x) project_function(t,x,delta,0,rho_p0,r_p0,H_p0,avg_flag,drag_model,...
                                                        Mjd_UTC_Epoch,pert_fac,PC),tspan,x0,options);
    %
    % Convert integrator output 'xx' from Milankovitch elements to classical orbital elements (coe)
    ward0 = milankovitch2coe(t_integrator,xx);
    clear t_integrator xx
end

% ---------------------------------------------------------------------------------------------- %
%                                    Non-averaged: v_atm = we                                    %
% ---------------------------------------------------------------------------------------------- %
if flag_naw
    %
    % Define initial values for integrator as [r0 ; v0]
    x0 = [r0;v0];
    %
    % Run numerical integrator 
    [t,x] = ode113(@(t,x) nonave(t,x,AMR,Cd,r_p0,rho_p0,H_p0,Mjd_UTC_Epoch,pert_fac,wa,PC),tspan,x0,options); 
    %
    % Convert integrator output 'x' from [r;v] to classical orbital elements (coe), and save to naw
    [ecc_p,a_p,incl_p,omega_p,argp_p,nu_p,p_p,eps_p] = rv2coe4vec(x(:,1:3),x(:,4:6),mu_earth);
    naw = [t,a_p,rad2deg(incl_p),rad2deg(omega_p),rad2deg(argp_p),ecc_p,nu_p,p_p,eps_p];
    clear t x
end

% ---------------------------------------------------------------------------------------------- %
%                                    Non-averaged: v_atm = 0                                     %
% ---------------------------------------------------------------------------------------------- %
if flag_na0
    %
    % Define initial values for integrator as [r0 ; v0]
    x0 = [r0;v0];
    %
    % Run numerical integrator                                                      v--- wa = v_atm = 0
    [t,x] = ode113(@(t,x) nonave(t,x,AMR,Cd,r_p0,rho_p0,H_p0,Mjd_UTC_Epoch,pert_fac,0,PC),tspan,x0,options); 
    %
    % Convert integrator output 'x' from [r;v] to classical orbital elements (coe), and save to naw
    [ecc_p,a_p,incl_p,omega_p,argp_p,nu_p,p_p,eps_p] = rv2coe4vec(x(:,1:3),x(:,4:6),mu_earth);
    na0 = [t,a_p,rad2deg(incl_p),rad2deg(omega_p),rad2deg(argp_p),ecc_p,nu_p,p_p,eps_p];
    clear t x
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                           Plot Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
project_plotting_GurfWardReal
% project_plotting_from_scratch

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                        Save Images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q = get(0,'children');

for f = 1:length(q)
    figure(f);
    set(gca,'FontSize',12);
    if strcmp(noradID,'Gurfil')
        figname = sprintf(['GurfWardNoAve' noradID ' ' num2str(no_yrs) 'yrs ' num2str(pert_fac) ' Figure ' num2str(length(q)+1-f)]);
    else
        figname = sprintf(['GurfWardNoAveReal' noradID ' ' num2str(no_yrs) 'yrs ' num2str(pert_fac) ' Figure ' num2str(length(q)+1-f)]);
    end
    if flag_save_figs == 1
        print(q(f),fullfile(pwd,'Figures',figname),'-dpng','-r300');
        savefig(q(f),fullfile(pwd,'Figures',figname));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                        Save Mat Files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Define string of active simulations (simulations with flag == 1)
str_whichsim = [1*flag_gurfil; 2*flag_wardw; 3*flag_ward0; 4*flag_naw; 5*flag_na0; 6*exist('realdata','var')];
str_whichsim = str_whichsim(str_whichsim~=0);
str_sim_all = {'Gurfil'; 'Wardw'; 'Ward0'; 'NAw'; 'NA0'; 'Real'};
mat_name = strjoin(str_sim_all(str_whichsim));
%
% Define string vector of variables to save in mat file
str_whichvar = [1*flag_gurfil; 2*flag_wardw; 3*flag_ward0; 4*flag_naw; 5*flag_na0];
str_whichvar = str_whichvar(str_whichvar~=0);
str_var_all = {'gurfil', 'wardw', 'ward0', 'naw', 'na0'};
str_var = str_var_all(str_whichvar);
%
% Define filename for .mat file
filename = sprintf([ mat_name ' IC_' noradID ' ' num2str(no_yrs,'%.0f') 'yrs ' num2str(pert_fac) ' data']);
%
% Store variables in .mat file, and save to the following folder
if flag_save_matfile
    save(fullfile(pwd,'Data',filename),'noradID','no_yrs','avg_flag','pert_fac',flags{:},str_var{:});
    fprintf('\n\nFile saved as "%s.mat" in ''Data'' folder.\n',filename);
else
    fprintf('\n\nOutput file "%s.mat" not saved.\n',filename);
end

toc

% ODE Stop Condition
function [value, isterminal, direction] = myEvent(t, x)
value      = (([x(1);x(2);x(3)]'*[x(1);x(2);x(3)])^0.5 - 6478137);
isterminal = 1;   % Stop the integration
direction  = 0;
end