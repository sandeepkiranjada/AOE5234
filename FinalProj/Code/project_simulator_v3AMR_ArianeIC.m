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
flag_gurfil = 0;    % Gurfil
flag_wardw  = 1;    % Ward wa = we
flag_ward0  = 0;    % Ward wa = 0
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

% noradID = 'Gurfil';     no_yrs = 10;3;   % Gurfil
% noradID = '00049';      no_yrs = 05;   % Echo 1A (LEO)
% noradID = '02253';      no_yrs = 51;   % PAGEOS-A (Polar)
% noradID = '02324';      no_yrs = 11;   % PasComSat/OV1-8 (LEO)
% noradID = '11659';      no_yrs = 03;   % Ariane 1
% noradID = '16657';      no_yrs = 05;   % Ariane 3 R/B
% noradID = '19218';      no_yrs = 06;   % Ariane 44LP R/B
noradID = '37239';      no_yrs = 04;10;   % Ariane 5 R/B

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Numerical integration parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tf = no_yrs*(365.25*(24*(60*60)));
tspan = [0 tf];
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
% options = odeset('RelTol',1e-12,'AbsTol',1e-12,'Events',@myEvent);

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
% Define which perturbations to consider
%
%   atm drag   lunisolar   J2
% [   1/0         1/0      1/0 ]
pert_fac = [1 1 1]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Numerically integrate equations of motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---------------------------------------------------------------------------------------------- %
%                                        Ward: v_atm = we                                        %
% ---------------------------------------------------------------------------------------------- %
if flag_wardw
    %
    % Define drag model: (1) for Gurfil, (2) for Ward
    drag_model = 2;
    %
    % Run numerical integrator for different AMRs
    % AMRvec = [0.03 0.02 0.0153 0.012 0.01 0.005]; % S/m: Area to mass ratio of the spacecraft [m^2/kg]
    % AMRvec = [0.02 0.01 0.005];                 % S/m: Area to mass ratio of the spacecraft [m^2/kg]
    % deltavec = linspace(0.0168,0.009,8);
    deltavec = linspace(0.0139,0.0109,9);
    wardw = cell(length(deltavec),1);
    %
    % Change B instead of AMR
    %
    %   for idx = 1:length(AMRvec)
    for idx = 1:length(deltavec)
        
        % delta = 0.5*AMRvec(idx)*Cd;                         % Ballistic coefficient;
        delta = deltavec(idx);
        [t_integrator,xx] = ode113(@(t,x) project_function(t,x,delta,wa,rho_p0,r_p0,H_p0,avg_flag,drag_model,...
            Mjd_UTC_Epoch,pert_fac,PC),tspan,x0,options);
        %
        % Convert integrator output 'xx' from Milankovitch elements to classical orbital elements (coe)
        wardw{idx} = milankovitch2coe(t_integrator,xx);
        %
        % Plot results
        project_plotting_AMR
        
        clear t_integrator xx
        
        AMR = 2*delta/Cd;
        B = 2*delta;
        % legendstr{idx} = sprintf(['AMR = ' num2str(AMR,'%.4f')]);
        legendstr{idx} = sprintf(['B = ' num2str(B,'%.4f') ' AMR = ' num2str(AMR,'%.4f')]);
    end
    
end
legend(legendstr,'Location','Best')
legendstr{idx+1} = 'TLE';
% %{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         Define Legend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
q = get(0,'children');
for f = 1:length(q)
    figure(f);
    grid on
    %%% Modify properties of figures for aesthetics  
    legend(legendstr,'Location','Best')
    set(gca,'FontSize',12);      
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                        Save Images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Define string of perturbations
str_whichprts = [1*pert_fac(1); 2*pert_fac(2); 3*pert_fac(3)];
str_whichprts = str_whichprts(str_whichprts~=0);
str_prts_all = {'drag','lunisolar','J2'};
str2 = strjoin(str_prts_all(str_whichprts));
if flag_save_figs == 1
    for f = 1:length(q)
        figure(f);
        figname = sprintf(['AMR Wardw ' str2 ' Figure ' num2str(length(q)+1-f,'%01.0f')]);
% %         saveas(q(f),fullfile(pwd,'Figures/Ward AMR/EPS',figname),'epsc'); % Save as .eps file
% %         print(q(f),fullfile(pwd,'Figures/Ward AMR/PNG',figname),'-dpng','-r300'); % Save as bitmap file
%         savefig(q(f),fullfile(pwd,'Figures/Ward AMR/FIG',figname));               % Save as Matlab figure
        savefig(q(f),fullfile(pwd,'Figures/Ward B/FIG',figname));               % Save as Matlab figure
    end
    
else
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

% %}