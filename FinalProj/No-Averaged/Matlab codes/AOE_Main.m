clear
clc
% close all
% addpath('C:\Users\ahassani\Documents\GitHub\AOE5234\FinalProj\Code\Data')
addpath('./../../Code/Data')

%% Flags
drag_flag = 1;
J2_flag = 1;
Lunisolar_flag = 1;
atmo_rotation_flag = 1;

%% Constant
mu_earth    = 3.986004418e14;
muu_Sun     = 132712440041.939400e9;      	  	 % [m^3/s^2]; DE430
muu_Moon    = mu_earth/81.30056907419062;        % [m^3/s^2]; DE430
Re          = 6378.137*1e3;                      %  m
J2=0.0010826267;

%% Load data 
fid = fopen('eop19620101.txt','r');
eopdata = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 inf]);
fclose(fid);
load DE430Coeff.mat
PC = DE430Coeff;

% Mjd_UTC = Mjday(2010, 11, 28, 09, 08, 13.23);

%%
%% Spacecraft Initial Conditions
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

Mjd_UTC=Mjd_UTC_Epoch;
JD = Mjd_UTC+2400000.5;
i_PC = find(PC(:,1)<=JD & JD<=PC(:,2),1,'first');
PC = PC(i_PC:end,:);
mjd = (floor(Mjd_UTC));
i_epo = find(mjd==eopdata(4,:),1,'first');
if isempty(i_epo)~=1
    eopdata = eopdata(:,i_epo:end);
end

%%  Calculating Initial Condition for Density 

% [rho_p,H_rho] = atmosphere_gurfil(hp0);  % Input perigee altitude in km
% [rho_p,H_rho] = atmosphere_Rosengren(hp0);
[rho_p,H_rho] = atmosphere_og(hp0); 

%% Solving ODE

no_yrs  = (4);
tf      = no_yrs*(365.25*(24*(60*60)));
tspan   = [0 tf];
options = odeset('RelTol',1e-12,'AbsTol',1e-12,'Events', @myEvent);
x0      = [r0;v0];
[t,x] = ode113(@(t,x) NA_orbit(t,x,AMR,Cd,r_p0,rho_p,Re,H_rho,mu_earth,J2,Mjd_UTC,drag_flag,J2_flag,Lunisolar_flag,atmo_rotation_flag,PC,eopdata,muu_Moon,muu_Sun),tspan,x0,options); 
rmag   = (sum(x(:,1:3).^2,2)).^(1/2);
vmag   = (sum(x(:,4:6).^2,2)).^(1/2);

%% Regenerating COE
[ecc_p,a_p,incl_p,omega_p,argp_p,nu_p,p_p,eps_p] = rv2coe4vec(x(:,1:3),x(:,4:6),mu_earth);
r_p=a_p.*(1-ecc_p);
r_a=a_p.*(1+ecc_p);
h_p=r_p-Re;
h_a=r_a-Re;

incl_p=rad2deg(incl_p);
omega_p=rad2deg(omega_p);
argp_p=rad2deg(argp_p);

%% Saving .mat file

% if drag_flag==1 && J2_flag==0 && Lunisolar_flag==0
%     if atmo_rotation_flag
%         save('Only_drag_AtmosR','ecc_p','a_p','incl_p','omega_p','argp_p','nu_p','r_p','h_p')
%     else
%         save('Only_drag_noAtmosR','ecc_p','a_p','incl_p','omega_p','argp_p','nu_p','r_p','h_p')
%     end
% elseif drag_flag==1 && J2_flag==1 && Lunisolar_flag==1
%     if atmo_rotation_flag
%         save('All_per_AtmosR','ecc_p','a_p','incl_p','omega_p','argp_p','nu_p','r_p','h_p')
%     else
%         save('All_per_noAtmosR','ecc_p','a_p','incl_p','omega_p','argp_p','nu_p','r_p','h_p')
%     end
% end

%% Plot

t_plot = (((t/60)/60)/24)/365.25;
% UTC_datetime = datetime(realdata(line_no:end,1:6));
% UTC_datetime = UTC_datetime-UTC_datetime(1);
% time = datenum(UTC_datetime)/365.25;

figure
plot(t_plot,a_p/1e3,'r','linewidth',1); grid on; hold on
% plot(time,realdata(line_no:end,7),'k','linewidth',1); grid on; hold on
xlabel('Elapsed time (yrs)'); ylabel('a (km)');
% xlim([0 time(end)]);
xlim([0 no_yrs]);

figure
plot(t_plot,ecc_p,'r','linewidth',1); grid on; hold on
% plot(time,realdata(line_no:end,11),'k','linewidth',1); grid on; hold on
xlabel('Elapsed time (yrs)'); ylabel('e ');
% xlim([0 time(end)]);
xlim([0 no_yrs]);

figure
plot(t_plot,incl_p,'r','linewidth',1); grid on; hold on
% plot(time,realdata(line_no:end,8)*180/pi,'k','linewidth',1); grid on; hold on
xlabel('Elapsed time (yrs)'); ylabel('i (deg)');
% xlim([0 time(end)]);
xlim([0 no_yrs]);


%% Ode Stop Condition
function [value, isterminal, direction] = myEvent(t, x)
value      = (([x(1);x(2);x(3)]'*[x(1);x(2);x(3)])^0.5 < 6478137);
isterminal = 1;   % Stop the integration
direction  = 0;
end









%% Removed Lines of Code

% Mjd_UTC = juliandate(2010, 11, 28);

%%% Initial Conditions (Old Versian)
%    h_p0     =  237.6219454075454*1e3;  %1540.577862052572*1e3;  %250*1e3;         %- perigee altitude               m
%    h_a0     =  35730.19283802125*1e3;  %35943*1e3;       %- apogiee altitute               m
%    incl0    =  0.0305510046643136;     %deg2rad(6);      %- inclination                    0.0  to pi rad
%    omega0   =  3.15161424084136;       %deg2rad(60);     %- longitude of ascending node    0.0  to 2pi rad
%    argp0    =  2.94259253726977;       %deg2rad(178);    %- argument of perigee            0.0  to 2pi rad
%    nu0      =  -2.946745376973741;     %0;                      %- true anomaly                   0.0  to 2pi rad
% %%%   
%    r_p0     = (h_p0+re);              %- perigee distance              m
%    r_a0     = h_a0+re;                %- apogiee distance              m
%    sma0     = (r_p0+r_a0)/2;          %- semi-major axis               m
%    ecc0     = 1-r_p0/sma0;            %- eccentricity
%    p0       = sma0*(1-ecc0^2);        %- semilatus rectum              m
%    
%    [r0,v0]  = coe2rv(p0/1e3,ecc0,incl0,omega0,argp0,nu0);
%    r0  = r0*1e3; v0  = v0*1e3;        % m & m/s

%%% Config Constants
% C_D=1.9;    %1.9;      %2.2;      % Drag Coefficient
% AMR=0.0153; %11.056;   %0.02;   % Area to mass ratio 

%%%
% value      = (t > 31536000);
%%

% figure
% plot(t_plot,h_p/1e3,'r','linewidth',1); grid on; hold on
% xlabel('Elapsed time (yrs)'); ylabel('Perigee Height (km)');
% xlim([0 no_yrs]);
% 
% figure
% plot(t_plot,omega_p,'r','linewidth',1); grid on; hold on
% xlabel('Elapsed time (yrs)'); ylabel('\Omega (deg)');
% xlim([0 no_yrs]);
% 
% figure
% plot(t_plot,argp_p,'r','linewidth',1); grid on; hold on
% xlabel('Elapsed time (yrs)'); ylabel('\omega (deg)');
% xlim([0 no_yrs]);

% figure
% plot(t_plot,a_pk/1e3,'b',t_plot,a_p/1e3,'r','linewidth',1); grid on; hold on
% xlabel('Elapsed time (yrs)'); ylabel('a (km)');
% xlim([0 no_yrs]);
% 
% figure
% plot(t_plot,h_pk/1e3,'b',t_plot,h_p/1e3,'r','linewidth',1); grid on; hold on
% xlabel('Elapsed time (yrs)'); ylabel('Perigee Height (km)');
% xlim([0 no_yrs]);
% 
% figure
% plot(t_plot,incl_pk,'b',t_plot,incl_p,'r','linewidth',1); grid on; hold on
% xlabel('Elapsed time (yrs)'); ylabel('i (deg)');
% xlim([0 no_yrs]);
% 
% figure
% plot(t_plot,omega_pk,'b',t_plot,omega_p,'r','linewidth',1); grid on; hold on
% xlabel('Elapsed time (yrs)'); ylabel('\Omega (deg)');
% xlim([0 no_yrs]);
% 
% figure
% plot(t_plot,argp_pk,'b',t_plot,argp_p,'r','linewidth',1); grid on; hold on
% xlabel('Elapsed time (yrs)'); ylabel('\omega (deg)');
% xlim([0 no_yrs]);
% 
% figure
% plot(t_plot,ecc_pk,'b',t_plot,ecc_p,'r','linewidth',1); grid on; hold on
% xlabel('Elapsed time (yrs)'); ylabel('e ');
% xlim([0 no_yrs]);

% figure
% plot(t_plot,h_a/1e3,'r','linewidth',1); grid on; hold on
% xlabel('Elapsed time (yrs)'); ylabel('apogee Height (km)');
% xlim([0 no_yrs]);

