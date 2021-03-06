clear
clc
% close all

% addpath('C:\Users\ahassani\Documents\GitHub\AOE5234\FinalProj\Code\Data')
addpath('./../../Code/Data')


global PC eopdata
drag_flag = 1;
J2_flag = 1;
Lunisolar_flag = 1;
atmo_rotation_flag = 0;

%% Load Global data 
fid = fopen('eop19620101.txt','r');
eopdata = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 inf]);
fclose(fid);
load DE430Coeff.mat
PC = DE430Coeff;

Mjd_UTC = Mjday(2010, 11, 28, 09, 08, 13.23);

% Step   = 60;   % [s]
% N_Step = 5256000; % 10 years
% num = fix(N_Step*Step/86400)+50;
JD = Mjd_UTC+2400000.5;
i_PC = find(PC(:,1)<=JD & JD<=PC(:,2),1,'first');
PC = PC(i_PC:end,:);
mjd = (floor(Mjd_UTC));
i_epo = find(mjd==eopdata(4,:),1,'first');
if isempty(i_epo)~=1
    eopdata = eopdata(:,i_epo:end);
end

%%
muu=3.986004418e14;
J2=0.0010826267;
%%% Initial Conditions
   h_p0     =  237.6219454075454*1e3;  %1540.577862052572*1e3;  %250*1e3;         %- perigee altitude               m
   h_a0     =  35730.19283802125*1e3;  %35943*1e3;       %- apogiee altitute               m
   incl0    =  0.0305510046643136;     %deg2rad(6);      %- inclination                    0.0  to pi rad
   omega0   =  3.15161424084136;       %deg2rad(60);     %- longitude of ascending node    0.0  to 2pi rad
   argp0    =  2.94259253726977;       %deg2rad(178);    %- argument of perigee            0.0  to 2pi rad
   nu0      =  -2.946745376973741;     %0;                      %- true anomaly                   0.0  to 2pi rad
%%%   
   re       = 6378.137*1e3;           %  m
   r_p0     = (h_p0+re);              %- perigee distance              m
   r_a0     = h_a0+re;                %- apogiee distance              m
   sma0     = (r_p0+r_a0)/2;          %- semi-major axis               m
   ecc0     = 1-r_p0/sma0;            %- eccentricity
   p0       = sma0*(1-ecc0^2);        %- semilatus rectum              m
   
   [r0,v0]  = coe2rv(p0/1e3,ecc0,incl0,omega0,argp0,nu0);
   r0  = r0*1e3; v0  = v0*1e3;        % m & m/s
%%  Config Constants
C_D=2.2;  %1.9;      %2.2;      % Drag Coefficient
AMR=0.0153; %11.056;   %0.02;   % Area to mass ratio 

%%  Calculating Initial Condition for Density 

[rho_p,H_rho] = atmosphere_gurfil(h_p0);  % Input perigee altitude in km
% [rho_p,H_rho] = atmosphere_Rosengren(h_p0);
%% Solving ODE

no_yrs  = (2);
tf      = no_yrs*(365.25*(24*(60*60)));
tspan   = [0 tf];
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
x0      = [r0;v0;r0;v0];
[t,x] = ode113(@(t,x) NA_orbit(t,x,AMR,C_D,r_p0,rho_p,re,H_rho,muu,J2,Mjd_UTC,drag_flag,J2_flag,Lunisolar_flag,atmo_rotation_flag),tspan,x0,options); 
rmag_k = (sum(x(:,1:3).^2,2)).^(1/2);
rmag   = (sum(x(:,7:9).^2,2)).^(1/2);
vmag_k = (sum(x(:,4:6).^2,2)).^(1/2);
vmag   = (sum(x(:,10:12).^2,2)).^(1/2);

% rho = rho_p.*exp((r_p0-rmag)./(H_rho));
%% Regenerating COE
[ecc_pk,a_pk,incl_pk,omega_pk,argp_pk,nu_pk,p_pk,eps_pk] = rv2coe4vec(x(:,1:3),x(:,4:6),muu);
[ecc_p,a_p,incl_p,omega_p,argp_p,nu_p,p_p,eps_p]         = rv2coe4vec(x(:,7:9),x(:,10:12),muu);
r_pk=a_pk.*(1-ecc_pk);
r_p=a_p.*(1-ecc_p);
h_p=r_p-re;
h_pk=r_pk-re;

incl_pk=rad2deg(incl_pk);     incl_p=rad2deg(incl_p);
omega_pk=rad2deg(omega_pk);   omega_p=rad2deg(omega_p);
argp_pk=rad2deg(argp_pk);     argp_p=rad2deg(argp_p);


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

t_plot = (((t/60)/60)/24)/365;


figure
plot(t_plot,a_pk/1e3,'b',t_plot,a_p/1e3,'r','linewidth',1); grid on; hold on
xlabel('Elapsed time (yrs)'); ylabel('a (km)');
xlim([0 no_yrs]);

figure
plot(t_plot,h_pk/1e3,'b',t_plot,h_p/1e3,'r','linewidth',1); grid on; hold on
xlabel('Elapsed time (yrs)'); ylabel('Perigee Height (km)');
xlim([0 no_yrs]);

figure
plot(t_plot,incl_pk,'b',t_plot,incl_p,'r','linewidth',1); grid on; hold on
xlabel('Elapsed time (yrs)'); ylabel('i (deg)');
xlim([0 no_yrs]);

figure
plot(t_plot,omega_pk,'b',t_plot,omega_p,'r','linewidth',1); grid on; hold on
xlabel('Elapsed time (yrs)'); ylabel('\Omega (deg)');
xlim([0 no_yrs]);

figure
plot(t_plot,argp_pk,'b',t_plot,argp_p,'r','linewidth',1); grid on; hold on
xlabel('Elapsed time (yrs)'); ylabel('\omega (deg)');
xlim([0 no_yrs]);

figure
plot(t_plot,ecc_pk,'b',t_plot,ecc_p,'r','linewidth',1); grid on; hold on
xlabel('Elapsed time (yrs)'); ylabel('e ');
xlim([0 no_yrs]);

% figure
% plot(t_plot,rmag_k,'b','linewidth',1); grid on; hold on
% xlabel('Elapsed time (yrs)'); ylabel('r magnitude (Keplerian)');
% xlim([0 no_yrs]);
% 
% figure
% plot(t_plot,rmag,'r','linewidth',1); grid on; hold on
% xlabel('Elapsed time (yrs)'); ylabel('r magnitude (Perturbed)');
% xlim([0 no_yrs]);
% 
% figure
% plot(t_plot,rmag_k-rmag,'k','linewidth',1); grid on; hold on
% xlabel('Elapsed time (yrs)'); ylabel('difference of r magnitude');
% xlim([0 no_yrs]);

% figure
% plot(t_plot,rmag_k,'b',t_plot,rmag,'r','linewidth',1); grid on; hold on
% xlabel('Elapsed time (yrs)'); ylabel('r magnitude (Keplerian and Purturbed)');
% xlim([0 no_yrs]);

% figure
% plot(t_plot,f1mag,'.b','linewidth',1); grid on; hold on
% xlabel('Elapsed time (yrs)'); ylabel('Drag perturbations ');
% xlim([0 no_yrs]);

% figure
% plot(t_plot,eps_pk,'b',t_plot,eps_p,'r','linewidth',1); grid on; hold on
% xlabel('Elapsed time (yrs)'); ylabel('Specific mechanical energy');
% xlim([0 no_yrs]);
% 
% figure
% plot(t_plot,x(:,13),'b','linewidth',1); grid on; hold on
% xlabel('Elapsed time (yrs)'); ylabel('a');
% xlim([0 no_yrs]);



%%

%    arglat  =                 %- argument of latitude      (ci) 0.0  to 2pi rad
%    truelon =                 %- true longitude            (ce) 0.0  to 2pi rad
%    lonper  =                 %- longitude of periapsis    (ee) 0.0  to 2pi rad

% [rho_p,H_rho] = atmosphere_gurfil(h_p0);
%[rho_p,H_rho,~,~] = atmosphere(h_p0);

% x0      = [r0;v0;r0;v0;sma0];

% f1mag=test_func(x,AMR,C_D,r_p0,rho_p,re,H_rho,muu);

% for i=1:size(x,1)
% [p_pk(i,1),a_pk(i,1),ecc_pk(i,1),incl_pk(i,1),omega_pk(i,1),argp_pk(i,1),nu_pk(i,1)] = rv2coe (x(i,1:3),x(i,4:6));
% [p_p(i,1),a_p(i,1),ecc_p(i,1),incl_p(i,1),omega_p(i,1),argp_p(i,1),nu_p(i,1)] = rv2coe (x(i,7:9),x(i,10:12));
% end
%  [e,a,inc,raan,argp,nu,p]