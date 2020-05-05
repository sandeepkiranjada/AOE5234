%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear
close all

flag_save_figs = 1;

addpath('./Perturbations v1')
addpath('./Data')
addpath('./Post')
addpath('./../No-Averaged/Matlab codes')
project_constants

%
% Load data
%
% Gurfil
ICGURFIL
filename = 'Gurfil Wardw Ward0 IC_Gurfil 10yrs 1  0  0 data'; % Drag only
% filename = 'Gurfil Wardw Ward0 IC_Gurfil 10yrs 1  1  1 data'; % All perturbations
% %
% % '00049' Echo 1A (LEO)
% IC00049
% %
% % '02253' PAGEOS-A (Polar)
% IC02253
% filename = 'Gurfil Wardw Ward0 Real IC_02253 51yrs 1  1  1 datal'; % All perturbations
% %
% % '02324' PasComSat/OV1-8 (LEO)
% IC02324
% filename = 'Gurfil Wardw Ward0 Real IC_02324 11yrs 1  1  1 data'; % All perturbations
% %
% % '11659' Ariane 1
% IC11659
% filename = 'Gurfil Wardw Ward0 Real IC_11659 3yrs 1  1  1 data'; % All perturbations
% %
% % '19218' Ariane 44LP R/B
% IC19218
% filename = 'Gurfil Wardw Ward0 Real IC_19218 6yrs 1  1  1 data'; % All perturbations
% % '37239' Ariane 5 R/B
% IC37239
% filename = 'Gurfil Wardw Ward0 Real IC_37239 5yrs 1  1  1 data'; % All perturbations

load(filename)

%{
Input: 
        gurfil = [t_integrator,a,inc,raan,argp,e]
        wardw = [t_integrator,a,inc,raan,argp,e]
        ward0 = [t_integrator,a,inc,raan,argp,e]
        realdata = [YYYY,MM,DD,hh,mm,ss,a,inc,RAAN,argp,ecc,M,rx,ry,rz,vx,vy,vz]
        naw = [t,a,inc,raan,argp,ecc,nu,p,eps];
        na0 = [t,a,inc,raan,argp,ecc,nu,p,eps];
%}

linvec = {'-','-.',':','-','--','-'};
colrvec = [0.8*[1,1,1] ; 0*[1,1,1] ; 0*[1,1,1] ; 0.5*[1,1,1] ; 0.1*[1,1,1] ; 0*[1,1,1]];
linwidvec = [2.5 1 1 4 1.5 2];
mrkrvec = [1 1 1 3 3 2];
picscale = 1;

%
% Define time vectors for each solution
%
if flag_naw,    t_naw    = (((   naw(:,1)/60)/60)/24)/365.25; end
if flag_na0,    t_na0    = (((   na0(:,1)/60)/60)/24)/365.25; end
if flag_wardw,  t_wardw  = ((( wardw(:,1)/60)/60)/24)/365.25; end
if flag_ward0,  t_ward0  = ((( ward0(:,1)/60)/60)/24)/365.25; end
if flag_gurfil, t_gurfil = (((gurfil(:,1)/60)/60)/24)/365.25; end
if exist('realdata','var')
    UTC_datetime = datetime(realdata(line_no:end,1:6));
    UTC_datetime = UTC_datetime-UTC_datetime(1,:);
    t_real = datenum(UTC_datetime)/365.25;
end


%
% Define legend
%
str_whichleg = [1*flag_naw; 2*flag_na0; 3*flag_gurfil; 4*flag_wardw; 5*flag_ward0; 6*exist('realdata','var')];
str_whichleg = str_whichleg(str_whichleg~=0);
str_legend_all = {'Non-averaged \omega = \omega_\oplus', ...
                 'Non-averaged \omega = 0',...
                 'Ward \omega = \omega_\oplus', ...
                 'Ward \omega = 0', ...
                 'Wang', ...
                 'TLE'};
             str_sim_all = {'Gurfil'; 'Wardw'; 'Ward0'; 'NAw'; 'NA0'; 'Real'};
str_legend = str_legend_all(str_whichleg);

%
% Plot semi-major axis
%
f = figure(1);
f.Units = 'centimeters';
f.Position = picscale*[1 22 18 12];
if flag_naw,               plot(t_naw,naw(:,2)*1e-3,           linvec{1},'color',colrvec(1,:),'linewidth',linwidvec(1)); hold on; end
if flag_na0,               plot(t_na0,na0(:,2)*1e-3,           linvec{2},'color',colrvec(2,:),'linewidth',linwidvec(2)); hold on; end
if flag_wardw,             plot(t_wardw,wardw(:,2)*1e-3,       linvec{3},'color',colrvec(3,:),'linewidth',linwidvec(3));  hold on; end
if flag_ward0,             plot(t_ward0,ward0(:,2)*1e-3 ,      linvec{4},'color',colrvec(4,:),'linewidth',linwidvec(4),'MarkerIndices',[1:100:3000,3000:150:length(t_ward0)],'MarkerSize',mrkrvec(2),'MarkerFaceColor',colrvec(2,:)); hold on; end
if flag_gurfil,            plot(t_gurfil,gurfil(:,2)*1e-3,     linvec{5},'color',colrvec(5,:),'linewidth',linwidvec(5),'MarkerIndices',[50:100:3000,3000:150:length(t_gurfil)],'MarkerSize',mrkrvec(3),'MarkerFaceColor',colrvec(3,:)); hold on; end
if exist('realdata','var'),plot(t_real,realdata(line_no:end,7),linvec{6},'color',colrvec(6,:),'linewidth',linwidvec(6)); end
grid on;
xlabel('Elapsed time (yrs)'); ylabel('a (km)');
xlim([0 no_yrs]);
legend(str_legend);

%
% Plot eccentricity
%
f = figure(2);
f.Units = 'centimeters';
f.Position = picscale*[1 13.5 18 12];
if flag_naw,               plot(t_naw,naw(:,6),                 linvec{1},'color',colrvec(1,:),'linewidth',linwidvec(1)); hold on; end
if flag_na0,               plot(t_na0,na0(:,6),                 linvec{2},'color',colrvec(2,:),'linewidth',linwidvec(2)); hold on; end
if flag_wardw,             plot(t_wardw,wardw(:,6),             linvec{3},'color',colrvec(3,:),'linewidth',linwidvec(3));  hold on; end
if flag_ward0,             plot(t_ward0,ward0(:,6),             linvec{4},'color',colrvec(4,:),'linewidth',linwidvec(4),'MarkerIndices',[1:100:3000,3000:150:length(t_ward0)],'MarkerSize',mrkrvec(2),'MarkerFaceColor',colrvec(2,:)); hold on; end
if flag_gurfil,            plot(t_gurfil,gurfil(:,6),           linvec{5},'color',colrvec(5,:),'linewidth',linwidvec(5),'MarkerIndices',[50:100:3000,3000:150:length(t_gurfil)],'MarkerSize',mrkrvec(3),'MarkerFaceColor',colrvec(3,:)); hold on; end
if exist('realdata','var'),plot(t_real,realdata(line_no:end,11),linvec{6},'color',colrvec(6,:),'linewidth',linwidvec(6)); end
grid on;
xlabel('Elapsed time (yrs)'); ylabel('e');
xlim([0 no_yrs]);
legend(str_legend);

%
% Plot inclination
%
f = figure(3);
f.Units = 'centimeters';
f.Position = picscale*[1 1 18 12];
if flag_naw,               plot(t_naw,naw(:,3),                         linvec{1},'color',colrvec(1,:),'linewidth',linwidvec(1)); hold on; end
if flag_na0,               plot(t_na0,na0(:,3),                         linvec{2},'color',colrvec(2,:),'linewidth',linwidvec(2)); hold on; end
if flag_wardw,             plot(t_wardw,wardw(:,3),                     linvec{3},'color',colrvec(3,:),'linewidth',linwidvec(3));  hold on; end
if flag_ward0,             plot(t_ward0,ward0(:,3),                     linvec{4},'color',colrvec(4,:),'linewidth',linwidvec(4),'MarkerIndices',[1:100:3000,3000:150:length(t_ward0)],'MarkerSize',mrkrvec(2),'MarkerFaceColor',colrvec(2,:)); hold on; end
if flag_gurfil,            plot(t_gurfil,gurfil(:,3),                   linvec{5},'color',colrvec(5,:),'linewidth',linwidvec(5),'MarkerIndices',[50:100:3000,3000:150:length(t_gurfil)],'MarkerSize',mrkrvec(3),'MarkerFaceColor',colrvec(3,:)); hold on; end
if exist('realdata','var'),plot(t_real,rad2deg(realdata(line_no:end,8)),linvec{6},'color',colrvec(6,:),'linewidth',linwidvec(6)); end
grid on;
xlabel('Elapsed time (yrs)'); ylabel('i (deg)');
xlim([0 no_yrs]);
legend(str_legend,'Location','Best');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                        Save Images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Define string of active simulations (simulations with flag == 1)
str_whichsim = [1*flag_gurfil; 2*flag_wardw; 3*flag_ward0; 4*flag_naw; 5*flag_na0; 6*exist('realdata','var')];
str_whichsim = str_whichsim(str_whichsim~=0);
str_sim_all = {'Gurfil'; 'Wardw'; 'Ward0'; 'NAw'; 'NA0'; 'Real'};
str1 = strjoin(str_sim_all(str_whichsim));
%
% Define filename for .mat file
str_whichprts = [1*pert_fac(1); 2*pert_fac(2); 3*pert_fac(3)];
str_whichprts = str_whichprts(str_whichprts~=0);
str_prts_all = {'drag','lunisolar','J2'};
str2 = strjoin(str_prts_all(str_whichprts));

q = get(0,'children');
for f = 1:length(q)
    figure(f);
    set(gca,'FontSize',12);
    set(gcf,'color','w');
        figname = sprintf([str1 ' ' str2 ' Figure ' num2str(length(q)+1-f)]);
    if flag_save_figs == 1
%         print(q(f),fullfile(pwd,'Figures',figname),'-deps');         % Save as .eps file
        print(q(f),fullfile(pwd,'Figures',figname),'-dpng','-r300'); % Save as bitmap file
%         savefig(q(f),fullfile(pwd,'Figures',figname));               % Save as Matlab figure
    end
end