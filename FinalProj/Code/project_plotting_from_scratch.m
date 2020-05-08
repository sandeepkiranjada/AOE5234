%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear
close all

flag_save_figs = 1;
bw = 0;

addpath('./Perturbations v1')
addpath('./Data')
addpath('./Post')
addpath('./../No-Averaged/Matlab codes')
project_constants

%
% Load data
%

%
% Gurfil
%
ICGURFIL
% filename = 'Gurfil Wardw Ward0 IC_Gurfil 10yrs 1  0  0 data'; load(filename);          % Ave    Drag only
% filename = 'Gurfil Wardw Ward0 NAw NA0 IC_Gurfil 10yrs 1  1  1 data'; load(filename); % Ave,NA All perturbations
% filename = 'Gurfil Wardw Ward0 NAw NA0 IC_Gurfil 10yrs 1  1  1 data'; load(filename);flag_na0 = 0; % Ave,NA All perturbations
filename = 'Gurfil Wardw Ward0 NAw NA0 IC_Gurfil 10yrs 1  1  1 data'; load(filename);flag_ward0 = 0; flag_na0 = 0; % Ave,NA All perturbations

% IC02253
% load('Gurfil Wardw Ward0 Real IC_02253 51yrs 1  1  1 data');

% % '02324'  % PasComSat/OV1-8 (LEO)
% IC02324  
% filename = 'Gurfil Wardw Ward0 Real IC_02324 11yrs 1  1  1 data'; load(filename);

%
% '37239' Ariane 5 R/B
%
% IC37239
% % filename = 'Gurfil Wardw Ward0 Real IC_37239 4yrs 1  1  1 data';  load(filename); flag_ward0 = 0; % All perturbations
% filename = 'Gurfil Wardw Ward0 Real IC_37239 4yrs 1  1  1 data';  load(filename); % All perturbations
% load('All_per_AtmosR_Arian5')
% flag_naw = 1;
% naw = [t,a_p,incl_p,omega_p,argp_p,ecc_p,nu_p];
% % filename = 'Gurfil Wardw Ward0 Real IC_37239 10yrs 1  0  0 data'; load(filename); clear realdata % Just Drag 


%{
Input: 
        gurfil = [t_integrator,a,inc,raan,argp,e]
        wardw = [t_integrator,a,inc,raan,argp,e]
        ward0 = [t_integrator,a,inc,raan,argp,e]
        realdata = [YYYY,MM,DD,hh,mm,ss,a,inc,RAAN,argp,ecc,M,rx,ry,rz,vx,vy,vz]
        naw = [t,a,inc,raan,argp,ecc,nu,p,eps];
        na0 = [t,a,inc,raan,argp,ecc,nu,p,eps];
%}

%
% If needed, determine which data to plot.
% The mat file will already have the following flags defined, but we can
% overwrite them here.
% 
% flag_naw = 1;
% flag_na0 = 0;
% flag_wardw = 1;
% flag_ward0 = 0;
% flag_gurfil = 1;

linvec = {'-','-',':','-','--','-'};
if bw
    colrvec = [0.8*[1,1,1] ; 0*[1,1,1] ; 0*[1,1,1] ; 0.4*[1,1,1] ; 0.1*[1,1,1] ; 0*[1,1,1]];
    % linwidvec = [1 1 1 4 1.5 2];
    linwidvec = [2.5 1 1 4 1.5 2];
else
    colour.BLUE = [0,0.447,0.741];
    colour.RED = [0.85,0.325,0.098];
    colour.YELLOW = [0.929,0.694,0.125];
    colour.VIOLET = [0.494,0.184,0.556];
    colour.GREEN = [0.466,0.674,0.188];
    colour.LBLUE = [0.301,0.745,0.933];
    colour.DRED = [0.635,0.078,0.184];
    colrvec = [ colour.BLUE ; colour.DRED ; 0*[1,1,1] ; 0.4*[1,1,1] ; 0.1*[1,1,1] ; 0*[1,1,1]];
    linwidvec = [1 1 1 4 1.5 2];
end
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

%%
%
% Define legend
%
str_whichleg = [1*flag_naw; 2*flag_na0; 3*flag_wardw; 4*flag_ward0; 5*flag_gurfil; 6*exist('realdata','var')];
str_whichleg = str_whichleg(str_whichleg~=0);
str_legend_all = {'Non-averaged \bfv\rm_{atm} = \bf\omega\rm_\oplus', ...
                 'Non-averaged \bfv\rm_{atm} = 0',...
                 'Ward \bfv\rm_{atm} = \bf\omega\rm_\oplus', ...
                 'Ward \bfv\rm_{atm} = 0', ...
                 'Wang', ...
                 'TLE'};
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
xlabel('Elapsed time (yrs)'); ylabel('Semi-major Axis, a (km)');
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
xlabel('Elapsed time (yrs)'); ylabel('Eccentricity, e');
xlim([0 no_yrs]);
legend(str_legend);

%
% Plot inclination
%
f = figure(3);
f.Units = 'centimeters';
% f.Position = picscale*[1 1 10 12];
f.Position = picscale*[1 1 18 12];
if flag_naw,               plot(t_naw,naw(:,3),                         linvec{1},'color',colrvec(1,:),'linewidth',linwidvec(1)); hold on; end
if flag_na0,               plot(t_na0,na0(:,3),                         linvec{2},'color',colrvec(2,:),'linewidth',linwidvec(2)); hold on; end
if flag_wardw,             plot(t_wardw,wardw(:,3),                     linvec{3},'color',colrvec(3,:),'linewidth',linwidvec(3));  hold on; end
if flag_ward0,             plot(t_ward0,ward0(:,3),                     linvec{4},'color',colrvec(4,:),'linewidth',linwidvec(4),'MarkerIndices',[1:100:3000,3000:150:length(t_ward0)],'MarkerSize',mrkrvec(2),'MarkerFaceColor',colrvec(2,:)); hold on; end
if flag_gurfil,            plot(t_gurfil,gurfil(:,3),                   linvec{5},'color',colrvec(5,:),'linewidth',linwidvec(5),'MarkerIndices',[50:100:3000,3000:150:length(t_gurfil)],'MarkerSize',mrkrvec(3),'MarkerFaceColor',colrvec(3,:)); hold on; end
if exist('realdata','var'),plot(t_real,rad2deg(realdata(line_no:end,8)),linvec{6},'color',colrvec(6,:),'linewidth',linwidvec(6)); end
grid on;
xlabel('Elapsed time (yrs)'); ylabel('Inclination, i (deg)');
xlim([0 no_yrs]);
legend(str_legend,'Location','Best');

%
% Plot argp
%
f = figure(4);
f.Units = 'centimeters';
% f.Position = picscale*[10 1 10 12];
f.Position = picscale*[10 1 18 12];
if flag_naw,               plot(t_naw,naw(:,5),                         linvec{1},'color',colrvec(1,:),'linewidth',linwidvec(1)); hold on; end
if flag_na0,               plot(t_na0,na0(:,5),                         linvec{2},'color',colrvec(2,:),'linewidth',linwidvec(2)); hold on; end
if flag_wardw,             plot(t_wardw,wardw(:,5),                     linvec{3},'color',colrvec(3,:),'linewidth',linwidvec(3));  hold on; end
if flag_ward0,             plot(t_ward0,ward0(:,5),                     linvec{4},'color',colrvec(4,:),'linewidth',linwidvec(4),'MarkerIndices',[1:100:3000,3000:150:length(t_ward0)],'MarkerSize',mrkrvec(2),'MarkerFaceColor',colrvec(2,:)); hold on; end
if flag_gurfil,            plot(t_gurfil,gurfil(:,5),                   linvec{5},'color',colrvec(5,:),'linewidth',linwidvec(5),'MarkerIndices',[50:100:3000,3000:150:length(t_gurfil)],'MarkerSize',mrkrvec(3),'MarkerFaceColor',colrvec(3,:)); hold on; end
if exist('realdata','var'),plot(t_real,rad2deg(realdata(line_no:end,10)),linvec{6},'color',colrvec(6,:),'linewidth',linwidvec(6)); end
grid on;
xlabel('Elapsed time (yrs)'); ylabel('Argument of Perigee, \omega (deg)');
xlim([0 no_yrs]);
legend(str_legend,'Location','Best');

%
% Plot raan
%
f = figure(5);
f.Units = 'centimeters';
% f.Position = picscale*[20 1 10 12];
f.Position = picscale*[20 1 18 12];
if flag_naw,               plot(t_naw,naw(:,4),                         linvec{1},'color',colrvec(1,:),'linewidth',linwidvec(1)); hold on; end
if flag_na0,               plot(t_na0,na0(:,4),                         linvec{2},'color',colrvec(2,:),'linewidth',linwidvec(2)); hold on; end
if flag_wardw,             plot(t_wardw,wardw(:,4),                     linvec{3},'color',colrvec(3,:),'linewidth',linwidvec(3));  hold on; end
if flag_ward0,             plot(t_ward0,ward0(:,4),                     linvec{4},'color',colrvec(4,:),'linewidth',linwidvec(4),'MarkerIndices',[1:100:3000,3000:150:length(t_ward0)],'MarkerSize',mrkrvec(2),'MarkerFaceColor',colrvec(2,:)); hold on; end
if flag_gurfil,            plot(t_gurfil,gurfil(:,4),                   linvec{5},'color',colrvec(5,:),'linewidth',linwidvec(5),'MarkerIndices',[50:100:3000,3000:150:length(t_gurfil)],'MarkerSize',mrkrvec(3),'MarkerFaceColor',colrvec(3,:)); hold on; end
if exist('realdata','var'),plot(t_real,rad2deg(realdata(line_no:end,9)),linvec{6},'color',colrvec(6,:),'linewidth',linwidvec(6)); end
grid on;
xlabel('Elapsed time (yrs)'); ylabel('Right Ascension of the Ascending Node, \Omega (deg)');
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
% Define string of perturbations
str_whichprts = [1*pert_fac(1); 2*pert_fac(2); 3*pert_fac(3)];
str_whichprts = str_whichprts(str_whichprts~=0);
str_prts_all = {'drag','lunisolar','J2'};
str2 = strjoin(str_prts_all(str_whichprts));

q = get(0,'children');
q = flip(q);
for f = 1:length(q)
    h = figure(f);
    set(gca,'FontSize',12);
    set(gcf,'color','w');
        figname = sprintf([str1 ' ' str2 ' ' noradID ' Figure ' num2str(f)]);
        if flag_save_figs == 1
            foldername = sprintf(['/IC ' noradID]);
            % print(h,fullfile(pwd,'Figures',foldername,figname),'-deps');         % Save as .eps file
            print(h,fullfile(pwd,'Figures',foldername,figname),'-dpng','-r300'); % Save as bitmap file
            % savefig(h,fullfile(pwd,'Figures',foldername,figname));               % Save as Matlab figure
            clear figname
            if f>2
                h.Position = picscale*[10 1 10 12];
                figname = sprintf([str1 ' ' str2 ' ' noradID ' Figure ' num2str(f) 'b']);
                % print(h,fullfile(pwd,'Figures',foldername,figname),'-deps');         % Save as .eps file
                print(h,fullfile(pwd,'Figures',foldername,figname),'-dpng','-r300'); % Save as bitmap file
                % savefig(h,fullfile(pwd,'Figures',foldername,figname));               % Save as Matlab figure
                clear figname
            end
        end
end