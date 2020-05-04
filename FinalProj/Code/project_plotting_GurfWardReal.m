%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
Input: 
        gurfil = [t_integrator,a,inc,raan,argp,e]
        wardw = [t_integrator,a,inc,raan,argp,e]
        ward0 = [t_integrator,a,inc,raan,argp,e]
        realdata = [YYYY,MM,DD,hh,mm,ss,a,inc,RAAN,argp,ecc,M,rx,ry,rz,vx,vy,vz]
%}

close all
strvec = {':k','-k','--k','k','--k',':k'};
colrvec = [0*[1,1,1] ; 0.5*[1,1,1] ; 0.1*[1,1,1] ; 0*[1,1,1]];
linwidvec = [1 3 1.5 1];
mrkrvec = [1 3 3 1];
picscale = 1;

%
% Define time vectors for each solution
%
t_gurfil = (((gurfil(:,1)/60)/60)/24)/365.25;
t_wardw = (((wardw(:,1)/60)/60)/24)/365.25;
t_ward0 = (((ward0(:,1)/60)/60)/24)/365.25;

UTC_datetime = datetime(realdata(line_no:end,1:6));
UTC_datetime = UTC_datetime-UTC_datetime(1);
t_real = datenum(UTC_datetime)/365.25;

f = figure(1);
f.Units = 'centimeters';
f.Position = picscale*[1 22 18 12];
plot(t_wardw,wardw(:,2)*1e-3       ,strvec{1},'color',colrvec(1,:),'linewidth',linwidvec(1)); grid on; hold on
plot(t_ward0,ward0(:,2)*1e-3       ,strvec{2},'color',colrvec(2,:),'linewidth',linwidvec(2),'MarkerIndices',[1:100:3000,3000:150:length(t_ward0)],'MarkerSize',mrkrvec(2),'MarkerFaceColor',colrvec(2,:));
plot(t_gurfil,gurfil(:,2)*1e-3     ,strvec{3},'color',colrvec(3,:),'linewidth',linwidvec(3),'MarkerIndices',[50:100:3000,3000:150:length(t_gurfil)],'MarkerSize',mrkrvec(3),'MarkerFaceColor',colrvec(3,:)); 
plot(t_real,realdata(line_no:end,7),strvec{4},'color',colrvec(4,:),'linewidth',linwidvec(4));
xlabel('Elapsed time (yrs)'); ylabel('a (km)');
xlim([0 no_yrs]);
% ylim([1.6e4 2.5e4]);
legend('Ward \omega = \omega_\oplus','Ward \omega = 0','Wang','TLE');

f = figure(2);
f.Units = 'centimeters';
f.Position = picscale*[1 13.5 18 12];
plot(t_wardw,wardw(:,6)             ,strvec{1},'color',colrvec(1,:),'linewidth',linwidvec(1)); grid on; hold on
plot(t_ward0,ward0(:,6)             ,strvec{2},'color',colrvec(2,:),'linewidth',linwidvec(2),'MarkerIndices',[1:100:3000,3000:150:length(t_ward0)],'MarkerSize',mrkrvec(2),'MarkerFaceColor',colrvec(2,:));
plot(t_gurfil,gurfil(:,6)           ,strvec{3},'color',colrvec(3,:),'linewidth',linwidvec(3),'MarkerIndices',[50:100:3000,3000:150:length(t_gurfil)],'MarkerSize',mrkrvec(3),'MarkerFaceColor',colrvec(3,:)); 
plot(t_real,realdata(line_no:end,11),strvec{4},'color',colrvec(4,:),'linewidth',linwidvec(4));
xlabel('Elapsed time (yrs)'); ylabel('e');
xlim([0 no_yrs]);
% ylim([1.6e4 2.5e4]);
legend('Ward \omega = \omega_\oplus','Ward \omega = 0','Wang','TLE');

f = figure(3);
f.Units = 'centimeters';
f.Position = picscale*[1 1 18 12];
plot(t_wardw,wardw(:,3)             ,strvec{1},'color',colrvec(1,:),'linewidth',linwidvec(1)); grid on; hold on
plot(t_ward0,ward0(:,3)             ,strvec{2},'color',colrvec(2,:),'linewidth',linwidvec(2),'MarkerIndices',[1:100:3000,3000:150:length(t_ward0)],'MarkerSize',mrkrvec(2),'MarkerFaceColor',colrvec(2,:));
plot(t_gurfil,gurfil(:,3)           ,strvec{3},'color',colrvec(3,:),'linewidth',linwidvec(3),'MarkerIndices',[50:100:3000,3000:150:length(t_gurfil)],'MarkerSize',mrkrvec(3),'MarkerFaceColor',colrvec(3,:)); 
realdatainc = realdata(:,8)*180/pi;
plot(t_real,realdatainc(line_no:end),strvec{4},'color',colrvec(4,:),'linewidth',linwidvec(4));
xlabel('Elapsed time (yrs)'); ylabel('i (deg)');
xlim([0 no_yrs]);
legend('Ward \omega = \omega_\oplus','Ward \omega = 0','Wang','TLE','Location','Southeast');