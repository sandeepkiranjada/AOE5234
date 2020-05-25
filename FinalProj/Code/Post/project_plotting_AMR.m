%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_integrator = wardw{idx,1}(:,1);
a = wardw{idx,1}(:,2);
inc = wardw{idx,1}(:,3);
raan = wardw{idx,1}(:,4);
argp = wardw{idx,1}(:,5);
e = wardw{idx,1}(:,6);
rp = a.*(1-e);                              % Radius of perigee
hp = rp-Re;                                 % Perigee altitude
ra = (2*a-rp);                              % Radius of apogee
ha = ra-Re;                                 % Apogee altitude

count = 0;

t_plot = (((t_integrator/60)/60)/24)/365;
strvec = {'k','--k',':k','r','--r',':r','b','--b',':b'};
str = strvec{idx};
picscale = 1;

f = figure(count + 1);
f.Units = 'centimeters';
f.Position = picscale*[1 26 16 12];
plot(t_plot,a*1e-3,str,'linewidth',1); grid on; hold on
xlabel('Elapsed Time (yrs)'); ylabel('Semi-major Axis, a (km)');
xlim([0 no_yrs]);
% axis equal
% axis([-2 2 -2 2]);

f = figure(count + 2);
f.Units = 'centimeters';
f.Position = picscale*[1 13.5 16 12];
plot(t_plot,e,str,'linewidth',1); grid on; hold on
xlabel('Elapsed Time (yrs)'); ylabel('Eccentricity, e');
xlim([0 no_yrs]);

f = figure(count + 3);
f.Units = 'centimeters';
f.Position = picscale*[1 1 16 12];
plot(t_plot,inc,str,'linewidth',1); grid on; hold on
xlabel('Elapsed Time (yrs)'); ylabel('Inclination, i (deg)');
xlim([0 no_yrs]);

f = figure(count + 4);
f.Units = 'centimeters';
f.Position = picscale*[21 26 16 12];
plot(t_plot,argp,str,'linewidth',1); grid on; hold on
xlabel('Elapsed Time (yrs)'); ylabel('Argument of Perigee, \omega (deg)');
xlim([0 no_yrs]);

f = figure(count + 5);
f.Units = 'centimeters';
f.Position = picscale*[21 13.5 16 12];
plot(t_plot,raan,str,'linewidth',1); grid on; hold on
xlabel('Elapsed Time (yrs)'); ylabel('Right Ascension of the Ascending Node, \Omega (deg)');
xlim([0 no_yrs]);

f = figure(count + 6);
f.Units = 'centimeters';
f.Position = picscale*[21 1 16 12];
plot(t_plot,hp*1e-3,str,'linewidth',1); grid on; hold on
xlabel('Elapsed Time (yrs)'); ylabel('Perigee Altitude, h_p (km)');
xlim([0 no_yrs]);
% ylim([150 330])

f = figure(count + 7);
f.Units = 'centimeters';
f.Position = picscale*[42 26 16 12];
plot(t_plot,ha*1e-3,str,'linewidth',1); grid on; hold on
xlabel('Elapsed Time (yrs)'); ylabel('Apogee Altitude, h_a (km)');
xlim([0 no_yrs]);

% f = figure(count + 8);
% f.Units = 'centimeters';
% f.Position = picscale*[42 13.5 18 12];
% plot(t_plot,H*1e-6,str,'linewidth',1); grid on; hold on
% xlabel('Elapsed Time (yrs)'); ylabel('H (km^2/s)');
% xlim([0 no_yrs]);

if exist('realdata','var') && idx == length(deltavec)
    UTC_datetime = datetime(realdata(line_no:end,1:6));
    UTC_datetime = UTC_datetime-UTC_datetime(1,:);
    t_real = datenum(UTC_datetime)/365.25;
    f = figure(count + 1);
    plot(t_real,realdata(line_no:end,7),'k','LineWidth',3); 
    
    f = figure(count + 2);
    plot(t_real,realdata(line_no:end,11),'k','LineWidth',3); 
    
    f = figure(count + 3);
    plot(t_real,rad2deg(realdata(line_no:end,8)),'k','LineWidth',3); 
    
    f = figure(count + 4);
    plot(t_real,rad2deg(realdata(line_no:end,10)),'k','LineWidth',3); 
end
