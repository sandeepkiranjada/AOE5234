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

t_plot = (((t_integrator/60)/60)/24)/365;
strvec = {'k','--k',':k','k','--k',':k'};
str = strvec{idx};
picscale = 0.8;

if idx <= 3
    count = 0;
else
    count = 8;
end

f = figure(count + 1);
f.Units = 'centimeters';
f.Position = picscale*[1 26 18 12];
plot(t_plot,a*1e-3,str,'linewidth',1); grid on; hold on
xlabel('Elapsed Time (yrs)'); ylabel('a (km)');
xlim([0 no_yrs]);
% axis equal
% axis([-2 2 -2 2]);

f = figure(count + 2);
f.Units = 'centimeters';
f.Position = picscale*[1 13.5 18 12];
plot(t_plot,e,str,'linewidth',1); grid on; hold on
xlabel('Elapsed Time (yrs)'); ylabel('e ');
xlim([0 no_yrs]);

f = figure(count + 3);
f.Units = 'centimeters';
f.Position = picscale*[1 1 12 12];
plot(t_plot,inc,str,'linewidth',1); grid on; hold on
xlabel('Elapsed Time (yrs)'); ylabel('i (deg)');
xlim([0 no_yrs]);

f = figure(count + 4);
f.Units = 'centimeters';
f.Position = picscale*[21 26 12 12];
plot(t_plot,argp,str,'linewidth',1); grid on; hold on
xlabel('Elapsed Time (yrs)'); ylabel('\omega (deg)');
xlim([0 no_yrs]);

f = figure(count + 5);
f.Units = 'centimeters';
f.Position = picscale*[21 13.5 12 12];
plot(t_plot,raan,str,'linewidth',1); grid on; hold on
xlabel('Elapsed Time (yrs)'); ylabel('\Omega (deg)');
xlim([0 no_yrs]);

f = figure(count + 6);
f.Units = 'centimeters';
f.Position = picscale*[21 1 18 12];
plot(t_plot,hp*1e-3,str,'linewidth',1); grid on; hold on
xlabel('Elapsed Time (yrs)'); ylabel('Perigee Altitude, h_p (km)');
xlim([0 no_yrs]);
ylim([150 330])

f = figure(count + 7);
f.Units = 'centimeters';
f.Position = picscale*[42 26 18 12];
plot(t_plot,ha*1e-3,str,'linewidth',1); grid on; hold on
xlabel('Elapsed Time (yrs)'); ylabel('Apogee Altitude, h_a (km)');
xlim([0 no_yrs]);

% f = figure(count + 8);
% f.Units = 'centimeters';
% f.Position = picscale*[42 13.5 18 12];
% plot(t_plot,H*1e-6,str,'linewidth',1); grid on; hold on
% xlabel('Elapsed Time (yrs)'); ylabel('H (km^2/s)');
% xlim([0 no_yrs]);