%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_plot = (((t_integrator/60)/60)/24)/365;
strvec = {'k','--k',':k','k','--k',':k'};
str = strvec{idx};
picscale = 1;

if idx <= 3
    count = 0;
else
    count = 8;
end

f = figure(count + 1);
f.Units = 'centimeters';
% f.Position = picscale*[1 26 18 12];
plot(t_plot,a*1e-3,str,'linewidth',1); grid on; hold on
xlabel('Elapsed time (yrs)'); ylabel('a (km)');
xlim([0 no_yrs]);
% ylim([1.6e4 2.5e4]);
% axis equal
% axis([-2 2 -2 2]);

f = figure(count + 2);
f.Units = 'centimeters';
% f.Position = picscale*[1 13.5 18 12];
plot(t_plot,e,str,'linewidth',1); grid on; hold on
xlabel('Elapsed time (yrs)'); ylabel('e ');
xlim([0 no_yrs]);
% ylim([0.58 0.74]);

f = figure(count + 3);
f.Units = 'centimeters';
% f.Position = picscale*[1 1 12 12];
plot(t_plot,inc,str,'linewidth',1); grid on; hold on
xlabel('Elapsed time (yrs)'); ylabel('i (deg)');
xlim([0 no_yrs]);

f = figure(count + 4);
f.Units = 'centimeters';
% f.Position = picscale*[21 26 12 12];
plot(t_plot,argp,str,'linewidth',1); grid on; hold on
xlabel('Elapsed time (yrs)'); ylabel('\omega (deg)');
xlim([0 no_yrs]);

f = figure(count + 5);
f.Units = 'centimeters';
% f.Position = picscale*[21 13.5 12 12];
plot(t_plot,raan,str,'linewidth',1); grid on; hold on
xlabel('Elapsed time (yrs)'); ylabel('\Omega (deg)');
xlim([0 no_yrs]);

f = figure(count + 6);
f.Units = 'centimeters';
% f.Position = picscale*[21 1 18 12];
plot(t_plot,hp*1e-3,str,'linewidth',1); grid on; hold on
xlabel('Elapsed time (yrs)'); ylabel('h_p (km)');
xlim([0 no_yrs]);

f = figure(count + 7);
f.Units = 'centimeters';
% f.Position = picscale*[42 26 18 12];
plot(t_plot,ha*1e-3,str,'linewidth',1); grid on; hold on
xlabel('Elapsed time (yrs)'); ylabel('h_a (km)');
xlim([0 no_yrs]);

f = figure(count + 8);
f.Units = 'centimeters';
% f.Position = picscale*[42 13.5 18 12];
plot(t_plot,H*1e-6,str,'linewidth',1); grid on; hold on
xlabel('Elapsed time (yrs)'); ylabel('H (km^2/s)');
xlim([0 no_yrs]);