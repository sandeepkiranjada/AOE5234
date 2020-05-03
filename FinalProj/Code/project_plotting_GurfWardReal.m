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

strvec = {'--k','k',':k','k','--k',':k'};
picscale = 1;

%
% Define time vectors for each solution
%
t_gurfil = (((gurfil(:,1)/60)/60)/24)/365.25;
t_wardw = (((wardw(:,1)/60)/60)/24)/365.25;
t_ward0 = (((ward0(:,1)/60)/60)/24)/365.25;

UTC_datetime = datetime(realdata(line_no:end,1:6));
UTC_datetime = UTC_datetime-UTC_datetime(1);
t_realdata = datenum(UTC_datetime)/365.25;

f = figure(count + 1);
f.Units = 'centimeters';
f.Position = picscale*[1 26 18 12];
plot(t_gurfil  ,gurfil(:,2)*1e-3       ,strvec{1},'linewidth',1); grid on; hold on
plot(t_gurfil  ,gurfil(:,2)*1e-3       ,strvec{1},'linewidth',1); grid on; hold on
plot(t_gurfil  ,gurfil(:,2)*1e-3       ,strvec{1},'linewidth',1); grid on; hold on
plot(t_realdata,realdata(line_no:end,7),strvec{2},'linewidth',1); grid on; hold on
xlabel('Elapsed time (yrs)'); ylabel('a (km)');
xlim([0 no_yrs]);
% ylim([1.6e4 2.5e4]);
% axis equal
% axis([-2 2 -2 2]);