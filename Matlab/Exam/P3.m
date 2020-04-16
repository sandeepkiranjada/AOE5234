clc; clear all; close all;

Phi_rr = @(n,t) [4-3*cos(n*t) 0 0;...
                6*(sin(n*t)-n*t) 1 0;...
                0 0 cos(n*t)];
            
Phi_vr = @(n,t) [3*n*sin(n*t) 0 0;...
                6*n*(cos(n*t)-1) 0 0;...
                0 0 -n*sin(n*t)];
            
Phi_rv = @(n,t) [sin(n*t) 2*(1-cos(n*t)) 0;...
                2*(cos(n*t)-1) 4*sin(n*t)-3*n*t 0;...
                0 0 sin(n*t)]./n;
            
Phi_vv = @(n,t) [cos(n*t) 2*sin(n*t) 0;...
                -2*sin(n*t) 4*cos(n*t)-3 0;...
                0 0 cos(n*t)];
            
Phi = @(n,t) [Phi_rr(n,t) Phi_rv(n,t);...
              Phi_vr(n,t) Phi_vv(n,t)];
          
r_sat = 6845*1000; % in mts
mu_earth = 3.986004418e14; % in m^3/s^2
t_toDock = 4.5*3600; % in sec
n = sqrt(mu_earth/r_sat^3);

dr_0 = [257 144 89]'; % in m

dv_0 = -inv(Phi_rv(n,t_toDock))*Phi_rr(n,t_toDock)*dr_0;

disp('The required initial velocities are, in m/s:')
disp(dv_0)

% Phi(n,t_toDock)*[dr_0;dv_0];