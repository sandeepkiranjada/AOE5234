clc; close all; clear;

A = 1.29*1e6; % Kg m^2
B = 9.68*1e6; % Kg m^2
C = 10.10*1e6; % Kg m^2
T_orb = 90*60; % Sec

Omega = 2*pi*(1/T_orb);

psi3 = sqrt(3*Omega^2*(B-A)/C);

%% Dynamics of psi1 and psi2
% X = [psi1 psi2 psi1_dot psi2_dot]

Sys_mat = [0 0 1 0;...
           0 0 0 1;...
           -(C-B)*Omega^2/A    0  0  -(C-B-A)*Omega/A;...
           0  -4*(C-A)*Omega^2/B  -(B+A-C)*Omega/B  0];
       
e = eig(Sys_mat);
w1 = imag(e(1));
w2 = imag(e(3));
tspan = [0 10000];

y0 = [0.1 1 0 0];

[t,y] = ode45(@(t,y) Sys_mat*y, tspan, y0);

figure;
plot(t,y(:,[1:2]));hold on;
plot(t,cos(w1.*t));hold on;
plot(t,cos(w2.*t));hold on;
legend('psi1','psi2','test1','test2')
% figure;
% plot(t,y(:,[3:4]))
% legend