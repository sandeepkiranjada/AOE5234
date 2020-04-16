%% Dynamics of Pitch
clc; close all; clear;
format long

A = 1.29*1e6; % Kg m^2
B = 9.68*1e6; % Kg m^2
C = 10.10*1e6; % Kg m^2
T_orb = 90*60; % Sec

Omega = 2*pi*(1/T_orb);

w3 = sqrt(3*Omega^2*(B-A)/C);

%% Dynamics of Roll and Yaw
% X = [psi1 psi2 psi1_dot psi2_dot]

Sys_mat = [0 0 1 0;...
           0 0 0 1;...
           -(C-B)*Omega^2/A    0  0  -(C-B-A)*Omega/A;...
           0  -4*(C-A)*Omega^2/B  -(B+A-C)*Omega/B  0];
       
e = eig(Sys_mat);
w1 = imag(e(1))
w2 = imag(e(3))
w3
