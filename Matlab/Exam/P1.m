clc; clear all; close all;

omega3_0 = 0.005;
omega3_f = 0.525;
I3 = 1320;
r_thruster = 1.25;
thruster_impulse = 12.5;

delta_omega3 = omega3_f - omega3_0;

integral_M = I3*delta_omega3;

torque_impulse = 2*thruster_impulse*r_thruster;

N_impluses = ceil(integral_M/torque_impulse);

disp(['The number of impulses needed are: ' num2str(N_impluses)]);