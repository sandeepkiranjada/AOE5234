clear all; clc; close all;


a1 = 1.1;
a2 = 2;
e1 = 0.05;
e2 = 0.2;
mu = 1;
nu_r = 120*pi/180;

b1 = a1*sqrt(1-e1^2);
b2 = a2*sqrt(1-e2^2);

rp1 = a1*(1-e1);
rp2 = a2*(1-e2);

p1 = a1 * (1-e1^2);
p2 = a2 * (1-e2^2);


k = (rp2/rp1)*(1+e2)/(1+e2*cos(nu_r)) ;

et = (k-1) / (1-k*cos(nu_r));
rpt = rp1;
at = a1*(1-e1)/(1-et);
bt = at*sqrt(1-et^2);
pt = at * (1-et^2);
% at*et

%% Delta vs

Dv_1t = sqrt(2*mu/rpt - mu/at) - sqrt(2*mu/rpt - mu/a1);

V2_atnu_r = sqrt(mu/p2).*[-sin(nu_r) e2+cos(nu_r)];
Vt_atnu_r = sqrt(mu/pt).*[-sin(nu_r) et+cos(nu_r)];
DV_t2 = V2_atnu_r - Vt_atnu_r;

Dv_t2 = norm(DV_t2);

Dv = Dv_1t+Dv_t2;

%% ToF

E_r = acos((et+cos(nu_r))/(1+et*cos(nu_r)));

M_r = E_r - et*sin(E_r);

ToF = M_r*sqrt(at^3/mu);

disp('The Total Delta v, in DU/TU, required is:')
disp(Dv)
disp('The Time of Flight , in TU, to Orbits is:')
disp(ToF)















