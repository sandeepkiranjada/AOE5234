clear; close all; clc; format long

mu_e = 1; % Mu earth in DU^3/TU^2

r = [-1.41338724296068;-0.161641525116409;0.623463878918531]'; % in DU
v = [0.701979913508663;0.621534620604373;0.123703451458036]'; % in DU/TU

r_mag = sqrt(r*r'); % magnitude of r in DU
v_mag = sqrt(v*v'); % magnitude of v in DU/TU


H = cross(r,v); % Specific Angular momentum in DU^2/TU
H_mag = sqrt(H*H'); % magnitude of v in DU^2/TU

energy = 0.5*v_mag^2 - mu_e/r_mag; % total spcific energy of the Orbit 
                                   % in DU^2/TU^2

a = -mu_e/(2*energy); % Semi-Major Axis in AU

e = cross(v,H)/mu_e - r/r_mag; % e_hat unit less

e_mag = sqrt(e*e'); % Eccentricity unit less

k = [0 0 1]; % unit vector of the earths spin axis

cos_i = k*H'/H_mag;

i = acos(cos_i); % inclination in rad

kcrossH = cross(k,H);
kcrossH_mag = sqrt(kcrossH*kcrossH');

n_hat = kcrossH/kcrossH_mag; % n_hat unit less

ni = n_hat(1);
nj = n_hat(2);

Omega = atan2(nj,ni); % Right Ascension of the ascending node in rad

omega = acos(n_hat*e'/e_mag); % Argument of Perigee in rad

nu = acos(r*e'/(r_mag*e_mag)); % True Anomoly in rad

disp('The six orbital elements are,');
disp('Semi-major axis, a, (DU):');
disp(a);
disp('Eccentricity, e, (Unit Less):')
disp(e_mag);
disp('Inclination, i, (deg):')
disp(rad2deg(i));
disp('Right Ascension of the ascending node, \Omega, (deg):')
disp(rad2deg(Omega));
disp('Argument of Perigee, \omega, (deg):')
disp(rad2deg(omega));
disp('True Anomoly, \nu, (deg):')
disp(rad2deg(nu));
