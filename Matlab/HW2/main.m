%% Problem 1

clear; close all; clc;

mu_e = 1; % Mu earth in DU^3/TU^2

r = [0.33; -0.40; -0.71]'; % in DU
v = [-0.20; -0.20; -0.25]'; % in DU/TU

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

disp('The six orbital elements are');
disp('Semi-major axis (DU):');
disp(a);
disp('Eccentricity (Unit Less):')
disp(e_mag);
disp('Inclination (deg):')
disp(rad2deg(i));
disp('Right Ascension of the ascending node (deg):')
disp(rad2deg(Omega));
disp('Argument of Perigee (deg):')
disp(rad2deg(omega));
disp('True Anomoly (deg):')
disp(rad2deg(nu));


%% Problem 2


mu_s = 1; % Mu sun in AU^3/TU^2
r_e = 1; % Radius of Earth's orbit in AU
r_j = 5.203; % Radius of Jupiter's orbit in AU

a1 = r_e; % circular Earth orbit
a2 = r_j; % circular Jupiter orbit

Delta_v_1 = sqrt(mu_s/a1)*( sqrt(2*a2/(a1+a2)) - 1 ); % to enter transfer orbit
Delta_v_2 = sqrt(mu_s/a2)*( 1 - sqrt(2*a1/(a1+a2)) ); % to circularize the orbit at Jupiter

TOF = pi*sqrt((a1+a2)^3 / (8*mu_s));

n_e = sqrt(mu_s/r_e^3); % mean motion at Earth's orbit
n_j = sqrt(mu_s/r_j^3); % mean motion at Jupiter's orbit

phi_dot = n_e - n_j;

disp(Delta_v_1)
disp(Delta_v_2)
disp(TOF*58.132821)



