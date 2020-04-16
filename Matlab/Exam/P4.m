clear; close all; clc; format long

mu_e = 1; % Mu earth in DU^3/TU^2

r = [-3.85, 1.25, 1.05]; % in DU
v = [0.15, 0.15, -0.2]; % in DU/TU

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

if e*k'< 0
    omega = 2*pi-omega;
end

nu = acos(r*e'/(r_mag*e_mag));% True Anomoly in rad

if r*v'< 0
    nu = 2*pi-nu;
end

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
%% Computing the Change in true anamoly.
% cosE = (e_mag + cos(nu))/(1+e_mag*cos(nu));
% sinE = sqrt(1-e_mag^2)*sin(nu)/(1+e_mag*cos(nu));
% E_e = atan2(sinE,cosE);
tan_Eover2 = sqrt((1-e_mag)/(1+e_mag))*tan(nu/2);
E_e = 2*atan(tan_Eover2);
M_e = E_e-e_mag*sin(E_e);
n = sqrt(mu_e/a^3);

oneTU = 806.7 ; % sec
delT_sec = 2.75*3600; %sec

delT_TU = delT_sec/oneTU;

M_f = M_e + delT_TU*n;

En = 0;
Enp1 = M_f;

while abs(En-Enp1)>0.00000001
    En = Enp1;
    Enp1 = En  - (En-e_mag*sin(En)-M_f)/(1-e_mag*cos(En));
end
E_f = Enp1;
% cosnu = (e_mag - cos(E_f))/(e_mag*cos(E_f)-1);
% sinnu = sin(E_f)*(1+e_mag*cosnu)/sqrt(1-e_mag^2);
% nu_f = atan2(sinnu,cosnu);
tan_Nuover2 = sqrt((1+e_mag)/(1-e_mag))*tan(E_f/2);
nu_f = 2*atan(tan_Nuover2);

%% Computing r and v

% r_f = a*(1-e_mag^2)/(1+e_mag*cos(nu_f));

% r_f_peri =  [r_f*cos(nu_f) r_f*sin(nu_f)];
% v_f_peri =  sqrt(mu_e/(a*(1-e_mag^2)))*[-sin(nu_f) e_mag+cos(nu_f)];
nu=nu_f
r_f = a*(1-e_mag^2)/(1+e_mag*cos(nu));
r_f_peri =  [r_f*cos(nu) r_f*sin(nu)];
v_f_peri =  sqrt(mu_e/(a*(1-e_mag^2)))*[-sin(nu) e_mag+cos(nu)];

R3 = @(x) [cos(x) -sin(x) 0; sin(x) cos(x) 0; 0 0 1]';
R1 = @(x) [1 0 0;0 cos(x) -sin(x);0 sin(x) cos(x)]';

R_pi = R3(-Omega)*R1(-i)*R3(-omega);

r_f_vec = R_pi*[r_f_peri';0]
v_f_vec = R_pi*[v_f_peri';0]


%% Check to see if COEs are presereved except Nu



mu_e = 1; % Mu earth in DU^3/TU^2

r = r_f_vec'; % in DU
v = v_f_vec'; % in DU/TU

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

if e*k'< 0
    omega = 2*pi-omega;
end

nu = acos(r*e'/(r_mag*e_mag));% True Anomoly in rad

if r*v'< 0
    nu = 2*pi-nu;
end

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
