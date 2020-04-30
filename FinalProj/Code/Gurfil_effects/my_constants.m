%% Oblateness

J2 = 1082.63e-6;



%% Luni-Solar Constants
inc = 23.44*pi/180;

a_sun = 1.495874268632075e11; % Mean value used, circular orbit assumption
V_mean_earth = 29.785731561564607e3;
H_sun_vec = V_mean_earth*a_sun.*[0 -sin(inc) cos(inc)]';
H_sun_hat_vec = [0 -sin(inc) cos(inc)]';
mu_sun = 1.32712440018e20;
h_sun_vec = H_sun_vec./sqrt(a_sun*mu_sun);
h_sun = norm(h_sun_vec);
n_sun = sqrt(mu_sun/a_sun^3);

a_moon = 385e6; % Mean value used, circular orbit assumption
V_mean_moon = 1.022e3;
H_moon_vec = V_mean_moon*a_moon.*[0 -sin(inc) cos(inc)]';
H_moon_hat_vec = [0 -sin(inc) cos(inc)]';
mu_earth = 3.986004418e14;
h_moon_vec = H_moon_vec./sqrt(a_moon*mu_earth);
h_moon = norm(h_moon_vec);
n_moon = sqrt(mu_earth/a_moon^3);


mu_moon = 4.9048695e12;

%% Luni-Solar orbital elements

M_sun_epoch = -((31+28+20)*24+22)*3600*n_sun;% Mean anamoly at epoch = 1st Jan 2015
M_moon_epoch = -((31+28+20)*24+22)*3600*n_moon;% Mean anamoly at epoch = 1st Jan 2015 (incidentally new moon, aligned with Sun)

R_peri2ECI =  [1  0 0;...
               0  cos(inc) -sin(inc);...
               0  sin(inc)  cos(inc)];
           
           
d_sun = a_sun; % Gurfil Notation for Eq 49-52
r_sun = d_sun;

d_moon = a_moon; % Gurfil Notation for Eq 49-52
r_moon = d_moon;



