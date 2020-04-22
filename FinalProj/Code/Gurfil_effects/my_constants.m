%% Oblateness

J2 = 1082.63e-6;



%% Luni-Solar Constants
inc = 23.44*pi/180;

a_sun = 149597870700;
V_mean_earth = 29.78e3;
H_sun_vec = V_mean_earth*a_sun.*[0 -sin(inc) cos(inc)]';
H_sun_hat_vec = [0 -sin(inc) cos(inc)]';
mu_sun = 1.32712440018e20;
h_sun_vec = H_sun_vec./sqrt(a_sun*mu_sun);
h_sun = norm(h_sun_vec);
n_sun = sqrt(mu_sun/a_sun^3);

a_moon = 385e6;
V_mean_moon = 1.022e3;
H_moon_vec = V_mean_moon*a_moon.*[0 -sin(inc) cos(inc)]';
H_moon_hat_vec = [0 -sin(inc) cos(inc)]';
mu_earth = 3.986004418e14;
h_moon_vec = H_moon_vec./sqrt(a_moon*mu_earth);
h_moon = norm(h_moon_vec);
n_moon = sqrt(mu_earth/a_moon^3);


mu_moon = 4.9048695e12;
