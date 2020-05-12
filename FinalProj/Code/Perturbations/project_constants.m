global Re J2 mu_earth zhat
global a_sun mu_sun H_sun_hat_vec h_sun
global a_moon mu_moon H_moon_hat_vec h_moon

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                        Earth Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Re = 6378*1e3;                              % Radius of the Earth [m]
J2 = 0.0010826267;                          % Geopotential perturbation
mu_earth = 3.986004418e14;                  % Gravitational parameter of the Earth [m^3/s^2]
inc_earth = 23.44*pi/180;                   % Earth's orbit's inclination from the ecliptic
xhat = [1;0;0];                             % Inertial X-direction
yhat = [0;1;0];                             % Inertial Y-direction
zhat = [0;0;1];                             % Inertial Z-direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         Sun Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_sun = 149597870700;                       % Sun sma, mean value used, circular orbit assumption [m]
mu_sun = 132712440041.939400e9;             % Gravitational parameter of the Sun [m^3/s^2]
V_mean_earth = 29.785731561564607e3;        % Mean Earth velocity around the Sun and vice versa [rad/s]
H_sun_vec = V_mean_earth*a_sun.*[0 -sin(inc_earth) cos(inc_earth)]';
H_sun_hat_vec = [0 -sin(inc_earth) cos(inc_earth)]';
h_sun_vec = H_sun_vec./sqrt(a_sun*mu_sun);
h_sun = norm(h_sun_vec);
n_sun = sqrt(mu_sun/a_sun^3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                          Moon Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_moon = 385e6;                             % Moon sma, mean value used, circular orbit assumption
mu_moon = mu_earth/81.30056907419062;       % Gravitational parameter of the Moon [m^3/s^2]
V_mean_moon = 1.022e3;                      % Mean moon velocity around the Earth [rad/s]
H_moon_vec = V_mean_moon*a_moon.*[0 -sin(inc_earth) cos(inc_earth)]';
H_moon_hat_vec = [0 -sin(inc_earth) cos(inc_earth)]';
h_moon_vec = H_moon_vec./sqrt(a_moon*mu_earth);
h_moon = norm(h_moon_vec);
n_moon = sqrt(mu_earth/a_moon^3);