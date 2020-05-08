function [tcoe] = milankovitch2coe(t_integrator,xx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Convert results of integration to Keplerian elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

project_constants

%
% Define H, Hhat, e, ehat from results of integration
%
Hvec = [xx(:,1)   xx(:,2)   xx(:,3)];       % Angular momentum vector
H = Mag([xx(:,1)' ; xx(:,2)' ; xx(:,3)']);  % Magnitude of angular momentum
Hhat = [xx(:,1) xx(:,2) xx(:,3)]./H';       % Angular momentum unit vector
evec = [xx(:,4) xx(:,5) xx(:,6)];           % Eccentricity vector
e = Mag([xx(:,4)' ; xx(:,5)' ; xx(:,6)']);  % Eccentricity
ehat = [xx(:,4) xx(:,5) xx(:,6)]./e';       % Eccentricity unit vector

%
% Preliminary calculations
% 
zhatvec = repmat(zhat',length(t_integrator),1);     % Replicate zhat into 3xn
n = cross(zhatvec,Hvec,2);                          % Line of nodes vector pointing in the direction of the ascending node
n_mag = (sum(n.^2,2)).^(1/2);                       % Calculate magnitude of n
n = n./n_mag;                                       % Line of nodes unit vector pointing in the direction of the ascending node

%
% Define other Keplerian elements
%
a = H.^2./(mu_earth*(1-e.^2));
inc = rad2deg(acos(dot(zhatvec,Hhat,2)));
raan = acosd(n(:,1));
raan(n(:,2) < 0) = 360 - raan(n(:,2) < 0);
argp = acosd(dot(n,evec,2)./e');
argp(evec(:,3) < 0) = 360 - argp(evec(:,3) < 0);

% rp = a.*(1-e);                              % Radius of perigee
% hp = rp-Re;                                 % Perigee altitude
% ra = (2*a-rp);                              % Radius of apogee
% ha = ra-Re;                                 % Apogee altitude

tcoe = [t_integrator,a',inc,raan,argp,e'];