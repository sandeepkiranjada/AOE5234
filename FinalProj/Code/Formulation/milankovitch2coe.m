function [tcoe] = milankovitch2coe(t_integrator,xx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Convert results of integration to Keplerian elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

project_constants

%
% Define H, Hhat, e, ehat from results of integration
%
H = Mag([xx(:,1)' ; xx(:,2)' ; xx(:,3)']);  % Magnitude of angular momentum
Hhat = [xx(:,1) xx(:,2) xx(:,3)]./H';       % Angular momentum unit vector
evec = [xx(:,4) xx(:,5) xx(:,6)];           % Eccentricity vector
e = Mag([xx(:,4)' ; xx(:,5)' ; xx(:,6)']);  % Eccentricity
ehat = [xx(:,4) xx(:,5) xx(:,6)]./e';       % Eccentricity unit vector

%
% Define some constants for convenience
% 
zhatvec = repmat(zhat',length(t_integrator),1)';       % Replicate zhat into 3xn
xhatvec = repmat(xhat',length(t_integrator),1)';       % Replicate xhat into 3xn
c1 = cross(zhatvec,Hhat');                             % zhat x Hhat
c2 = Mag(c1);                                          % |zhat x Hhat|

%
% Define other Keplerian elements
%
a = H.^2./(mu_earth*(1-e.^2));
inc = rad2deg(acos(dot(zhatvec,Hhat')));
raan = rad2deg(asin(dot(xhatvec,Hhat')./c2));
argp = rad2deg(acos(dot(evec',c1)./(e.*c2)));

% rp = a.*(1-e);                              % Radius of perigee
% hp = rp-Re;                                 % Perigee altitude
% ra = (2*a-rp);                              % Radius of apogee
% ha = ra-Re;                                 % Apogee altitude

tcoe = [t_integrator,a,inc,raan,argp,e];