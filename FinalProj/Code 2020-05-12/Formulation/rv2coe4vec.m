function [e,a,inc,raan,argp,nu,p,eps] = rv2coe4vec(r,v,muu)
%
% Preliminary calculations
%
muu = repmat(muu,size(r,1),1);
rmag = (sum(r.^2,2)).^(1/2);        % Calculate magnitude of r
vmag = (sum(v.^2,2)).^(1/2);        % Calculate magnitude of v
H = cross(r,v,2);                   % Angular momentum vector
Hmag = (sum(H.^2,2)).^(1/2);        % Calculate magnitude of H
k = repmat([0, 0, 1],size(r,1),1);  % Unit vector in the z-direction for ECI
n = cross(k,H,2);                   % Line of nodes vector pointing in the direction of the ascending node
n_mag = (sum(n.^2,2)).^(1/2);       % Calculate magnitude of n
n = n./n_mag;                       % Line of nodes unit vector pointing in the direction of the ascending node
eps = 0.5*vmag.^2 - muu./rmag;      % Specific mechanical energy [DU^2/TU^2]
%
% (1) Eccentricity
%
e_vec = (cross(v,H,2)-muu.*r./rmag)./muu;
e = (sum(e_vec.^2,2)).^(1/2);           % Calculate magnitude of e
%
% (2) Semi-major axis (a) [DU]
%
a = -0.5*muu./eps;
a(e==1) = inf;
%
% (3) Inclination [deg]
%
inc = acosd(dot(k,H,2)./Hmag);
%
% (4) Right ascension of the ascending node
%
raan = acosd(n(:,1));
raan(n(:,2) < 0) = 360 -raan(n(:,2) < 0);
%
% (5) Argument of perigee [deg]
%
argp = acosd(dot(n,e_vec,2)./e);
argp(e_vec(:,3) < 0) = 360 - argp(e_vec(:,3) < 0);
%
% (6) True anomaly [deg]
%
nu = acosd(dot(r,e_vec,2)./(rmag.*e));
nu(dot(r,v,2)<0) = 360 - nu(dot(r,v,2)<0);
%
% (7) Calculate the semilatus rectum
%
p = a.*(1-e.^2);
%
% Convert angles to radians
%
inc = deg2rad(inc);         % Inclination [rad]
raan = deg2rad(raan);       % Right ascension of the ascending node [rad]
argp = deg2rad(argp);       % Argument of periapsis [rad]
nu = deg2rad(nu);           % True anomaly [rad]
end