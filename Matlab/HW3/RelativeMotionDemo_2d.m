clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All units are in meters and meters/sec, instead of mu = 3.986e5 km^3/s^3,
% it is 3.986e14 m^3/s^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is provided 'as-is' for use with AOE 5234: Orbital Mechanics.
% Numerical accuracy is guaranteed only to within the bounds specified by
% The MathWorks Inc.
%
% Author: Andrew Rogers
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gravitational parameter
mu = 3.986e14; 

epsilon = @(xn,xl) norm(xn-xl)/norm(xl); % function computing epsilon.

%% 
Npts = 50;
xdl = linspace(0,10,Npts);
zdl = linspace(0,10,Npts);
[xd0_mat,zd0_mat]=meshgrid(xdl,zdl);
eps_mat = zeros(Npts,Npts);
base_pow = linspace(1,4,Npts);
baseline_length = 10.^base_pow;

for nn=1:Npts
for mm=1:Npts

%% Set up Chief orbit with orbital elements
caseNum   = 'Case 2';
if strcmp(caseNum,'Case 1') == 1
    a      = 7000e3;
    ecc    = 0.0;
    inc    = 45*pi/180;
    raan   = pi/6;
    argper = pi/6;
    f0     = 0;
elseif strcmp(caseNum,'Case 2') == 1
    a      = 15000e3;
    ecc    = 0.4;
    inc    = 45*pi/180;
    raan   = pi/6;
    argper = pi/6;
    f0     = 0;
elseif strcmp(caseNum,'Case 3') == 1
    a      = 45000e3;
    ecc    = 0.8;
    inc    = 45*pi/180;
    raan   = pi/6;
    argper = pi/6;
    f0     = 0;
else
    a      = 6678e3;
    ecc    = 0.025;
    inc    = 28.5*pi/180;
    raan   = 0;
    argper = 0;
    f0     = 0;
end

%% Set up analysis time
% Mean motion
n = sqrt(mu/a^3);

% Initial time
t0 = 0;

% Compute orbital period
period = 2*pi/n;

% Number of orbital periods
numPeriod = 1;

% Total analysis time, given in seconds
tf = numPeriod*period;

% Number of analysis time-steps (defaults to 100 if not specified)
N = 500;

% Time vector for ODE solver
t = linspace(t0,tf,N);

%% Set up deputy initial conditions
% Bounded relative motion constraint factor (dimensionless)
eccFactor = -n*(2+ecc)/(sqrt((1+ecc)*(1-ecc)^3));

% Relative offsets from Chief (meters)
x0 = -300;
y0 = -300;
z0 = 100;

% Relative velocities from Chief (meters/second)
xd0 = xd0_mat(nn,mm);
yd0 = eccFactor*x0;
zd0 = zd0_mat(nn,mm);

% Initial condition vector for deputy satellite
XT_INIT = [x0 y0 z0 xd0 yd0 zd0]';

% Integrator options
options = odeset('RelTol',3e-12,'AbsTol',1e-15,'Stats','on');

%% Numerically integrate equations of motion

% Linear Equations (LERM)
[T1,XL] = ode113(@LinearFormationFlyingEquations,t,XT_INIT,options,mu,a,ecc);

% Nonlinear Equations (NERM)
[T2,XN] = ode113(@NonlinearFormationFlyingEquations,t,XT_INIT,options,mu,a,ecc);

%% Epsilon Calc

eps_mat(nn,mm) = epsilon(XN,XL);

end
end
%% 
figure
surf(xd0_mat,zd0_mat,eps_mat);
xlabel('$\dot{x}_0$','interpreter','latex')
ylabel('$\dot{z}_0$','interpreter','latex')
zlabel('Error, \epsilon')