clc; clear; close all;

mu = 0.01213;
oneDU = 384400 * 1000; % in m
oneTU =  4.348*24*3600; % in s
mu_earth = 1;
% mu_earth = 3.986004418e14; % in m^3/s^2

r1 = @(x,y) sqrt(  (x-mu)^2 + y^2  );
r2 = @(x,y) sqrt(  (x+1-mu)^2 + y^2  );

myfun = @(X) X-[((1-mu)*(X(1)-mu))/r1(X(1),X(2))^3 + (mu*(X(1)+1-mu))/r2(X(1),X(2))^3 ;...
                ((1-mu)*X(2))/r1(X(1),X(2))^3 + (mu*X(2))/r2(X(1),X(2))^3];
            
clc
X_sol = fsolve(myfun,[-1;1])

% r_final = sqrt((X_sol(1)-mu)^2+X_sol(2)^2);

%% Dv for transfer

r_park = 1/60;
r_final = sqrt((X_sol(1)-mu)^2+X_sol(2)^2);



dv_1 = sqrt(mu_earth/r_park) * ( sqrt((2*r_final)/(r_final+r_park)) -1 ) *oneDU/oneTU
dv_2 = sqrt(mu_earth/r_final) * ( 1- sqrt((2*r_park)/(r_final+r_park)) ) *oneDU/oneTU;

dv = dv_1+dv_2;


%% Available Dv

m0 = 1300;
m_star = 1000;
pi_final = m_star/m0;
eps = 0.05;
I_sp = 400;
g_sealevel = 9.806; 
c = I_sp*g_sealevel;

dv_rock = -c*log(eps+pi_final*(1-eps))


