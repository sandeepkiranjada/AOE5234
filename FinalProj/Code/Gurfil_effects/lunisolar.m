function [ dxdt_lunisol ] = lunisolar( x,mu,~,~,~,~)
%LUNISOLAR Summary of this function goes here
%   Detailed explanation goes here
my_constants

%%% Define H, Hhat, e, ehat from results of integration step
H_vec = [x(1) x(2) x(3)]'; % Angular momentum  
H = norm(H_vec);           % Magnitude of angular momentum
e_vec = [x(4) x(5) x(6)]'; % Eccentricity vector
e = norm(e_vec);           % Eccentricity 

%%% Define other parameters in terms of state
a = H^2/(mu*(1-e^2));           % Keplerian semi-major axis
h_vec = H_vec./sqrt(mu*a);
n = sqrt(mu/a^3);

H_dot_sol_vec = -(3*a^2*mu_sun)/( 4*a_sun^3 * h_sun^3 ) * ...
                 ( 5*(e_vec'*H_sun_hat_vec)*cross(e_vec,H_sun_hat_vec) - ...
                  (h_vec'*H_sun_hat_vec)*cross(h_vec,H_sun_hat_vec));

H_dot_lun_vec = -(3*a^2*mu_moon)/( 4*a_moon^3 * h_moon^3 ) * ...
                 ( 5*(e_vec'*H_moon_hat_vec)*cross(e_vec,H_moon_hat_vec) - ...
                  (h_vec'*H_moon_hat_vec)*cross(h_vec,H_moon_hat_vec));

e_dot_sol_vec  = -(3*mu_sun)/( 4*n*a_sun^3*h_sun^3) * ...
                  ( 5*(H_sun_hat_vec'*e_vec)*cross(h_vec,H_sun_hat_vec) - ...
                    (h_vec'*H_sun_hat_vec)*cross(e_vec,H_sun_hat_vec) - ...
                    2*cross(h_vec,e_vec) );
                
e_dot_lun_vec  = -(3*mu_moon)/( 4*n*a_moon^3*h_moon^3) * ...
                  ( 5*(H_moon_hat_vec'*e_vec)*cross(h_vec,H_moon_hat_vec) - ...
                    (h_vec'*H_moon_hat_vec)*cross(e_vec,H_moon_hat_vec) - ...
                    2*cross(h_vec,e_vec) );


%%% Differential equations for angular momentum (H) and eccentricity vector (e)
dHvec = H_dot_sol_vec+H_dot_lun_vec;
devec = e_dot_sol_vec+e_dot_lun_vec;

%%% Save differential equations to output variable for function
dxdt_lunisol = [dHvec;devec];


end

