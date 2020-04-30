function [ dxdt_lunisol ] = lunisolar(t,x,mu,~,~,~,~,avg_flag)
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


M_sun = M_sun_epoch + n_sun*t; % nu is M is E for circular orbit. t is since Epoch
nu_sun = M_sun;
d_sun_vec_peri = [r_sun*cos(nu_sun) r_sun*sin(nu_sun) 0]';
d_sun_vec_ECI = R_peri2ECI*d_sun_vec_peri;
d_sun_hat_vec_ECI = d_sun_vec_ECI/norm(d_sun_vec_ECI);


M_moon = M_moon_epoch + n_moon*t; % nu is M is E for circular orbit. t is since Epoch
nu_moon = M_moon;
d_moon_vec_peri = [r_moon*cos(nu_moon) r_moon*sin(nu_moon) 0]';
d_moon_vec_ECI = R_peri2ECI*d_moon_vec_peri;
d_moon_hat_vec_ECI = d_moon_vec_ECI/norm(d_moon_vec_ECI);


if avg_flag == 1
    
    
    
H_dot_sol_vec =  (3*a^2*mu_sun)/( 2*d_sun^3) * ...
                 ( 5*(d_sun_hat_vec_ECI'*e_vec)*cross(e_vec,d_sun_hat_vec_ECI) - ...
                  (d_sun_hat_vec_ECI'*h_vec)*cross(h_vec,d_sun_hat_vec_ECI));

H_dot_lun_vec =  (3*a^2*mu_moon)/( 2*d_moon^3) * ...
                 ( 5*(d_moon_hat_vec_ECI'*e_vec)*cross(e_vec,d_moon_hat_vec_ECI) - ...
                  (d_moon_hat_vec_ECI'*h_vec)*cross(h_vec,d_moon_hat_vec_ECI));

e_dot_sol_vec =  (3*mu_sun)/( 2*n*d_sun^3) * ...
                  ( 5*(d_sun_hat_vec_ECI'*e_vec)*cross(h_vec,d_sun_hat_vec_ECI) - ...
                    (d_sun_hat_vec_ECI'*h_vec)*cross(e_vec,d_sun_hat_vec_ECI) - ...
                    2*cross(h_vec,e_vec) );
                
e_dot_lun_vec =  (3*mu_moon)/( 2*n*d_moon^3) * ...
                  ( 5*(d_moon_hat_vec_ECI'*e_vec)*cross(h_vec,d_moon_hat_vec_ECI) - ...
                    (d_moon_hat_vec_ECI'*h_vec)*cross(e_vec,d_moon_hat_vec_ECI) - ...
                    2*cross(h_vec,e_vec) );
                
                
elseif avg_flag == 2
    
    
    
H_dot_sol_vec =  (3*a^2*mu_sun)/( 2*d_sun^3) * ...
                 ( 5*(d_sun_hat_vec_ECI'*e_vec)*cross(e_vec,d_sun_hat_vec_ECI) - ...
                  (d_sun_hat_vec_ECI'*h_vec)*cross(h_vec,d_sun_hat_vec_ECI));

H_dot_lun_vec = -(3*a^2*mu_moon)/( 4*a_moon^3 * h_moon^3 ) * ...
                 ( 5*(e_vec'*H_moon_hat_vec)*cross(e_vec,H_moon_hat_vec) - ...
                  (h_vec'*H_moon_hat_vec)*cross(h_vec,H_moon_hat_vec));

e_dot_sol_vec =  (3*mu_sun)/( 2*n*d_sun^3) * ...
                  ( 5*(d_sun_hat_vec_ECI'*e_vec)*cross(h_vec,d_sun_hat_vec_ECI) - ...
                    (d_sun_hat_vec_ECI'*h_vec)*cross(e_vec,d_sun_hat_vec_ECI) - ...
                    2*cross(h_vec,e_vec) );
                
e_dot_lun_vec  = -(3*mu_moon)/( 4*n*a_moon^3*h_moon^3) * ...
                  ( 5*(H_moon_hat_vec'*e_vec)*cross(h_vec,H_moon_hat_vec) - ...
                    (h_vec'*H_moon_hat_vec)*cross(e_vec,H_moon_hat_vec) - ...
                    2*cross(h_vec,e_vec) );    
                
elseif avg_flag == 3
    
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
                
                
                
else
    error('Incorrect Value of Wang et al. averaging flag. Set avg_flag to: 1, 2, or 3')
end


%%% Differential equations for angular momentum (H) and eccentricity vector (e)
dHvec = H_dot_sol_vec+H_dot_lun_vec;
devec = e_dot_sol_vec+e_dot_lun_vec;

%%% Save differential equations to output variable for function
dxdt_lunisol = [dHvec;devec];


end

