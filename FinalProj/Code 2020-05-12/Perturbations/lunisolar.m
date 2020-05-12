function [ dxdt_lunisol ] = lunisolar(t,x,avg_flag,Mjd_UTC_Epoch,PC)

% LUNISOLAR Summary of this function goes here
%   Detailed explanation goes here

global eopdata
global mu_earth
global a_sun mu_sun H_sun_hat_vec h_sun
global a_moon mu_moon H_moon_hat_vec h_moon

%
% Get sun and moon position vector
%
MJD_UTC = Mjd_UTC_Epoch+t/86400;
[UT1_UTC,TAI_UTC] = IERS(eopdata,MJD_UTC,'l');
[~,~, ~, TT_UTC, ~] = timediff(UT1_UTC,TAI_UTC);
MJD_TT = MJD_UTC+TT_UTC/86400;
MJD_TDB = Mjday_TDB(MJD_TT);
%
[~,d_moon_vec_ECI,d_sun_vec_ECI] = JPL_Eph_DE430(MJD_TDB,PC);
%
d_sun = norm(d_sun_vec_ECI);
d_sun_hat_vec_ECI = d_sun_vec_ECI/d_sun;
d_moon = norm(d_moon_vec_ECI);
d_moon_hat_vec_ECI = d_moon_vec_ECI/d_moon;

%
% Define H, Hhat, e, ehat from results of integration step
%
H_vec = [x(1) x(2) x(3)]'; % Angular momentum  
H = norm(H_vec);           % Magnitude of angular momentum
e_vec = [x(4) x(5) x(6)]'; % Eccentricity vector
e = norm(e_vec);           % Eccentricity

%
% Define other parameters in terms of state
%
a = H^2/(mu_earth*(1-e^2));           % Keplerian semi-major axis
h_vec = H_vec./sqrt(mu_earth*a);
n = sqrt(mu_earth/a^3);

%
% Singly Averaged Dynamics
%

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

%
% Doubly Averaged Dynamics
%

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

%
% Triply Averaged Dynamics
%

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


%
% Differential equations for angular momentum (H) and eccentricity vector (e)
%
dHvec = H_dot_sol_vec+H_dot_lun_vec;
devec = e_dot_sol_vec+e_dot_lun_vec;

%
% Save differential equations to output variable for function
%
dxdt_lunisol = [dHvec;devec];


end

