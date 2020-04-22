function [ dxdt_J2 ] = oblate( x,mu,~,Re,~,~ )
%OBLATE Summary of this function goes here
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
h=norm(h_vec);
n = sqrt(mu/a^3);

p_hat_vec = [0 0 1]';

H_dot_J2_vec = - (3*mu*J2*Re^2) / (2*a^3*h^5) * ...
                  (p_hat_vec'*h_vec)*cross(p_hat_vec,h_vec);
              
e_dot_J2_vec = - (3*n*J2*Re^2) / (4*a^2*h^5) * ...
                 ( (1 - (5/h^2)*(p_hat_vec'*h_vec)^2) * cross(h_vec,e_vec) + ...
                  2*(p_hat_vec'*h_vec)*cross(p_hat_vec,e_vec) );
              
              
%%% Differential equations for angular momentum (H) and eccentricity vector (e)
dHvec = H_dot_J2_vec;
devec = e_dot_J2_vec;

%%% Save differential equations to output variable for function
dxdt_J2 = [dHvec;devec];

end

