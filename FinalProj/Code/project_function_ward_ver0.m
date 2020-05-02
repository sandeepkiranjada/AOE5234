function  [dxdt] = project_function_ward_ver0(t,x,delta,wa,zhat,rho_p0,rp0,H_p0,flag_J2,flag_lunisolar,jd0)

global Re J2 mu_earth mu_sun mu_moon

%%% Define H, Hhat, e, ehat from results of integration step
H_vec = [x(1); x(2); x(3)];     % Angular momentum vector
e_vec = [x(4); x(5); x(6)];     % Eccentricity vector
%
H = norm(H_vec);                % Magnitude of angular momentum  
e = norm(e_vec);                % Eccentricity
%
Hhat = H_vec/H;                 % Angular momentum unit vector
ehat = e_vec/e;                 % Eccentricity unit vector
    
%%% Define other parameters in terms of the states
a = H^2/(mu_earth*(1-e^2));     % Semi-major axis [m]
rp = a*(1-e);                   % Perigee radius [m]
hp = rp-Re;                     % Perigee altitude [m]
inc = acos(dot(zhat,Hhat));     % Inclination [rad]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  Atmospheric Drag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%% (A) Simply use atmosphere function from Dr. Rosengren
% %%% Compute for atmospheric density (perigee altitude must be < 1000 km)
% [rho,H_p] = atmosphere_og(hp);          % Input perigee altitude in km
% beta = 1/(H_p*1e3);                     % Inverse of scale height [1/m]

% %%% (B) H_p0 is constant
% %%% Compute for atmospheric density (perigee altitude must be < 1000 km)
% rho = rho_p0*exp((rp0-a)/(H_p0));
% beta = 1/(H_p0);                        % Inverse of scale height [1/m]

%%% (C) H_p0 changes
%%% Compute for atmospheric density (perigee altitude must be < 1000 km)
[~,H_p] = atmosphere_og(hp);            % Input perigee altitude in km and compute new scale height
rho = rho_p0*exp((rp0-a)/(H_p));        %
beta = 1/(H_p);                         % Inverse of scale height [1/m]

%%% Compute for Bessel coefficients
K = beta*a*e;
I0 = besseli(0,K);
I1 = besseli(1,K);
I2 = besseli(2,K);

%%% Compute for terms in cross product needed for differential equations
ehatperp = cross(Hhat,ehat);
c1 = ((1+e^2)*I0-2*e*I1)*dot(ehatperp,zhat)*ehat-0.5*(1-e^2)*(I0-I2)*dot(ehat,zhat)*ehatperp;
c2 = 2*beta*a*(1-e^2);
c3 = 0.5*(1-e^2)*(I0-I2)*dot(ehat,zhat)*ehatperp;

%%% Terms needed for the differential equations
dH = -((delta*rho*H^2)/a)*(I0 + e/(2*beta*a*(1-e^2))*I1 - (2*wa*a^2*cos(inc)/H)*((1+e^2)*I0-2*e*I1));
dHhat = delta*rho*wa*a*cross(c1,Hhat);
de = -2*(delta*rho*H/a)*((1-(2-e^2)/c2)*I1 + (1-1/c2)*e*I0 - ((2*wa*a^2*(1-e^2)*cos(inc))/H)*(I1-e*I0));
dehat = -delta*rho*wa*a*cross(c3,Hhat);

%%% Differential equations for angular momentum (H) and eccentricity vector (e)
dHvec = dH*Hhat + H*dHhat;
devec = de*ehat + e*dehat;



if flag_J2 || flag_lunisolar
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                  Parameters for J2 and Lunisolar
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Define other parameters needed for differential equations
    h_vec = H_vec/sqrt(mu_earth*a);     % Scaled angular momentum vector [m]
    h = norm(h_vec);                    % Magnitude of scaled angular momentum vector [m]
    n = sqrt(mu_earth/a^3);             % Mean motion [rad/s]
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                    Earth Oblateness (J2)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if flag_J2
        %
        %%% Compute for J2 perturbation
        dHvec_J2 = -(3*mu_earth*J2*Re^2)/(2*a^3*h^5)*dot(zhat,h_vec)*cross(zhat,h_vec);
        devec_J2 = -(3*n*J2*Re^2)/(4*a^2*h^5)*( ((1-(5/h^2)*dot(zhat,h_vec)^2))*cross(h_vec,e_vec) - (2*(dot(zhat,h_vec)))*cross(zhat,e_vec) );
        %
        %%% Compute for new dHvec, devec
        dHvec = dHvec + dHvec_J2;
        devec = devec + devec_J2;
        %
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                     Lunisolar Third Body
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if flag_lunisolar
        %
        %%% Compute for current time in Julian days
        jd = jd0 + t/86400;
        %
        %%% Compute for unit vectors from Earth to the Sun and Moon
        [r_s,~,~] = sun(jd);            % Ouputs r_sun in AU
        [r_m,~,~] = moon(jd);           % Outputs r_moon in ER
        AU2m = 149598073*1e3;
        ER2m = Re;
        r_s = r_s*AU2m;
        r_m = r_m*ER2m;
        d_s = norm(r_s);
        d_m = norm(r_m);
        d_shat = r_s'/norm(r_s);
        d_mhat = r_m'/norm(r_m);
        %
        %%% Compute for third-body perturbation due to the Sun and Moon
        dHvec_sun =   (3*a^2*mu_sun)/(2*d_s^3)*(5*dot(d_shat,e_vec)*cross(e_vec,d_shat)-dot(d_shat,h_vec)*cross(h_vec,d_shat));
        dHvec_moon = (3*a^2*mu_moon)/(2*d_m^3)*(5*dot(d_mhat,e_vec)*cross(e_vec,d_mhat)-dot(d_mhat,h_vec)*cross(h_vec,d_mhat));
        devec_sun =   (3*mu_sun)/(2*n*d_s^3)*(5*dot(d_shat,e_vec)*cross(h_vec,d_shat)-dot(d_shat,h_vec)*cross(e_vec,d_shat)-2*cross(h_vec,e_vec));
        devec_moon = (3*mu_moon)/(2*n*d_m^3)*(5*dot(d_mhat,e_vec)*cross(h_vec,d_mhat)-dot(d_mhat,h_vec)*cross(e_vec,d_mhat)-2*cross(h_vec,e_vec));
        %
        %%% Compute for new dHvec, devec
        dHvec = dHvec + dHvec_sun + dHvec_moon;
        devec = devec + devec_sun + devec_moon;
    end
end

%%% Save differential equations to output variable for function
dxdt = [dHvec;devec];

end