function dx=NA_orbit(t,x,AMR,C_D,r_p0,rho_p,re,H_rho,muu,J2,Mjd_UTC,drag_flag,J2_flag,Lunisolar_flag,atmo_rotation_flag)

global PC eopdata

if atmo_rotation_flag
    W_earth = 7.292115e-5;         % Earth (atmoshere) angular rotation
else
    W_earth = 0;
end
k       = [0;0;1];
w_earth = W_earth*k;
del     = 0.5*AMR*C_D;
r_k     = [x(1);x(2);x(3)];
rdot_k  = [x(4);x(5);x(6)];
r       = [x(7);x(8);x(9)];
rdot    = [x(10);x(11);x(12)];

rhat_k  = r_k/norm(r_k);
rhat    = r/norm(r);
v_k     = norm(rdot_k);
v       = norm(rdot);
H_k     = cross(r_k,rdot_k,1);
H       = cross(r,rdot,1);

%% Perturbed accelerations
a_drag=0; a_j2=0; a_LuniSolar=0; 

%%% Drag
if drag_flag
rho     = rho_p*exp((r_p0-norm(r))/(H_rho));                                           % Density
a_drag  = 1*-(del)*rho*v*(1-dot(w_earth,H,1)/(v^2))*(rdot-cross(w_earth,r,1));         % Drag Perturbation
end

%%% J2
if J2_flag
a_j2    = (-3*J2*muu*(re^2)/(2*(norm(r)^4)))*((1-5*(dot(k,rhat))^2)*rhat+(2*dot(k,rhat))*k);
end

%%% LuniSolar
if Lunisolar_flag
MJD_UTC = Mjd_UTC+t/86400;
[UT1_UTC,TAI_UTC] = IERS(eopdata,MJD_UTC,'l');
[~,~, ~, TT_UTC, ~] = timediff(UT1_UTC,TAI_UTC);
MJD_TT = MJD_UTC+TT_UTC/86400;
MJD_TDB = Mjday_TDB(MJD_TT);

[~,r_Moon,r_Sun] = JPL_Eph_DE430(MJD_TDB);

muu_Sun     = 132712440041.939400e9;      	  	 % [m^3/s^2]; DE430
muu_Moon    = muu/81.30056907419062;              % [m^3/s^2]; DE430
aS = AccelPointMass(r,r_Sun,muu_Sun);
aM = AccelPointMass(r,r_Moon,muu_Moon);

a_LuniSolar=aS+aM;
end

%% States
dx1     = rdot_k;
dx2     = -(muu/(norm(r_k)^3))*r_k;
dx3     = rdot;
dx4     = -(muu/(norm(r)^3))*r + drag_flag*a_drag + J2_flag*a_j2 + Lunisolar_flag*a_LuniSolar;
dx      = [dx1;dx2;dx3;dx4];
end

% w_earth=[0;0;0];

% e_vec = (cross(rdot_k,H_k,1)-muu*r_k/norm(r_k))/muu;
% e = (sum(e_vec.^2,1)).^(1/2);           % Calculate magnitude of e
% nu = acosd(dot(r_k,e_vec,1)./(norm(r_k).*e));
% nu(dot(r_k,rdot_k,1)<0) = 360 - nu(dot(r_k,rdot_k,1)<0);
% if (nu > 359 && nu<=360) || (nu >= 0 && nu<1)
%     h = norm(r_k)-re;                 
%     [rho_p,H_rho,~,~] = atmosphere(h);
% end

% v       = norm(rdot);
% [~,H_p,~,~] = atmosphere(h);      % Input perigee altitude in km
% [~,H_p] = atmosphere_gurfil(h);
% h       = norm(r_k)-re;          
% [ rho , ~ ] = atmosphere_Rosengren(h);                                                   % Density% height [m]
% rho = rho_p*exp((r_p0-norm(r_k))/(H_p));
% rho = rho_p*exp((r_p0-norm(r_k))/(H_rho));
% rho = rho_p;

% f1=0;
% f1=-del*rho*norm((rdot_k-cross(w_earth,r_k,1)))*(rdot_k-cross(w_earth,r_k,1)); 
% dx5=-((x(13)^2)*(1/(2*del))*rho*norm((rdot_k-cross(w_earth,r_k,1)))*dot((rdot_k-cross(w_earth,r_k,1)),rdot_k))/muu;
% f1      =  -del*rho*norm(rdot)*rdot;% Drag Perturbation
% f1=-del*rho*norm((rdot-cross(w_earth,r,1)))*(rdot-cross(w_earth,r,1)); 

% dx = [dx1;dx2;dx3;dx4;dx5];