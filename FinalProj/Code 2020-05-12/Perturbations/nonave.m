function dx = nonave(t,x,AMR,Cd,r_p0,rho_p,H_rho,Mjd_UTC,pert_fac,W_earth,PC)

global eopdata
global mu_earth J2 Re
global mu_sun
global mu_moon

drag_flag = pert_fac(1);
Lunisolar_flag = pert_fac(2);
J2_flag = pert_fac(3);

k       = [0;0;1];
w_earth = W_earth*k;
del     = 0.5*AMR*Cd;
r       = [x(1);x(2);x(3)];
rdot    = [x(4);x(5);x(6)];

r_mag   = (r'*r)^0.5;
rhat    = r/r_mag;

%% Perturbed accelerations
a_drag=0; a_j2=0; a_LuniSolar=0;

%%% Drag
if drag_flag
    rho     = rho_p*exp((r_p0-r_mag)/(H_rho));
    a_drag  = -del*rho*norm((rdot-cross(w_earth,r,1)))*(rdot-cross(w_earth,r,1));
end

%%% J2
if J2_flag
    a_j2    = (-3*J2*mu_earth*(Re^2)/(2*(r_mag^4)))*((1-5*(dot(k,rhat))^2)*rhat+(2*dot(k,rhat))*k);
end

%%% LuniSolar
if Lunisolar_flag
    MJD_UTC = Mjd_UTC+t/86400;
    [UT1_UTC,TAI_UTC] = IERS(eopdata,MJD_UTC,'l');
    [~,~, ~, TT_UTC, ~] = timediff(UT1_UTC,TAI_UTC);
    MJD_TT = MJD_UTC+TT_UTC/86400;
    MJD_TDB = Mjday_TDB(MJD_TT);
    [~,r_Moon,r_Sun] = JPL_Eph_DE430(MJD_TDB,PC);
    
    
    aS = AccelPointMass(r,r_Sun,mu_sun);
    aM = AccelPointMass(r,r_Moon,mu_moon);
    
    a_LuniSolar = aS+aM;
end

%% States

dx1     = rdot;
dx2     = -(mu_earth/(r_mag^3))*r + drag_flag*a_drag + J2_flag*a_j2 + Lunisolar_flag*a_LuniSolar;
dx      = [dx1;dx2];
end
