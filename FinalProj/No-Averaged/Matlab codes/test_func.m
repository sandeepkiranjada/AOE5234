function f1mag=test_func(x,AMR,C_D,r_p0,rho_p,re,H_rho,muu)

w_earth=[0;0;7.292115e-5];
w_earth=repmat(w_earth',length(x),1);
% w_earth=[0;0;0];
del=0.5*AMR*C_D;
r_k=x(:,1:3);
rdot_k=x(:,4:6);
rmag=(sum(r_k.^2,2)).^(1/2);
v_k=(sum(rdot_k.^2,2)).^(1/2);
H_k=cross(r_k,rdot_k,2);
% h = x(1)-re;                      % height [m]
% [~,H_p,~,~] = atmosphere(h);      % Input perigee altitude in km
% [~,H_p] = atmosphere_gurfil(h);
% rho = rho_p*exp((r_p0-norm(r_k))/(H_p));
rho = rho_p.*exp((repmat(r_p0,length(x),1)-rmag)./(H_rho));
% rho = rho_p;
f1=-del.*rho.*v_k.*(1-sum((w_earth.*H_k),2)./(v_k.^2)).*(rdot_k-cross(w_earth,r_k,2)); % Drag Perturbation
f1mag=(sum(f1.^2,2)).^(1/2);
% f1=-del*rho*v*(rdot_k-cross(w_earth,r_k)); 

end