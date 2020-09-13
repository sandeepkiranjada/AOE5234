function  [dxdt] = project_function(t,x,delta,wa,rho_p0,r_p0,H_p0,avg_flag,drag_model,Mjd_UTC_Epoch,pert_fac,PC)
%
% Atmospheric Drag Perturbation
%
if drag_model == 1
    dxdt_drag = atmdrag_gurfil(x,delta,rho_p0,r_p0,H_p0);
elseif  drag_model == 2
    dxdt_drag = atmdrag_ward_lim(x,delta,wa,rho_p0,r_p0,H_p0);
else
    error('Incorrect Value drag model. Set drag_model to: 1, or 2')
end

%
% Sun and Moon Third Body Perturbation
%
dxdt_lunisol = lunisolar(t,x,avg_flag,Mjd_UTC_Epoch,PC);

%
% Earth Oblateness (J2)
%
dxdt_J2 = oblate(x);

%
% Save differential equations to output variable for function
%
dxdt = pert_fac(1)*dxdt_drag + pert_fac(2)*dxdt_lunisol + pert_fac(3)*dxdt_J2;

end