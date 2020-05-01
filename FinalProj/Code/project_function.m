function  [dxdt] = project_function(t,x,delta,wa,rho_p0,r_p0,H_p0,avg_flag,drag_model,Mjd_UTC_Epoch)
%
% Atmospheric Drag Perturbation
%
if drag_model == 1
    dxdt_drag = atmdrag_gurfil(x,delta,rho_p0,r_p0,H_p0);
elseif  drag_model == 2
    dxdt_drag = atmdrag_ward(x,delta,wa,rho_p0,r_p0,H_p0);
else
    error('Incorrect Value drag model. Set drag_model to: 1, or 2')
end

%
% Sun and Moon Third Body Perturbation
%
dxdt_lunisol = lunisolar(t,x,avg_flag,Mjd_UTC_Epoch);

%
% Earth Oblateness (J2)
%
dxdt_J2 = oblate(x);

%
% Save differential equations to output variable for function
%
dxdt = dxdt_drag + dxdt_lunisol + dxdt_J2;

end