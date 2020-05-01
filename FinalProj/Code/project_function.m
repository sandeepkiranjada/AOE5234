function  [dxdt] = project_function(t,x,mu,delta,~,~,Re,rho_p0,r_p0,avg_flag,drag_model,Mjd_UTC_Epoch)

if drag_model == 1
    dxdt_drag = atmdrag_gurfil(x,mu,delta,Re,rho_p0,r_p0);
elseif  drag_model == 2
    dxdt_drag = atmdrag_ward(x,mu,delta,Re,rho_p0,r_p0);
else
    error('Incorrect Value drag model. Set drag_model to: 1, or 2')
end
dxdt_lunisol = lunisolar(t,x,mu,delta,Re,rho_p0,r_p0,avg_flag,Mjd_UTC_Epoch);
dxdt_J2 = oblate( x,mu,delta,Re,rho_p0,r_p0 );

%%% Save differential equations to output variable for function
dxdt = dxdt_drag + dxdt_lunisol + dxdt_J2;
% dxdt = dxdt_drag;

end