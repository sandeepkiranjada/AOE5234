%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [rho,H] = atmosphere_gurfil(h)
%
% ATMOSPHERE calculates density for altitudes from sea level
% through 1000 km using exponential interpolation.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%...Handle altitudes outside of the range:
h=h/1000;
if h > 500
    h = 500;
elseif h < 200
    h = 200;
end

    
lpg_rho = 7.0725e-6*(h-200)*(h-400) - 9.7875e-3*(h-200)-9.595;
rho = 10^lpg_rho;

X= [200 250 300 350 400 500];
Y = [38.70 41.38 44.37 47.82 51.87 62.40];

H = interp1(X,Y,h,'spline');
H=H*1000;
% H = -1/log(10)/(1.4145e-5*h-5.544e-3);
% H = 41.38;

end  %atmopshere
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~