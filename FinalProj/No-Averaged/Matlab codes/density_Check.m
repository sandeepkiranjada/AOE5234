clear
clc

i=1;
for h=1e3:100:1800*1e3
    hh(i)=h;
%     [ density(i) , H(i) ] = atmosphere_Rosengren(h);
        density(i)  = atmosphere_Rosengren_modified(h);
    i=i+1;
end

figure
plot(hh,log(density))
ylabel('\rho')
% figure
% plot(hh,H)
% ylabel('H')
    