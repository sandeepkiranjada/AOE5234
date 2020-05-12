function nu = trueanomaly(M,e)

En = 0;
Enp1 = M;


while abs(En-Enp1)>0.00000001
    En = Enp1;
    Enp1 = En  - (En-e*sin(En)-M)/(1-e*cos(En));
end


E_f = Enp1;

tan_Nuover2 = sqrt((1+e)/(1-e))*tan(E_f/2);
nu = 2*atan(tan_Nuover2);

end
