clear;
clc;

M_f = -2.31947023336981;
e_mag = 0.728439910910860;
En = 0;
Enp1 = M_f;


while abs(En-Enp1)>0.00000001
    En = Enp1;
    Enp1 = En  - (En-e_mag*sin(En)-M_f)/(1-e_mag*cos(En));
end
E_f = Enp1;
% cosnu = (e_mag - cos(E_f))/(e_mag*cos(E_f)-1);
% sinnu = sin(E_f)*(1+e_mag*cosnu)/sqrt(1-e_mag^2);
% nu_f = atan2(sinnu,cosnu);
tan_Nuover2 = sqrt((1+e_mag)/(1-e_mag))*tan(E_f/2);
nu_f = 2*atan(tan_Nuover2)