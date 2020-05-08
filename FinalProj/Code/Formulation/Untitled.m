% close all; 
clear; clc;

Max = 1000;
step = 1;
arg = 0:step:Max;

neg_exp = exp(-arg);

besseIs = [besseli(0,arg);...
           besseli(1,arg);...
           besseli(2,arg)];...
           
besseIs_negex = [besseli(0,arg).*neg_exp;...
           besseli(1,arg).*neg_exp;...
           besseli(2,arg).*neg_exp];...

% figure;
% subplot(2,1,1)
% semilogy(arg,neg_exp);
% ylabel('Negative Exponent, exp(-x )')
% legend('e^{-x}')
% axis tight
% % title()
% subplot(2,1,2)
% semilogy(arg,besseIs);
% xlabel('Orbital parameter, x = ae/H_P')
% ylabel('Bessels Functions')
% legend('I_0(x)','I_1(x)','I_2(x)')
% axis tight
% 
% xb = 0:0.01:5;
% besseIs_small = [besseli(0,xb);...
%            besseli(1,xb);...
%            besseli(2,xb)];...
% figure
% plot(xb,besseIs_small);
% xlabel('Argument, x')
% ylabel('Bessels Functions, I_\nu(x)')
% legend('I_0(x)','I_1(x)','I_2(x)')
% title('Modified Bessels Function of First Kind')
% axis tight
% 
% figure;
% plot(arg,besseIs_negex);
% xlabel('Orbital parameter, x = ae/H_P')
% ylabel('Bessels Functions')
% legend('I_0(x)','I_1(x)','I_2(x)')
% axis tight

figure;
subplot(3,1,1:2)
semilogy(arg,besseIs(1,:),'--k','Linewidth',1);hold on;
semilogy(arg,besseIs_negex(1,:),'-k','Linewidth',1);
semilogy(arg,neg_exp,'-.k','Linewidth',1);

ylabel('Terms in Averaged Drag Model by Ward')
legend('I_0(x)','exp(-x)I_0(x)','exp(-x)')
axis tight
subplot(3,1,3)
plot(arg,besseIs_negex(1,:),'-k','Linewidth',1);
xlabel('Orbital parameter, x = ae/H_\rho')
ylabel('Product exp(-x)I_0(x)')
legend('exp(-x)I_0(x)')
axis tight


param = @(a,e,H_p) (a*e/H_p);
z = 425;
e = 0.74
a =  260001e3
% H_p = 43342;
[~,H] = atmosphere_og(z*1e3);
