close all    
clearvars
clc
format long

%% Constants 
Eo = 1;
rad = 10^(-6);
w = 1.37*10^16;
g = 5.32*10^13;
c = physconst('LightSpeed');

%% Vector of frequencies
freq = logspace(log10(0.01*c/rad),log10(10*c/rad),100);        
lamda = c./freq*(2*pi);
ko = 2*pi./lamda;
e_r = (1-w^2./(freq.*(freq+1i*g)));
n_s = sqrt(e_r);

%% Power Calculation
Power_Scat = zeros(length(n_s),1);
Power_Ext = zeros(length(n_s),1);

for i=1:length(n_s)
    a = MieScattering(40,rad,Eo,ko(i),n_s(i));
    [Power_Scat(i),Power_Ext(i)] = a.Power_Efficiency();
end
Power_abs = Power_Ext + Power_Scat;

%% Results
fig=figure(1);
loglog(freq*rad/c,Power_Scat,'ro');hold on;
loglog(freq*rad/c,Power_abs,'bx');
axis tight;
grid on;
xlabel('\omega R/c');
ylabel('Q')
title('Power Efficiency');
legend('Q_{sca}','Q_{abs}','Location','best');
print(fig,'../../../Writing/Documentations/MieSeries/Images/power','-depsc');
