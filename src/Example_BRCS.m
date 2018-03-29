close all    
clearvars
clc
format long

%% Constants 
Eo = 1;
rad = 1;
lamda = 1;
ko = 2*pi/lamda;
j = linspace(0,100,5);
e_r = 1.72-1i*0.0028*j;
n_s = sqrt(e_r);

%% Location
phi = 0;
theta = -pi:pi/180:pi; % Input theta is in radians

%% Power Calculation
BRCS = zeros(length(theta),length(n_s)+1);

for i=1:length(n_s)
    a = MieScattering(40,rad,Eo,ko,n_s(i));
    BRCS(:,i) = a.BistaticRCS(phi,theta);
end
a = MieScattering(40,rad,Eo,ko,0); % n_s = 0 means PEC
BRCS(:,end) = a.BistaticRCS(phi,theta);

%% Results
fig=figure(1);
plot(theta*180/pi,BRCS(:,1),'r');hold on;
plot(theta*180/pi,BRCS(:,2),'b');hold on;
plot(theta*180/pi,BRCS(:,3),'m');hold on;
plot(theta*180/pi,BRCS(:,4),'c');hold on;
plot(theta*180/pi,BRCS(:,5),'g');hold on;
plot(theta*180/pi,BRCS(:,6),'k');
axis tight;
xlabel('\theta^o');
ylabel('BRCS')
title('Bistatic Radar Cross Section');
legend('\epsilon_r = 1.72','\epsilon_r = 1.72-j0.07', ...
          '\epsilon_r = 1.72-j0.14','\epsilon_r = 1.72-j0.21', ...
          '\epsilon_r = 1.72-j0.28','PEC','Location','best');
print(fig,'../../../Writing/Documentations/MieSeries/Images/brcs','-depsc');
