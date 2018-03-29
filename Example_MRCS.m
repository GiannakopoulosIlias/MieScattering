close all    
clearvars
clc
format long

%% Constants 
Eo = 1;
rad_vector = linspace(0.01,20,1000);
ko = 1;

MRCS = zeros(length(rad_vector),1);

for i = 1:length(rad_vector)
    a = MieScattering(40,rad_vector(i),Eo,ko,0); % n_s = 0 means PEC
    MRCS(i) = a.MonostaticRCS();
end

%% Results
fig=figure(1);
plot(rad_vector,MRCS,'r'); hold on;
plot(rad_vector,ones(length(rad_vector),1),'b')
axis tight;
xlabel('Radius');
ylabel('MRCS')
title('Monostatic Radar Cross Section');
print(fig,'../../../Writing/Documentations/MieSeries/Images/mrcs','-depsc');
