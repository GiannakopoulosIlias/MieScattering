close all    
clearvars
clc
format long

%% Constants
Eo = 1;
ko = 1;
R = 0.5;
e_r = 20-1i*29.76;
n_s = sqrt(e_r);

%% Create Object
a = MieScattering(10,R,Eo,ko,n_s);

%% Axis
x = linspace(-0.48,0.48,100);
y = linspace(-0.48,0.48,100);
z = linspace(-0.48,0.48,100);

%% Calculation of Fields
[E,H] = a.Interior_Fields_Line(x,0,0);
Exx = abs(E(:,1));
Eyx = abs(E(:,2));
Ezx = abs(E(:,3));
Hxx = abs(H(:,1));
Hyx = abs(H(:,2));
Hzx = abs(H(:,3));
[E,H] = a.Interior_Fields_Line(0,y,0);
Exy = abs(E(:,1));
Eyy = abs(E(:,2));
Ezy = abs(E(:,3));
Hxy = abs(H(:,1));
Hyy = abs(H(:,2));
Hzy = abs(H(:,3));
[E,H] = a.Interior_Fields_Line(0,0,z);
Exz = abs(E(:,1));
Eyz = abs(E(:,2));
Ezz = abs(E(:,3));
Hxz = abs(H(:,1));
Hyz = abs(H(:,2));
Hzz = abs(H(:,3));

%% Results
set(gcf, 'Position', get(0, 'Screensize'));
fig=figure(1);
subplot(2,3,1)
plot(x,Exx,'r*');hold on;
plot(x,Eyx,'ro');hold on;
plot(x,Ezx,'r+');
axis tight;
xlabel('x');
ylabel('|E|');
title('|E| on x axis');
legend('E_x','E_y','E_z','Location','best');
subplot(2,3,4)
plot(x,Hxx,'b*');hold on;
plot(x,Hyx,'bo');hold on;
plot(x,Hzx,'b+');
title('|E| on x axis');
axis tight;
xlabel('x');
ylabel('|H|');
title('|H| on x axis');
legend('H_x','H_y','H_z','Location','best');
subplot(2,3,2)
plot(y,Exy,'r*');hold on;
plot(y,Eyy,'ro');hold on;
plot(y,Ezy,'r+');
axis tight;
xlabel('y');
ylabel('|E|');
title('|E| on y axis');
legend('E_x','E_y','E_z','Location','best');
subplot(2,3,5)
plot(y,Hxy,'b*');hold on;
plot(y,Hyy,'bo');hold on;
plot(y,Hzy,'b+');
axis tight;
xlabel('y');
ylabel('|H|');
title('|H| on y axis');
legend('H_x','H_y','H_z','Location','best');
subplot(2,3,3)
plot(z,Exz,'r*');hold on;
plot(z,Eyz,'ro');hold on;
plot(z,Ezz,'r+');
axis tight;
xlabel('z');
ylabel('|E|');
title('|E| on z axis');
legend('E_x','E_y','E_z','Location','best');
subplot(2,3,6)
plot(z,Hxz,'b*');hold on;
plot(z,Hyz,'bo');hold on;
plot(z,Hzz,'b+');
axis tight;
xlabel('z');
ylabel('|H|');
title('|H| on z axis');
legend('H_x','H_y','H_z','Location','best');
print(fig,'../../../Writing/Documentations/MieSeries/Images/lines','-depsc');