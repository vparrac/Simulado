function A3P1
%Preparación del código
clc
close all
clear all
%%
%Método de euler (RK orden 1)
dt=0.01;
t=0:dt:3000;
y=zeros(1,length(t));
y(1)=0.085;
for i=1:length(t)-1
   y(i+1)=y(i)+dt*fun(t(i),y(i)); 
end
V=y.^(3).*pi().*(4/3);


figure (1)
plot(t,V)
grid on
xlabel('Tiempo [s]')
ylabel('Volumen [m^3]')
title('Solución Método Euler');

%%
%RK orden 4

z=zeros(1,length(t));
z(1)=0.085;
for i=1:length(t)-1
    k1=fun(t(i),z(i));
    k2=fun(t(i)+dt/2,z(i)+(k1/2)*dt);
    k3=fun(t(i)+dt/2,z(i)+(k2/2)*dt);
    k4=fun(t(i)+dt,z(i)+k3*dt);
    z(i+1)=z(i)+(dt/6)*(k1+2*k2+2*k3+k4);
end
Vol=(z.^3).*(4/3).*pi();
figure (2)
plot(t,Vol)
grid on
xlabel('Tiempo [s]')
ylabel('Volumen [m^3]')
title('Solución Método RK4');


end

function resp=fun(t,var)
R=var(1);
%Parámetros
h=300;
Lf=333000;
rho=917;
T=273.15;
Tamb=291.15;
%Función
resp=-(h/(rho*Lf))*(Tamb-T);
end