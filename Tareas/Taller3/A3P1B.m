function A3P1B
%Preparación del código
clc
close all
clear all
%Parámetros
h=[50 100 150 200 250 300];
dt=0.01;
t=0:dt:3000;
y0=0.085;
%Método númerico
for j=1:length(h)
    y=zeros(1,length(t));
    y(1)=y0;
    for i=1:length(t)-1
        T2=fun(t(i),y(i),h(j))+(dt/2)*der(t(i),h(j));
        y(i+1)=y(i)+dt*T2;
    end
    sol(:,j)=y;
end
Vol=sol.^3.*(4/3).*pi();
plot(t,Vol);
legend('h=50','h=100','h=150','h=200','h=250','h=300');
title('Solución para los diferentes h');
xlabel('Tiempo [s]')
ylabel('Volumen [m^3]')
end


function resp=fun(t,var,h)
%Parámetros
Lf=333000;
rho=917;
T=273.15;
Tamb=291.15;
%Función
resp=-(h/(rho*Lf))*(Tamb-T);
end

function resp=der(point,hc)
h=1e-6;
resp=(fun(point+h,0,hc)-fun(point,0,hc))/h;
end