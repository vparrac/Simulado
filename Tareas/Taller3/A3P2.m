function A3P2
clc
clear all
close all
%%
%ODE45
span=linspace(0,96.43,100);
cI=[2;0;0;0];
[xspan,ysol]=ode45(@funode,span,cI);

figure (1)
plot(xspan,ysol');
legend('C','P','A','O')
xlabel('Espacio-Tiempo [min]')
ylabel('Flujo mol/min');
title('Concentracion en función del espacio-tiempo- ODE45')
%%
%RK Orden 4

a=0;
b=96.43;
N=1000;
h=(b-a)/N;

tspan=a:h:b;

y(1,1)=2;
y(2,1)=0;
y(3,1)=0;
y(4,1)=0;

for i=1:length(tspan)-1
    k1=funode(tspan(i),y(:,i));
    k2=funode(tspan(i)+h/2,y(:,i)+k1/2*h);
    k3=funode(tspan(i)+h/2,y(:,i)+k2/2*h);
    k4=funode(tspan(i)+h/2,y(:,i)+k3*h);
    y(:,i+1)=y(:,i)+(k1+2*k2+2*k3+k4)*(h/6);
end


figure(2)
plot(tspan,y);
legend('C','P','A','O')
xlabel('Espacio-Tiempo [min]')
ylabel('Flujo mol/min');
title('Concentracion en función del espacio-tiempo-RK4')

%% RK orden 1
a1=0;
b1=96.43;
N1=1000;
h1=(b1-a1)/N1;

tspan1=a1:h1:b1;

y1(1,1)=2;
y1(2,1)=0;
y1(3,1)=0;
y1(4,1)=0;

for i=1:length(tspan1)-1
     y1(:,i+1)=y1(:,i)+h*funode(tspan1(i),y1(:,i));
end

figure (3)
plot(tspan1,y1)
legend('C','P','A','O')
xlabel('Espacio-Tiempo [min]')
ylabel('Flujo mol/min');
title('Concentracion en función del espacio-tiempo-Euler')


%% Taylor de orden 2
dt=0.01;
a2=0;
b2=96.43;
time=a2:dt:b2;
y2(1,1)=0.085;
y2(2,1)=0;
y2(3,1)=0;
y2(4,1)=0;

for i=1:length(time)-1
    T2=funode(time(i),y2(:,i))+(dt/2)*der(time(i),y2(:,i));  
    y2(:,i+1)=y2(:,i)+dt*T2;
end
figure (4)
plot(time,y2)
legend('C','P','A','O')
xlabel('Espacio-Tiempo [min]')
ylabel('Flujo mol/min');
title('Concentracion en función del espacio-tiempo- T2')

end


function resp=funode(t,var)
%Flows
Fc=var(1);
Fp=var(2);
Fa=var(3);
Fo=var(4);

%Constants
k1=0.12;
k2=0.046;
k3=0.02;
k4=0.034;
k5=0.04;
v0=8.5;

%Concentrations
Cc=Fc/v0;
Cp=Fp/v0;
Ca=Fa/v0;
Co=Fo/v0;

%Equations
resp(1,1)=v0*(-k1*Cc);
resp(2,1)=v0*(k1*Cc+k4*Ca-k3*Cp);
resp(3,1)=v0*(k2*Cc+k3*Cp-k4*Ca-k5*Ca);
resp(4,1)=v0*(k5*Ca);
end

function resp=der(t,point)
h=1e-7;
resp=(funode(t,point+h)-funode(t,point))/h;
end