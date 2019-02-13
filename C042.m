function dae()
clc
M=[1 0 0
   0 1 0
   0 0 0];
y0=[1;0;1e-3]; %Orden de mis variables

tspan=[0 4*logspace(-6,6)];%Linspace en escala logaritmica
options=odeset('Mass',M); %Importante el nombre, porque cambia el valor por parámetro
[t,y]=ode15s(@ecs,tspan,y0,options);

y(:,2)=1e4*(y(:,2));
semilogx(t,y);
end


function resp=ecs(t,var)

y1=var(1);
y2=var(2);
y3=var(3);

dy1=-0.04*y1+1e4*y2*y3;
dy2=0.04*y1-1e4*y2*y3-3e7*y2^2;
ec1=y1+y2+y3-1;

resp=[dy1;dy2;ec1];
end