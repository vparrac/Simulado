function RK4
clc
a=0;
b=1;

N=1000;

h=(b-a)/N;

tspan=a:h:b;

y(1,1)=0;
y(2,1)=0;

for i=1:length(tspan)-1
    k1=ecs(tspan(i),y(:,i));
    k2=ecs(tspan(i)+h/2,y(:,i)+(k1/2)*h);
    k3=ecs(tspan(i)+h/2,y(:,i)+(k2/2)*h);
    k4=ecs(tspan(i)+h,y(:,i)+k3*h);
    
    y(:,i+1)=y(:,i)+(k1+2*k2+2*k3+k4)*(h/6);
end

[t,sol]=ode45(@ecs,tspan,[0;0]);

plot(tspan,y(1,:),tspan,y(2,:),t,sol(:,1),t,sol(:,2));

end
function resp=ecs(t,var)
x1=var(1);
x2=var(2);



dydt1=-20*x1+10*x2+100;
dydt2=10*x1-20*x2;

resp=[dydt1;dydt2];

end
