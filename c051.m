function c051
clc
clear all

global hd tspan h
error=200;
tole=1e-6;
iter=0;
itermax=1000;
hd=1e-7;

a=0;
b=1;
N=10;
h=(b-a)/N;
tspan=a:h:b;

y1f=3;
y2f=-0.09697;
y30=3;
y10g=2;
y20g=2;
y0g=[y10g;y20g;y30];
yfconocido=[y1f;y2f];


while error>tole && iter<itermax
   sol=RK4(y0g);
   yf=sol(1:2,end);%[sol(1,end);sol(2,end)]
   fkx=(yf-yfconocido);
   
   J=jacobian(y0g);
   
   Dk=-inv(J)*fkx;
    
   y0gnew=y0g(1:2,:)+Dk;
   
   error=norm(Dk); %También la norma de fkx
   y0g(1:2,:)=y0gnew;
   
   iter=iter+1;
   
end
   
sol=RK4(y0g);
plot(tspan,sol(1,:),tspan,sol(2,:),tspan,sol(3,:));


end



function resp=ecs(t,var)
y1=var(1);
y2=var(2);
y3=var(3);
dy1=2*t;
dy2=-y1*y2;
dy3=3*y1-5*t*y2;
resp=[dy1;dy2;dy3];
end


function J=jacobian(y)
global hd
n=length(y);
for i=1:n-1
   vech=zeros(n,1);
   vech(i)=hd;   
   RungeK4h=RK4(y+vech);
   RungK4=RK4(y);
   der=(RungeK4h(1:2,:)-RungK4(1:2,:))/hd;
   J(:,i)=der(:,end);    
end
end



function resp=RK4(x0)
global h tspan 

y(:,1)=x0;

for i=1:length(tspan)-1
   yi=y(:,i);
   k1=ecs(tspan(i),yi);
   k2=ecs(tspan(i)+h/2,yi+(k1/2)*h);
   k3=ecs(tspan(i)+h/2,yi+(k2/2)*h);
   k4=ecs(tspan(i)+h,yi+k3*h);
   
   y(:,i+1)=yi+(k1+2*k2+2*k3+k4)*(h/6);
end


resp=y;
end