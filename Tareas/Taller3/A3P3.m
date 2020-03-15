function A3P3
clc
clear all
close all

global tspan hd
error=200;
tole=1e-6;
iter=0;
itermax=1000;
hd=1e-7;

m4f=1;
m1f=1;
m4g=0;
m1g=0;
m30=0;
m20=0;

y0g=[m1g;m4g;m20;m30];

yfconocido=[m1f;m4f];
a=0.000001;
b=1;
N=100;
h=(b-a)/N;
tspan=a:h:b;
hd=1e-7;


while error>tole && iter<itermax
    [xit,yit]=ode45(@fun,tspan,y0g);
    yf=yit(end,1:2)';
    fxy=(yf-yfconocido);
    J=jacobian(y0g);
    Dk=-inv(J)*fxy;
    y0gnew=y0g(1:2,:)+Dk;
    error=norm(Dk);
    y0g(1:2,:)=y0gnew;
    iter=iter+1;    
end
[titer,soliter]=ode45(@fun,tspan,y0g);
figure (1)
plot(titer,soliter(:,1),titer,soliter(:,2));
legend('y','T');
title('T,y vs x')
xlabel('Coordenada x')
figure (2)
plot(titer,soliter(:,3),titer,soliter(:,4));
legend("y'","T'");
title("T',y' vs x");
xlabel('Coordenada x')
end

function resp=fun(x,var)
m1=var(1);
m4=var(2);
m2=var(3);
m3=var(4);
phi=0.38;
gamma=25;
betha=0.4;
resp(1,1)=m2;
resp(2,1)=m3;
resp(3,1)=(phi^2)*m1*exp(gamma*(1-1/m4))-(2/x)*m2;
resp(4,1)=-betha*(phi^2)*m1*exp(gamma*(1-1/m4))-(2/x)*m3;
end


function J=jacobian(x)
global hd tspan
n=length(x);
for i=1:n-2
    vech=zeros(n,1);
    vech(i)=hd;
    [x1,funh]=ode45(@fun,tspan,x+vech);
    [x2,funx]=ode45(@fun,tspan,x);
    der=(funh(:,1:2)-funx(:,1:2))/hd;
    %El jacobiano sólo es de las condiciones 
    J(:,i)=der(end,:)';    
end

end