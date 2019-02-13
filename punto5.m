function T1P2A
clc
clear all
rho=998.2;
d=0.05;
miu=1002e-6;
vsol=fsolve(@(v) fun(v),[10;10])
(rho*d*vsol(1,1))/(miu)
end

function resp=fun(var)
v1=var(1);
v2=var(2);
L1=60;
D1=0.05;
L2=55;
D2=0.04;
g=9.8;
Kl=1.5;
Q=0.036;
f1=fsolve(@(f) friccion(v1,f,D1),0.026);
f2=fsolve(@(f) friccion(v2,f,D2),0.026);
f1
f2
resp(1,1)=f2*(L2/D2)*(v2^2)/(2*g)+Kl*(v2^2)/(2*g)-((v1^2)/(2*g))*(L1/D1)*f1;
resp(2,1)=Q-v1*pi*(D1^2)/4-v2*pi*(D2^2)/4;
end

function aux=friccion(v,f,d)
e=0.00015;
rho=998.2;
miu=1002e-6;
Re=(rho*d*v)/(miu);
aux=-2*log10(((e/d)/3.7)+(2.51/(Re*sqrt(f))))-(1/sqrt(f));
end