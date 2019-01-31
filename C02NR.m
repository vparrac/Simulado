%%Complementaria 2
function C02
clc
global T P R
T=450;
P=56;
R=0.08206;
Vo=(R*T)/P;

%Mètodo como tal
tol=1e-6;
error=100;
iter=0;
itermax=1000;
h=1e-7;
while error>tol && iter<itermax
    fun=ecs(Vo);
    %derfun=dersyms(Vo);
    derfun=(ecs(Vo+h)-ecs(Vo))/h;
    V=Vo-fun/derfun;
    error=abs((V-Vo)/V);
    Vo=V;
    iter=iter+1;
end

fprintf('El volumen especifico del gas amonio es;');
fprintf('%7.5f\n\n',Vo);
fprintf('El numero de iteraciones');
fprintf('%2.0f\n\n',iter);
disp(['El error final del metodo es: ',num2str(error),' ']);
end


function der=dersyms(var)
V=var;
global P R T
%Definición de la derivada simbolica
syms Vs;
dersym=diff((R*T)/(Vs-b)-(a)/(T^(1/2)*Vs*(Vs+b))-P,Vs);
der=double(subs(dersym,Vs,V));
end

function resp= ecs(V)
global T P R
b=0.02590;
a=4.2527;
resp=((R*T)/(V-b))-(a/(sqrt(T)*V*(V+b)))-P;
end

