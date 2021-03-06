%%Complementaria 2, Quasi  Newton
function C02
clc
global T P R
T=450;
P=56;
R=0.08206;


%Mètodo como tal
tol=1e-6;
err=100;
iter=0;
itermax=1000;
h=1e-7;
l=0.04;
fl=ecs(l);
u=1;
fu=ecs(u);

if fu*fl>0
    error('No se cumple la condicion del metodo');
end

while err>tol && iter<itermax
    m=(u+l)/2;
    fl=ecs(l);
    fm=ecs(m);
    
    if fl*fm>0
       l=m;
    else
        u=m;
    end
    
    err=abs(fm);
    iter=iter+1;
end

fprintf('El volumen especifico del gas amonio es ');
fprintf('%7.5f\n\n',m);
fprintf('El numero de iteraciones ');
fprintf('%2.0f\n\n',iter);

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
