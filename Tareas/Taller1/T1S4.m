function T1S4
%Preparar el codigo
clc
%Declara variables globales
global P T Tc Pc R
P=2000000;
T=300;
Tc=190.6;
Pc=4599000;
R=8.314;

%Variables aux
h=1e-15;
tol=1e-6;
itermax=1000;
iter=1;
error=100;
Vo=0.00119744666;


%Metodo
while error>tol&&iter<itermax       %Mientras que deba iterar
    fun=ecs(Vo);                    %Evaluo la funcion en el punto
    derfun=(ecs(Vo+h)-ecs(Vo))/(h);   %Evaluo la derivada de la funcion en el punto
    Vm=Vo-(fun)/(derfun);               %Realizo el criterio
    error=abs((Vm-Vo)/Vm);           %Miro el error
    Vo=Vm;                          %Asigno el nuevo valor para la iteracion
    iter=iter+1;                    %Avanzo el numero de iteraciones
end

z=P*Vo/(R*T);
fprintf('Z es: ');
fprintf('%7.5f\n\n',z);

end 


function res=ecs(V)
global P T Tc Pc R
a=(27*(R^2)*(Tc^2))/(64*Pc);
b=(R*Tc)/(8*Pc);
%res=V^3-((b+(R*T/P))*V^2)+(a/P)*V-((a*b)/P);
res=(V^3)-(b+(R*T)/P)*V^2+((a/P)*V)-(a*b/P);
end