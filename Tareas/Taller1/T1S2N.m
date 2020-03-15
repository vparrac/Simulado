function T1S2NR
%Preparacion del codigo
clc
%Declaracion de las variables globales
global R V pi
R=3; %Radio en m
V=30;%Volumen en m3
pi=3.1416;

%Variables
Vo=1;           %Valor iniciar para el método
tol=1e-6;       %Tolerancia del metodo
error=100;      %Error inicial
iter=0;         %Variable para no caer en un loop infinito
itermax=1000;   %Boundary
h=1e-3;         %Numero pequeño para la derivada
tic
while error>tol&&iter<itermax       %Mientras que deba iterar
    fun=ecs(Vo);                    %Evaluo la funcion en el punto
    derfun=(ecs(Vo+h)-ecs(Vo))/h;   %Evaluo la derivada de la funcion en el punto
    Vm=Vo-fun/derfun;               %Realizo el criterio
    error=abs((Vm-Vo)/Vm);          %Miro el error
    Vo=Vm;                          %Asigno el nuevo valor para la iteracion
    iter=iter+1;                    %Avanzo el numero de iteraciones
end
toc
%Salida de la consola
fprintf('La altura del tanque es: ');
fprintf('%7.5f\n\n',Vo);
fprintf('El numero de iteraciones');
fprintf('%2.0f\n\n',iter);
disp(['El error final del metodo es: ',num2str(error),' ']);

end


function res=ecs(h)
global R V pi
res=(pi*(h^2)*((3*R-h)/(3)))-V;
end