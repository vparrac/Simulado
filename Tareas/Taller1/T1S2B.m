function T1S2B
%Preparacion del codigo
clc
%Declaracion de variables globales
global R V pi
R=3;
V=30;
pi=3.1416;

%Variables necesarias para el metodo de biseccion
tol=1e-6;
error=100;
iter=0;
itermax=1000; 

%Declarar un low y un up
u=0;
l=3;
%Evaluar la ecuacion en u y l
fu=ecs(u);
fl=ecs(l);

%Verificacion inicial para el metodo de biseccion
if fu*fl>0
   error('No se cumplen las condiciones para el método');    
end

%Metodo
tic
while error>tol && iter<itermax
   m=(u+l)/2;       %Hallar un m
   fl=ecs(l);       %Evaluar la ecuacion en l   
   fm=ecs(m);       %Evaluar la ecuacion en m
   
   if fl*fm>0       %Si fl y fm tiene el mismo signo
      l=m;          %Entonces el limite inferior sube a m
   else
       u=m;         %Sino, el superior baja a m
   end
   error=abs(fm);   %Encontrar el error
   iter=iter+1;     %Incrementar el numero de iteraciones
end
toc
fprintf('La altura del tanque es: ');
fprintf('%7.5f\n\n',m);
fprintf('El numero de iteraciones ');
fprintf('%2.0f\n\n ',iter);
disp(['El error final del metodo es: ',num2str(error),' ']);

end


function res=ecs(h)
global R V pi
res=pi*(h^2)*((3*R-h)/(3))-V;
end