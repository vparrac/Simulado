function T1P2A
%Preparaci�n
clc
clear all
%Declaraci�n de la variable global a iterar
global v
v= 1:0.1:10;
% M representa la inicializaci�n que se le dar�.
%El factor de fricci�n es muy peque�o por lo que se
%inicializar� siempre en 0.003
m=1:0.1:10;

for i=1:size(m)
   m(i)=0.003; 
end

%Llamado del fsol
fsol=fsolve(@(X) friccion(X),m);

plot(v,fsol);
xlabel('Velocidad en m/s')
ylabel('Factor de fricci�n')
title('Factor de friccion vs velocidad')

end

function resp=friccion(f)
%Importar variables globales
global v
for i=1:91
%Calculos intermedios
Re=(1000*v(i)*0.01)/(0.001);
%Ecuaciones en forma residual
    resp(i)=-2*log10(((0.000002/0.01)/3.7)+(2.51/(Re*sqrt(f(i)))))-(1/sqrt(f(i)));
end
end