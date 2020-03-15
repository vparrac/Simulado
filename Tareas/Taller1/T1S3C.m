function T1P2A
%Preparación
clc
clear all
%Declaración de la variable global a iterar
global e
e= 0.000001:0.000001:0.0001;
% M representa la inicialización que se le dará.
%El factor de fricción es muy pequeño por lo que se
%inicializará siempre en 0.003
m=1:1:100;
for i=1:size(m)
   m(i)=0.003; 
end

%Llamado del fsol
fsol=fsolve(@(X) friccion(X),m);

plot(e,fsol);
xlabel('Rugosidad en m')
ylabel('Factor de fricción')
title('Factor de friccion vs rugosidad')

end

function resp=friccion(f)
%Importar variables globales
global e
for i=1:100
%Calculos intermedios
Re=(1000*10*0.01)/(0.001);
%Ecuaciones en forma residual
    resp(i)=-2*log10(((e(i)/0.01)/3.7)+(2.51/(Re*sqrt(f(i)))))-(1/sqrt(f(i)));
end
end