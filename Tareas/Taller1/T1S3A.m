function T1P2A
%Preparar el codigo
clc
clear all
%Declaracion global de d
global d
d= 0.001:0.00009:0.01;
%Declaracion vector inicial
u=zeros(1,101);
%Inicializacion
for i=1:101
   u(i)=0.001;
end
%Llamado al fsolve
fsol=fsolve(@(X) friccion(X),u);
%Graficar
plot(d,fsol);
xlabel('Diametro en m')
ylabel('Factor de fricción')
title('Factor de friccion vs diametro')

end

function resp=friccion(f)
global d
for i=1:101
%Calculos intermedios
Re=(1000*10*d(i))/(0.001);
%Ecuaciones
    resp(i)=-2*log10(((0.000002/d(i))/3.7)+(2.51/(Re*sqrt(f(i)))))-(1/sqrt(f(i)));
end
end