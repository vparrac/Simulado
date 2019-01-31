function Friccion_fsolve()
%Limpiar todo
clc
%Parametros
tic
global eD D hf l nu rho g
eD=0.0002;
D=0.3;
hf=8;
rho=950;
nu=2e-5;
l=100;
g=9.81;
%M. FSOLVE


sol=fsolve(@ecu,0.1);
disp('Factor fricción')
disp(sol(1))
disp('Velocidad en m/s')
disp((hf*(D/l)*(2*9.81/sol(1)))^0.5);
toc
end
%Funcion auxiliar
function y=ecu(x)
global eD D hf nu l

f=x(1);
%Calculos intermedios
v=(hf*(D/l)*(2*9.81/f))^0.5;
Re=D*v/nu;
%Ecuaciones
y(1)=(1.74-2*log10(2*eD+18.7/(Re*f^0.5)))^-2-f;
end