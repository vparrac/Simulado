function Friccion_SS()
%Limpiar todo
clc
%Parametros
global eD D hf l nu rho g
eD=0.0002;
D=0.3;
hf=8;
rho=950;
nu=2e-5;
l=100;
g=9.81;
%M. SS
cond=true;
tol=1e-6;
f0=0.05;

while cond==true
   fnew=ecuSS(f0);
   error=abs((fnew-f0)/(fnew));
   if error < tol
       break
   else
      f0=fnew; 
   end    
end

sol=f0;
disp('Factor fricción')
disp(sol)
disp('Velocidad en m/s')
disp((hf*(D/l)*(2*9.81/sol))^0.5);

end
%Funcion auxiliar
function y=ecu(x)
global eD D hf nu l

f=x(1);
%Calculos intermedios
v=(hf*(D/l)*(2*9.81/f))^0.5;
Re=D*v/nu;
%Ecuaciones
y(1)=(1.74-2*log10(2*eD+18.7/(Re*f^0.5)))^-2;  %%Ya no esta en residual
end