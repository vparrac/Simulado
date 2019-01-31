function M02
clc
clear all
close all
%parametros
error=100;
tol=1e-6;
x0=[1;1;1]
A=[120 -20 0 
-80 80 0
-40 -60 120];
b=[400;0;200];
z=A\b;
%Descomposicion
D=diag(diag(A));%diag solo me da el vector diagonal, diagonal de la diagonal me da una matriz de diagonal dicha
L=tril(A,-1);%tril(A,0)=Coge la parte inferior incluyendo la diagonal, -1 inferior sin la diagonal y 1 superior sin diagonal
U=triu(A,1);
%Gauss-Seidel
Tg=inv(D+L)*U;
Cg=inv(D+L)*b;

while error>tol
   xnew=-Tg*x0+Cg; 
   error= norm(xnew-x0);
   x0=xnew;
end

disp('Resultado');
disp(xnew);
disp('Resultado 2')
disp(z);
end


