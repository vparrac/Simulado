function M021
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
D=diag(diag(A));

end
