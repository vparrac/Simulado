function txy
clc
global P
P=101325/1000; %Presión en kPa

x0=linspace(0,1,100);% Composicion del agua

%Ciclo que llama al fsolve

for i=1:100
paraSol=fsolve(@(X) aux(X,x0(i)),x0);
y(i)=paraSol(3);
T(i)=paraSol(4);
end

%Graficar
figure (1)
plot(y,T,x0,T);
legend('Curva vapor','Curva líquido');
xlabel('Composicion del agua')
ylabel('Temperatura [C]');
title('Txy agua-Tolueno');

end


function resp=aux(var,x1)
global P
%Extraccion de variables a resolver
gamma1=var(1); %Gamma del agua
gamma2=var(2); %Gamma del Tolueno
y1=var(3);     %Composicion del vapor del agua
T=var(4);      %Temperatura

%Presiones de Sat
%Constantes del Agua
A1=16.3872; 
B1=3885.7;
C1=230.17;
%Constantes del tolueno
A2=13.932;
B2=3056.96;
C2=217.625;

%Parametro de Margules 2 ctes
A=3.5-0.006*(T+273);       %La temperatura se convierte para que quede en Kelvin
P1sat=exp(A1-(B1/(T+C1))); %Psat del agua
P2sat=exp(A2-(B2/(T+C2))); %Psat del tolueno
x2=1-x1;                   %Relación de composicion de agua y tolueno en liquido
y2=1-y1;                   %Relación de composicion de agua y tolueno en vapor

resp(1)=exp(A*(x2^2))-gamma1; %Ecuacion 1: residual del gamma del agua
resp(2)=exp(A*(x1^2))-gamma2; %Ecuacuin 2: residual del gamma del tolueno
resp(3)=gamma1*x1*(P1sat)-y1*P; %Ecuacion 3: residual de la presion del agua
resp(4)=gamma2*x2*(P2sat)-y2*P; %Ecuacion 4: residual de la presion del tolueno
end