function C03

tol= 1e-6;
error=200;
iter=0;
itermax=2000;
x0=[0.7;2000];

while iter<itermax&&error>tol
   fuxk=ecs(x0);
   J= Jacobian(x0);
   Dk=-inv(J)*fuxk;
   xnew=x0+Dk;
   error=norm(Dk);
   iter=iter+1;
   x0=xnew;
end

display(xnew(1));
display(xnew(2));

end 




function resp = ecs(var)
 
X = var(1);
T = var(2);
 
% Parametros
nO=0.4;                               % Moles de Oxigeno
P=1;                                  % Presion [atm]
R=8.314;                              % Constante G.I.
To=25+273.15;                         % Temperatura de referencia [K]
DGo=-257e3;                           % Delta G estandar [J/mol CO]
DHo=-283e3;                           % Delta H estandar [J/mol CO]
A=[26.16 25.66 28.67 26.37];          % Constante A para Cp [CO O2 CO2 N2]
B=[8.75 12.52 35.72 7.61]*10^-3;      % Constante B para Cp [CO O2 CO2 N2]
C=[-1.92 -3.37 -10.39 -1.44]*10^-6;   % Constante C para Cp [CO O2 CO2 N2]
 
% Calculos Intermedios
n=[1-X,nO-X/2,X,(79/21)*nO];                     % Moles_i [CO O2 CO2 N2]
Am=n*A';                                         % Cte A de mezcla para Cp
Bm=n*B';                                         % Cte B de mezcla para Cp
Cm=n*C';                                         % Cte C de mezcla para Cp
CpdT=Am*(T-To)+Bm*(T^2-To^2)/2+Cm*(T^3-To^3)/3;  % Cp de la mezcla
 
v=[-1,-0.5,1];                                   % Coef. Esteq. [CO O2 CO2]
Ae=v*A(1:3)';                                    % Cte A de rxn para Cp
Be=v*B(1:3)';                                    % Cte B de rxn para Cp
Ce=v*C(1:3)';                                    % Cte C de rxn para Cp
CpdT_rxn=Ae*(T-To)+Be*(T^2-To^2)/2+Ce*(T^3-To^3)/3;         % CpdT rxn
CpTdT_rxn=Ae*(log(T)-log(To))+Be*(T-To)+Ce*(T^2-To^2)/2;    % (Cp/T)dT
DGT=(DHo/T)+(CpdT_rxn/T)+(DGo-DHo)/To-CpTdT_rxn;            % dG(T)/T
Ka=exp(-DGT/R);                                             % Keq en Termo
Kp=((X)/(1-X))*(((1+(100/21)*nO-X/2)/(P*(nO-X/2)))^(1/2));  % Keq en P
 
%Funci?n objetivo a resolver
 
resp(1,1)=Kp-Ka;
resp(2,1)=X*DHo+CpdT;
 
end

function J=Jacobian(x)
n=length(x);
h=1e-7;
 for i=1:n
    vech=zeros(n,1);
     vech(i)=h;     der=(ecs(x+vech)-ecs(x))/h;
     J(:,i)=der;

 end
 
end