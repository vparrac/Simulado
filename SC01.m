function principal
global T P R
P=[1:40]*101.325;
T=500;
R=8314;
Vsol=fsolve(@(X) aux(X),P.*R./T)
plot(P,Vsol)
end

function resp=aux(V)
global T P R
Tc=425.2;
Pc=3797;
w=0.1931;
S=0.4850+1.55171*w-1.15613*(w^2);
alpha=(1+S*(1-(T/Tc)^0.5))^2;
a=(0.4278*R^2*Tc^2)/Pc;
b=(0.0867*R*Tc)/Pc;
for i=1:40
    resp(i)=((R*T)/(V(i)-b))-((a*alpha)/(V(i)*(V(i)+b)))-P(i);
end
end