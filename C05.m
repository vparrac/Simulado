function C05
%dy/dt, dy1/t, dy2/dt ..... n, tengo y0 para y1 pero n-1 condiciones de frontera. 
%Entonces quiero mover mis condiciones hacia atrás. Combinamos Newton
%Raphson con RK, el método diferencial con el método iterativo
%Te inventas una inicialización del valor que no conozco.--- Te lo inventas
%con cierto sentido.
%Jacobiano reducido de las cosas que estamos buscando

error=200;
tole=1e-6;
iter=0;
itermax=1000;
hd=1e-7;

a=0;
b=1;
N=10;
h=(b-a)/N;
tspan=a:h:b;


y20=-1;
y1f=3;
y10g=0;


y0g=[y10g;y20];



while error>tole && iter<itermax
   [t,sol]=ode45(@ecs,tspan,y0g);
   yf=sol(end,1);
   y0gh= y0g+[hd;0];
   [t,solh]=ode45(@ecs,tspan,y0gh);
   yfh=solh(end,1);   
   der=(yfh-yf)/hd;
   dif=yf-y1f;   
   y10gnew=y10g-dif/der;
   error=abs(dif);
   y10g=y10gnew;
   y0g=[y10g;y20];   
   iter=iter+1;   
end


 [t,sol]=ode45(@ecs,tspan,y0g);
 plot(t,sol(:,1),t,sol(:,2));

end

function resp=ecs(t,var)
%Esta ecución es normal, no tiene nada de los métodos númericos
y1=var(1);
y2=var(2);
dy1=2*t;
dy2=-y1*y2;
resp=[dy1;dy2];
end