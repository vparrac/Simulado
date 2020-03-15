function A2P2
clc
tol=1e-6;
iter=0;
itermax=1000;
x0=ones(17,1);
fsolve(@fun,x0)
error=100;
while iter<itermax && error>tol
    funx=fun(x0);
    jacobiano=jacobian(x0);
    Dk=-inv(jacobiano)*funx;
    xnew=x0+Dk;
    error=norm(Dk);
    iter=iter+1;
    x0=xnew;
end



end
function resp=fun(var)
x3a=var(1);
x3o2=var(2);
x3n2=var(3);
x4a=var(4);
x4b=var(5);
x4co2=var(6);
x4n2=var(7);
x4h2o=var(8);
x6n2=var(9);
x6co2=var(10);
f1=var(11);
f2=var(12);
f3=var(13);
f4=var(14);
f5=var(15);
r1=var(16);
r2=var(17);

x1o2=0.75;
x1n2=0.25;
x2a=1;
x5a=0.5;
x5b=0.5;
x6h2o=0.75;
f6=100;


resp(1,1)=f1*x1o2-f3*x3o2;
resp(2,1)=f1*x1n2-f3*x3n2;
resp(3,1)=f2*x2a-f3*x3a;
resp(4,1)=f3*x3a-r1-r2-f4*x4a;
resp(5,1)=f3*x3o2-r1-r2;
resp(6,1)=f3*x3n2-f4*x4n2;
resp(7,1)=r1-f4*x4b;
resp(8,1)=r2-f4*x4co2;
resp(9,1)=r1+r2-f4*x4h2o;
resp(10,1)=x4a+x4b+x4co2+x4n2+x4h2o-1;
resp(11,1)=x3a+x3o2+x3n2-1;
resp(12,1)=f5*x5a-f4*x4a;
resp(13,1)=f5*x5b-f4*x4b;
resp(14,1)=f6*x6co2-f4*x4co2;
resp(15,1)=f6*x6n2-f4*x4n2;
resp(16,1)=f6*x6h2o-f4*x4h2o;
resp(17,1)=f5+f6-f4;


end

function j=jacobian(x)
n=length(x);
h=1e-7;
for i=1:n
   hvect=zeros(n,1);
   hvect(i)=h;
   der= (fun(x+hvect)-fun(x))/h;
   j(:,i)=der;
end
end