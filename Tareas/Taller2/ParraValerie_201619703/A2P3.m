function A2P3
clc
k21=0.2;
k12=0.1;
k31=0.8;
k13=0.5;
k24=0.4;
k42=0.6;
k34=0.4;
k43=0.25;
k54=0.35;
k45=0.11;
k65=0.21;
k56=0.45;
k76=0.15;
k67=0.23;
k75=0.41;
k57=0.69;

A0=1;
F0=1;

A=[
    k12+k13 -k21 -k31 0 0 0 0;    
    -k12 k21+k24 0 -k42 0 0 0;
    -k13 0 k31+k43 -k34 0 0 0;
    0 -k24 -k43 k42+k34 0 0 0;
    0 0 0 -k45 k65+k54+k57 -k56 -k75;
    0 0 0 0 -k65 k56+k76 -k67;
    1 1 1 1 1 1 1
    ];

b=[0;0;0;0;0;0;A0+F0];
M=[A,b];

%%Gauss-Seidel

tole=1e-6;
error=100;
D=diag(diag(A));
L=tril(A,-1);
U=triu(A,1);
Tg=inv(D+L)*U;
Cg=inv(D+L)*b;
x0=zeros(7,1);
while error>tole
    xnew=-Tg*x0+Cg;
    error= norm(xnew-x0);
    x0=xnew;
end

display('El vector encontrado con la función predeterminada de MatLab')
display(A\b);

display('El vector encontrado por Gauss-Seidel');
x0


n=length(b);
for i=1:n
    for j=[1:i-1,i+1:n]
        k=M(j,i)/M(i,i);
        M(j,:)=M(j,:)-k*M(i,:);
    end
end

display('El vector encontrado con Gauss');
b=M(:,n+1);
sol=b./diag(M)

A=[
    1 1 1 1 1 1 1;
    k12+k13 -k21 -k31 0 0 0 0;    
    -k12 k21+k24 0 -k42 0 0 0;
    -k13 0 k31+k43 -k34 0 0 0;
    0 -k24 -k43 k42+k34 0 0 0;
    0 0 0 -k45 k65+k54+k57 -k56 -k75;
    0 0 0 0 -k65 k56+k76 -k67    
    ];

b=[A0+F0;0;0;0;0;0;0];


G=[A,b];
[f,c]=size(G);
for i=1:c-1
    G(i,:)=G(i,:)/G(i,i);
    for j=i+1:f
        G(j,:)=G(j,:)-G(i,:)*G(j,i);
    end
end
for i=2:f
    ii=f-i+2;
    for j=1:ii-1
        jj=ii-j;
        G(jj,:)=G(jj,:)-G(ii,:)*G(jj,ii);
    end
end


display('El vector encontrado con Gauss-Jordan');

sol=G(:,end)

end

