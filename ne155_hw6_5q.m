clear all;
clc; 
sigma=0.7;
nsigma=0.6;
D=1;
a=4;
S=8;
h=0.1;
n=(2*a)/h;
phi0=zeros(n,1);
k0=1;
tol1= 10^(-4);
tol2=tol1;
alpha=0; beta=0; p=0; q=(sigma/D); r=-S/D;

A=zeros(n, n);

for i = 1: n-1
    A(i, i) = (2+h^2)*q;
    A(i, i+1) = ((h/2)*p)-1;
    A(i+1, i) = ((-h/2)*p)-1;
end
A(n, n)=(2+h^2)*q;

[phi]=ne155_hw6_5b(a,h,D,S,sigma,phi0);

Q0=zeros(n,1);
Q0(1)=nsigma*phi(1);
for i =2:n-1
    Q0(i)=nsigma*phi(i);
end

[phin,k]=ne155_hw6_5(n, h, A, k0, Q0, nsigma, tol1, tol2);
eigen=max(abs(k))

x=linspace(-4,4,length(phin));
plot(x,phin)
title('phi vs. x')
xlabel('x')
ylabel('phi')