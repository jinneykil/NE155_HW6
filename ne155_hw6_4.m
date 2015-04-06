%% h = 1cm
clear all;
clc;
a=4; D=1; sigma=0.2; S=8; h=1; w=1.2;

alpha=0; beta=0; p=0; q=(sigma/D); r=-S/D;

n=((2*a)/h)+1;
A=zeros(n-1, n-1);

for i = 1: n-2
    A(i, i) = (2+h^2)*q;
    A(i, i+1) = ((h/2)*p)-1;
    A(i+1, i) = ((-h/2)*p)-1;
end
A(n-1, n-1)=(2+h^2)*q;

b=zeros(n-1,1);
b(1)=((h^2)*r)+(-alpha(((h/2)*p)+1));
for i = 2:n-2
    b(i)=(h^2)*r;
end
b(n-1)=((h^2)*r)+(-(((h/2)*p)+1)*beta);

x0=zeros(n-1,1); tol1 = 10^(-3); tol2 = 10^(-5);

[phiJh1_1,kJh1_1]=Jacobi2(A, b, x0, tol1);
[phiGh1_1,kGh1_1]=GaussSeidel2(A, b, x0, tol1);
[phiSh1_1, kSh1_1]=SOR2(A, b, x0, tol1, w);

[phiJh1_2,kJh1_2]=Jacobi2(A, b, x0, tol2);
[phiGh1_2,kGh1_2]=GaussSeidel2(A, b, x0, tol2);
[phiSh1_2, kSh1_2]=SOR2(A, b, x0, tol2, w);

Kh1=[kJh1_1, kGh1_1, kSh1_1; kJh1_2, kGh1_2, kSh1_2];

%% h = 0.5 cm

a=4; D=1; sigma=0.2; S=8; h=0.5; w=1.2;

alpha=0; beta=0; p=0; q=(sigma/D); r=-S/D;

n=((2*a)/h)+1;
A=zeros(n-1, n-1);

for i = 1: n-2
    A(i, i) = (2+h^2)*q;
    A(i, i+1) = ((h/2)*p)-1;
    A(i+1, i) = ((-h/2)*p)-1;
end
A(n-1, n-1)=(2+h^2)*q;

b=zeros(n-1,1);
b(1)=((h^2)*r)+(-alpha(((h/2)*p)+1));
for i = 2:n-2
    b(i)=(h^2)*r;
end
b(n-1)=((h^2)*r)+(-(((h/2)*p)+1)*beta);

x0=zeros(n-1,1); tol1 = 10^(-3); tol2 = 10^(-5);

[phiJh2_1,kJh2_1]=Jacobi2(A, b, x0, tol1);
[phiGh2_1,kGh2_1]=GaussSeidel2(A, b, x0, tol1);
[phiSh2_1, kSh2_1]=SOR2(A, b, x0, tol1, w);

[phiJh2_2,kJh2_2]=Jacobi2(A, b, x0, tol2);
[phiGh2_2,kGh2_2]=GaussSeidel2(A, b, x0, tol2);
[phiSh2_2, kSh2_2]=SOR2(A, b, x0, tol2, w);

Kh2=[kJh2_1, kGh2_1, kSh2_1; kJh2_2, kGh2_2, kSh2_2];

%% h = 0.1 cm

a=4; D=1; sigma=0.2; S=8; h=0.1; w=1.2;

alpha=0; beta=0; p=0; q=(sigma/D); r=-S/D;

n=((2*a)/h)+1;
A=zeros(n-1, n-1);

for i = 1: n-2
    A(i, i) = (2+h^2)*q;
    A(i, i+1) = ((h/2)*p)-1;
    A(i+1, i) = ((-h/2)*p)-1;
end
A(n-1, n-1)=(2+h^2)*q;

b=zeros(n-1,1);
b(1)=((h^2)*r)+(-alpha(((h/2)*p)+1));
for i = 2:n-2
    b(i)=(h^2)*r;
end
b(n-1)=((h^2)*r)+(-(((h/2)*p)+1)*beta);

x0=zeros(n-1,1); tol1 = 10^(-3); tol2 = 10^(-5);

[phiJh3_1,kJh3_1]=Jacobi2(A, b, x0, tol1);
[phiGh3_1,kGh3_1]=GaussSeidel2(A, b, x0, tol1);
[phiSh3_1, kSh3_1]=SOR2(A, b, x0, tol1, w);

[phiJh3_2,kJh3_2]=Jacobi2(A, b, x0, tol2);
[phiGh3_2,kGh3_2]=GaussSeidel2(A, b, x0, tol2);
[phiSh3_2, kSh3_2]=SOR2(A, b, x0, tol2, w);

Kh3=[kJh3_1, kGh3_1, kSh3_1; kJh3_2, kGh3_2, kSh3_2];

%% h = 0.05 cm

a=4; D=1; sigma=0.2; S=8; h=0.05; w=1.2;

alpha=0; beta=0; p=0; q=(sigma/D); r=-S/D;

n=((2*a)/h)+1;
A=zeros(n-1, n-1);

for i = 1: n-2
    A(i, i) = (2+h^2)*q;
    A(i, i+1) = ((h/2)*p)-1;
    A(i+1, i) = ((-h/2)*p)-1;
end
A(n-1, n-1)=(2+h^2)*q;

b=zeros(n-1,1);
b(1)=((h^2)*r)+(-alpha(((h/2)*p)+1));
for i = 2:n-2
    b(i)=(h^2)*r;
end
b(n-1)=((h^2)*r)+(-(((h/2)*p)+1)*beta);

x0=zeros(n-1,1); tol1 = 10^(-3); tol2 = 10^(-5);

[phiJh4_1,kJh4_1]=Jacobi2(A, b, x0, tol1);
[phiGh4_1,kGh4_1]=GaussSeidel2(A, b, x0, tol1);
[phiSh4_1,kSh4_1]=SOR2(A, b, x0, tol1, w);

[phiJh4_2,kJh4_2]=Jacobi2(A, b, x0, tol2);
[phiGh4_2,kGh4_2]=GaussSeidel2(A, b, x0, tol2);
[phiSh4_2,kSh4_2]=SOR2(A, b, x0, tol2, w);

Kh4=[kJh4_1, kGh4_1, kSh4_1; kJh4_2, kGh4_2, kSh4_2];

%% h = 0.01 cm

a=4; D=1; sigma=0.2; S=8; h=0.01; w=1.2;

alpha=0; beta=0; p=0; q=(sigma/D); r=-S/D;

n=((2*a)/h)+1;
A=zeros(n-1, n-1);

for i = 1: n-2
    A(i, i) = (2+h^2)*q;
    A(i, i+1) = ((h/2)*p)-1;
    A(i+1, i) = ((-h/2)*p)-1;
end
A(n-1, n-1)=(2+h^2)*q;

b=zeros(n-1,1);
b(1)=((h^2)*r)+(-alpha(((h/2)*p)+1));
for i = 2:n-2
    b(i)=(h^2)*r;
end
b(n-1)=((h^2)*r)+(-(((h/2)*p)+1)*beta);

x0=zeros(n-1,1); tol1 = 10^(-3); tol2 = 10^(-5);

[phiJh5_1,kJh5_1]=Jacobi2(A, b, x0, tol1);
[phiGh5_1,kGh5_1]=GaussSeidel2(A, b, x0, tol1);
[phiSh5_1,kSh5_1]=SOR2(A, b, x0, tol1, w);

[phiJh5_2,kJh5_2]=Jacobi2(A, b, x0, tol2);
[phiGh5_2,kGh5_2]=GaussSeidel2(A, b, x0, tol2);
[phiSh5_2,kSh5_2]=SOR2(A, b, x0, tol2, w);

Kh5=[kJh5_1, kGh5_1, kSh5_1; kJh5_2, kGh5_2, kSh5_2];


% As total number of meshes increases, the number of iterations needed
% decreases. 




