function [phin,kn]=ne155_hw6_5(n, h, A, k0, Q0, nsigma, tol1, tol2)

phin=inv(A).*(1/k0)*Q0;
Qn=nsigma.*(phin);
Q01=nsigma*phin(1);

summ = 0;
summ2 = 0;

for i = 2:n-1
    summ=summ+(Qn(i)*h/2);
end
for i = 2:n-1
    summ2=summ2+(Q0(i)*h/2);
end

kn=k0*((Q01*(h/2))+summ)/((Q0*(h/2))+summ2);

k0=kn;
Q0=Qn;

while tol1<abs((kn-k0)/kn) && tol2<abs((phin-phi0)/phin)
    phin=inv(A).*(1/k0)*Q0;
    Qn=nsigma.*(phin);
    Q01=nsigma*phin(1);
    
    summ = 0;
    summ2 = 0;
    
    for i = 2:n-1
        summ=summ+(Qn(i)*h/2);
    end
    for i = 2:n-1
        summ2=summ2+(Q0(i)*h/2);
    end
    
    kn=k0*((Q01*(h/2))+summ)/((Q0*(h/2))+summ2);
    k0=kn;
    Q0=Q;
end

