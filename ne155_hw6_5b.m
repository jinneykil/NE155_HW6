function [phi]=ne155_hw6_5b(a,h,D,S,sigma,phi0)

n=(2*a)/h;
phi=zeros(1,n);
for i=2:n-1
    phi(i)=(((h^2)*(S/D))+phi(i+1)+phi(i-1))/((2+((h^2)*sigma/D)));
end
end
