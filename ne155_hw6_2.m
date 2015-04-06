function [phi,e]=ne155_hw6_2(a,h,D,S,sigma,phi0)

n=-a:h:a;
e=zeros(1, length(n));
phi=zeros(1,length(n));
for i=2:length(n)-1
    phi(i)=(((h^2)*(S/D))+phi(i+1)+phi(i-1))/((2+((h^2)*sigma/D)));
end
for i = 1: length(n)
    e(i)=abs((phi(i)-phi0(i))/phi0(i));
end
end
