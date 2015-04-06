clear all;
close all;
clc;

% Problem 2/3
% Constants
a=4; D=1; sigma=0.2; S=8;
% h values (Problem 3)
h1=1; h2=0.5; h3=0.1; h4=0.05; h5=0.01;
% x values (Problem 3)
x1=linspace(-4,4, (2*a/h1)+1); x2=linspace(-4,4, (2*a/h2)+1); x3=linspace(-4,4, (2*a/h3)+1); x4=linspace(-4,4, (2*a/h4)+1); x5=linspace(-4,4, (2*a/h5)+1);

phi01=exp((sigma-sqrt((sigma^2)-(4*S*D))).*x1)+exp((sigma+sqrt((sigma^2)-(4*S*D))).*x1);
phi02=exp((sigma-sqrt((sigma^2)-(4*S*D))).*x2)+exp((sigma+sqrt((sigma^2)-(4*S*D))).*x2);
phi03=exp((sigma-sqrt((sigma^2)-(4*S*D))).*x3)+exp((sigma+sqrt((sigma^2)-(4*S*D))).*x3);
phi04=exp((sigma-sqrt((sigma^2)-(4*S*D))).*x4)+exp((sigma+sqrt((sigma^2)-(4*S*D))).*x4);
phi05=exp((sigma-sqrt((sigma^2)-(4*S*D))).*x5)+exp((sigma+sqrt((sigma^2)-(4*S*D))).*x5);


[phi1,e1]=ne155_hw6_2(a,h1,D,S,sigma,phi01);
[phi2,e2]=ne155_hw6_2(a,h2,D,S,sigma,phi02);
[phi3,e3]=ne155_hw6_2(a,h3,D,S,sigma,phi03);
[phi4,e4]=ne155_hw6_2(a,h4,D,S,sigma,phi04);
[phi5,e5]=ne155_hw6_2(a,h5,D,S,sigma,phi05);


figure(1)
plot(x1, phi1, x1, phi01)
title('phi vs. x')
xlabel('x')
ylabel('phi')
legend('phi', 'phi actual')


% Problem 3 Figure
mesh=[(2*a)/h1, (2*a)/h2, (2*a)/h3, (2*a)/h4, (2*a)/h5];
error=[max(e1), max(e2), max(e3), max(e4), max(e5)];

figure(2)
plot(mesh, error)
title('error vs. total number of meshes')
xlabel('total number of meshes')
ylabel('max error')

% As the total number of meshes increases, the maximum error decreases. 


















