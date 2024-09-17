clear

n1 = 255;
h1 = 1/(n1+1);
m1 = 20;
omega = 2/3;
tol = 1e-09;
kMax = 20;

x = linspace(h1,1.0-h1,n1)';
A1 = gallery('tridiag',n1,-1.0/h1,2.0/h1,-1.0/h1);
uExact = exp(sin(3.0*pi*x))-1.0;
f = A1*uExact;

[u,rate] = twoGrid(f,n1,m1,omega,kMax,tol,uExact);

clf
plot(x,u,'k-','LineWidth',1.5)
xlabel('$x$','Interpreter','latex');
ylabel('$u_1(x)$','Interpreter','latex');
title('Finite Element Approximation of the Poisson Equation');
axis([0.0,1.0,-1.0,2.0])
printstr = strcat('Solution_n1_',num2str(n1),'.pdf');
exportgraphics(gca, printstr)