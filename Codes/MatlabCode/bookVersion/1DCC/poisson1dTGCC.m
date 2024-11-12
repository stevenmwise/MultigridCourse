clear

n1 = 1024;
h1 = 1/n1;
m1 = 3;
omega = 2/3;
tol = 1e-09;
kMax = 40;

x = linspace(h1/2,1.0-h1/2,n1)';
A1 = gallery('tridiag',n1,-1.0,2.0,-1.0);
A1( 1, 1) = 3.0;
A1(n1,n1) = 3.0;
uExact = exp(sin(3.0*pi*x))-1.0;
%f = A1*uExact/(h1*h1);
f = A1*uExact;

[u,rate] = twoGridCC(f,n1,m1,omega,kMax,tol,uExact);

clf
plot(x,u(2:n1+1),'k-','LineWidth',1.5)
xlabel('$x$','Interpreter','latex');
ylabel('$u_1(x)$','Interpreter','latex');
title('Cell-Centered FD Approximation of the Poisson Equation');
axis([0.0,1.0,-1.0,2.0])
printstr = strcat('Solution_n1_',num2str(n1),'.pdf');
exportgraphics(gca, printstr)