clear

nL = 255;
hL = 1/(nL+1);
L = 7;
p = 1;
m1 = 3;
m2 = 3;
omega = 2/3;
tol = 1e-09;
kMax = 20;

x = linspace(hL,1.0-hL,nL)';
A1 = gallery('tridiag',nL,-1.0/hL,2.0/hL,-1.0/hL);
uExact = exp(sin(3.0*pi*x))-1.0;
f = A1*uExact;

MGParam.nL = nL;
MGParam.L = L;
MGParam.p = p;
MGParam.m1 = m1;
MGParam.m2 = m2;
MGParam.omega = omega;

[u,rate] = multiGrid(f,MGParam,kMax,tol,uExact);

clf
plot(x,u,'k-','LineWidth',1.5)
xlabel('$x$','Interpreter','latex');
ylabel('$u_L(x)$','Interpreter','latex');
title('Finite Element Approximation of the Poisson Equation');
axis([0.0,1.0,-1.0,2.0])
printstr = strcat('Solution_nL_',num2str(nL),'.pdf');
exportgraphics(gca, printstr)