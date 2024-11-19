clear

n1 = 255;
h1 = 1/(n1+1);
m1 = 3;
omega = 2/3;
tol = 1e-09;
kMax = 20;

%
% Seed a random intial guess.
SEED = 1234;
rng(SEED);
u = rand(n1,1)-0.5;

%
% Construct the exact solution and necessary forcing vector.
x = linspace(h1,1.0-h1,n1)';
A1 = gallery('tridiag',n1,-1.0/h1,2.0/h1,-1.0/h1);
uExact = exp(sin(3.0*pi*x))-1.0;
f = A1*uExact;

[u,errVals,kStop] = twoGridSolver(u,f,n1,m1,omega,kMax,tol,uExact);

figure(1)
clf

plot(x,u,'k-','LineWidth',1.5)
xlabel('$x$','Interpreter','latex');
ylabel('$u_1(x)$','Interpreter','latex');
title('Finite Element Approximation of the Poisson Equation');
axis([0.0,1.0,-1.0,2.0])
printstr = strcat('Solution_n1_',num2str(n1),'.pdf');
exportgraphics(gca, printstr)

figure(2)
clf

semilogy(errVals(1:kStop,3),'bo','LineWidth',1.5)
hold on
semilogy(errVals(1:kStop,2),'rs','LineWidth',1.5)
hold on

%
% Estimate the rate of contraction:
rate = 1.0;
if kStop >= 4
  kv = [kStop-3:kStop];
  le = log(errVals(kStop-3:kStop,3));
  p1 = polyfit(kv,le,1);
  rate = exp(p1(1));
  p1k = polyval(p1,[1:kStop]);
  semilogy([1:kStop],exp(p1k),'-k')
end

xlabel('$k$','Interpreter','latex');
title('Multigrid Iteration Errors', ...
  'Interpreter','latex');
legend('$\left\|{\bf u}_1^{\rm E}-{\bf u}_1^k\right\|_1$', ...
  '$\left\|{\bf u}_1^k-{\bf u}_1^{k-1} \right\|_1$', ...
  'log-linear fit','Interpreter','latex');
text(1.5,128*tol,strcat( ...
  '$\gamma_{\rm comp} =\hspace{.1cm}$', ...
  num2str(rate,'%10.5e')),'FontSize',14, ...
  'Interpreter','latex')
text(1.5,32*tol,strcat('$m_1 =\hspace{.1cm}$', ...
  num2str(m1)),'FontSize',14, ...
  'Interpreter','latex')
text(1.5,8*tol,strcat('$\omega =\hspace{.1cm}$', ...
  num2str(omega)),'FontSize',14, ...
  'Interpreter','latex')
printstr = strcat('Err_n1_',num2str(n1), ...
  '_m1_',num2str(m1), '_omega_', ...
  num2str(omega),'.pdf');
exportgraphics(gca, printstr)
hold off