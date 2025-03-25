%% Initialization and Parameter

% Clear command line window and workspace.
clear; clc;

% Tolerance of residual to stop multigrid iterations.
tol = 1.0e-10;

% Number of interior cells at finest level.
nL = 1023;

% Number of multigrid levels.
L = 9;

% Define the grid spacing at finest level.
hL = 1/(nL+1);

% Parameters for the multigrid solver.
p     = 1;
m1    = 3;
m2    = 3;
omega = 2/3;
kMax  = 100;

%% Setting up the multigrid solver

% Check mesh coarsening.
nl = nL;
for k = L-1:-1:0
  if mod(nl+1,2) ~= 0
    error('Coarsening error: L is too large.');
  end
  nl = (nl+1)/2-1;
end

% Give a random initial guess.
rng(1234);
u = rand(nL,1)-0.5;

% Construct the exact solution and necessary forcing vector.
x = linspace(hL,1.0-hL,nL)';
AL = gallery('tridiag',nL,-1.0/hL,2.0/hL,-1.0/hL);
uExact = exp(sin(2.0*pi*x))-1.0;
f = AL*uExact;

% Store parameters in a struct.
MGParam.nL = nL;
MGParam.L = L;
MGParam.p = p;
MGParam.m1 = m1;
MGParam.m2 = m2;
MGParam.omega = omega;
MGParam.tol = tol;
MGParam.kMax = kMax;

%% Calling the multigrid solver

tic
[u,errVals,kStop] = multiGridSolver(u,f,MGParam,uExact);
toc

finalErr = norm(uExact-u,2);
fprintf('\nFinal error in L2 sense = %g\n',finalErr);

%% Plotting the numerical solution and error plots

figure(1)
clf
plot(x,u,'k-','LineWidth',1.5)
xlabel('$x$','Interpreter','latex');
ylabel('$u_L(x)$','Interpreter','latex');
title('Finite Element Approximation of the Poisson Equation');
axis([0.0,1.0,-1.0,2.0])
printstr = strcat('Solution_nL_',num2str(nL),'.pdf');
exportgraphics(gca,printstr)

figure(2)
clf
semilogy(errVals(1:kStop,3),'bo','LineWidth',1.5)
hold on
semilogy(errVals(1:kStop,2),'rs','LineWidth',1.5)
hold on

% Estimate the rate of contraction.
rate = 1.0;
if kStop >= 4
  kv = kStop-3:kStop;
  le = log(errVals(kStop-3:kStop,3));
  p1 = polyfit(kv,le,1);
  rate = exp(p1(1));
  p1k = polyval(p1,1:kStop);
  semilogy(1:kStop,exp(p1k),'-k')
end

xlabel('$k$','Interpreter','latex');
title('Multigrid Iteration Errors','Interpreter','latex');
legend('$\left\|{\bf u}_L^{\rm E}-{\bf u}_L^k\right\|_L$',...
  '$\left\|{\bf u}_L^k-{\bf u}_L^{k-1} \right\|_L$',...
  'log-linear fit','Interpreter','latex');
text(1.5,128*tol,...
    strcat('$\gamma_{\rm comp} =\hspace{.1cm}$',...
  num2str(rate,'%10.5e')),'FontSize',14,'Interpreter','latex')
text(1.5,32*tol,strcat('$m_1 =\hspace{.1cm}$',...
  num2str(MGParam.m1)),'FontSize',14,'Interpreter','latex')
text(1.5,8*tol,strcat('$m_2 =\hspace{.1cm}$',...
  num2str(MGParam.m2)),'FontSize',14,'Interpreter','latex')
text(1.5,2*tol,strcat('$p =\hspace{.1cm}$',...
  num2str(MGParam.p)),'FontSize',14,'Interpreter','latex')
printstr = strcat('Err_nL_',num2str(nL),'_m1_',...
    num2str(MGParam.m1),'.pdf');
exportgraphics(gca,printstr)
hold off