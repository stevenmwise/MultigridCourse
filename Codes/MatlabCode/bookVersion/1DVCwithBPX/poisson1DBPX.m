%% Initialization and Parameters

% Clear command line window and workspace.
clear; clc;

% Tolerance of residual to stop PCG iterations.
tol = 1.0e-10;

% Number of interior points at finest level.
nL = 1023;

% Number of multigrid levels.
L = 9;

% Define the domain of the equation.
xLower = 0.0;
xUpper = 1.0;

% Grid spacing at finest level.
hL = 1/(nL+1);

% Check that nL+1 is divisible by 2^(L-1).
nl = nL;
for k = L-1:-1:0
  if mod(nl+1,2) ~= 0
    error('Coarsening error: L is too large.');
  end
  nl = (nl+1)/2-1;
end

% Initialize our approximate solution (e.g. random).
rng(1234);
u0 = rand(nL,1)-0.5;
px = xUpper-xLower;

% Build exact solution and define f accordingly.
x = linspace(xLower+hL,xUpper-hL,nL)';
AL = gallery('tridiag',nL,-1/hL,2/hL,-1/hL);
uExact = exp(sin(2*pi*x/px))-1;
f = AL*uExact;

% MGParam structure.
MGParam.nL = nL;
MGParam.L = L;
MGParam.kMax = 100;

% Preconditioner function handle.
precond = @(r) BPXPreconditioner(r,MGParam);

%% Calling the PCG solver

% Solve using PCG with BPX preconditioner.
tic;
[u,flag,relres,iter,resvec] = myPCG(AL,f,tol,MGParam.kMax, ...
    precond,u0);
toc;

% We can do a quick final L2 error check.
finalErr = norm(uExact-u,2);
fprintf('\nFinal error in L2 sense = %g\n',finalErr);

%% Plotting

figure(1);
clf;
plot(x,u,'k-','LineWidth',1.5);
xlabel('$x$','Interpreter','latex');
ylabel('$u_L(x)$','Interpreter','latex');
title('Numerical Solution with BPX Preconditioned PCG');
axis([xLower,xUpper,-1.0,2.0]);

figure(2);
clf;
semilogy(resvec,'ro','LineWidth',1.5);
xlabel('$k$','Interpreter','latex');
title('PCG-BPX Iteration Residuals','Interpreter','latex');

% Calculate convergence rate from last few iterations.
rate = 1.0;
if length(resvec) >= 4
  kv = length(resvec)-3:length(resvec);
  le = log(resvec(length(resvec)-3:length(resvec)));
  p1 = polyfit(kv,le,1);
  rate = exp(p1(1));
  p1k = polyval(p1,1:length(resvec));
  hold on
  semilogy(1:length(resvec),exp(p1k),'b','LineWidth',2)
end

ylabel('$\left\|{\bf r}_L^k\right\|_L$','Interpreter','latex');
text(1.5,256*tol, ...
    strcat('$\gamma_{\rm comp} =\hspace{.1cm}$', ...
    num2str(rate,'%10.5e')),'FontSize',14,'Interpreter','latex')
text(1.5,32*tol,strcat('$L =\hspace{.1cm}$', ...
    num2str(L)),'FontSize',14,'Interpreter','latex')
text(1.5,8*tol,strcat('$n_L =\hspace{.1cm}$', ...
    num2str(nL)),'FontSize',14,'Interpreter','latex')

printstr = strcat('Err_nL_', ...
    num2str(nL),'_L_',num2str(L),'.pdf');
exportgraphics(gca,printstr)
hold off