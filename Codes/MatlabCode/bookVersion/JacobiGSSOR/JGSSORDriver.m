clear

n = 31;
h = 1/(n+1);
h2 = h*h;
kMax = 30000;
tol = 1e-07;

errVals = zeros(kMax,3);
kStop = zeros(3,1);
%
% Store variables in a parameter structure.
iterParam.n = n;
iterParam.h = h;
iterParam.h2 = h2;
iterParam.kMax = kMax;
iterParam.tol = tol;
%
%Construct the right-hand-side vector, f.
x = linspace(h,1.0-h,n)';
A = gallery('tridiag',n,-1.0/h2,2.0/h2,-1.0/h2);
wExact = exp(sin(3.0*pi*x))-1.0;
f = A*wExact;
%
% Initialize w.
SEED = 1234;
rng(SEED);
w = rand(n,1)-0.5;
%
% Jacobi iteration loop.
[w,errVals(:,1),kStop(1)] = Jacobi(f,w,wExact,iterParam);
%
% Plot Jacobi errors.
figure(1)
clf
semilogy(errVals(1:kStop(1),1),'b-','LineWidth',1.5)
hold on
%
% Re-initialize w.
SEED = 1234;
rng(SEED);
w = rand(n,1)-0.5;
%
% Gauss-Seidel iteration loop.
[w,errVals(:,2),kStop(2)] = GaussSeidel(f,w,wExact,iterParam);
%
% Plot Gauss-Seidel errors.
figure(1)
semilogy(errVals(1:kStop(2),2),'r-','LineWidth',1.5)
hold on
%
% Re-initialize w.
SEED = 1234;
rng(SEED);
w = rand(n,1)-0.5;
%
% SOR iteration loop.
[w,errVals(:,3),kStop(3)] = SOR(f,w,wExact,iterParam);
%
% Plot SOR errors.
figure(1)
semilogy(errVals(1:kStop(3),3),'k-','LineWidth',1.5)
hold on
%
% Get contraction rates.
[rate] = contractionRates(errVals,kStop);
%
% Make a nice figure.
xlabel('$k$','Interpreter','latex');
title('Jacobi, Gauss-Seidel, and SOR Iteration Errors', ...
  'Interpreter','latex');
legend( ...
  '$\left\|{\bf w}-{\bf w}_{\rm J}^k\right\|$', ...
  '$\left\|{\bf w}-{\bf w}_{\rm GS}^k \right\|$', ...
  '$\left\|{\bf w}-{\bf w}_{\rm SOR}^k \right\|$', ...
  'Interpreter','latex');
text(0.7*kStop(1),13500*tol,strcat( ...
  '$n =\hspace{.1cm}$', ...
  num2str(n)),'FontSize',14, ...
  'Interpreter','latex')
text(0.7*kStop(1),4500*tol,strcat( ...
  '$\gamma_{\rm J} =\hspace{.1cm}$', ...
  num2str(rate(1),'%10.5e')),'FontSize',14, ...
  'Interpreter','latex')
text(0.7*kStop(1),1500*tol,strcat( ...
  '$\gamma_{\rm GS} =\hspace{.1cm}$', ...
  num2str(rate(2),'%10.5e')),'FontSize',14, ...
  'Interpreter','latex')
text(0.7*kStop(1),500*tol,strcat( ...
  '$\gamma_{\rm SOR} =\hspace{.1cm}$', ...
  num2str(rate(3),'%10.5e')),'FontSize',14, ...
  'Interpreter','latex')
printstr = strcat('Err_JGSSOR_',num2str(n),'.pdf');
exportgraphics(gca, printstr)
hold off