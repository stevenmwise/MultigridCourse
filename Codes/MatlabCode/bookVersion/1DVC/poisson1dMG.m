clear

%
% Set multigrid parameters.
nL = 255;
hL = 1/(nL+1);
L = 7;
p = 1;
m1 = 3;
m2 = 3;
omega = 2/3;
tol = 1e-09;
kMax = 20;

%
% Check mesh coarsening.
nl = nL;
for k = L-1:-1:0
  if mod(nl+1,2) ~= 0
    disp('Coarsening error: L is too large.')
    return
  end
  nl = (nl+1)/2-1;
end

%
% Give a random initial guess.
SEED = 1234;
rng(SEED);
u = rand(nL,1)-0.5;
rate = 1.0;

%
% Construct the exact solution and necessary forcing vector.
x = linspace(hL,1.0-hL,nL)';
AL = gallery('tridiag',nL,-1.0/hL,2.0/hL,-1.0/hL);
uExact = exp(sin(3.0*pi*x))-1.0;
f = AL*uExact;

MGParam.nL = nL;
MGParam.L = L;
MGParam.p = p;
MGParam.m1 = m1;
MGParam.m2 = m2;
MGParam.omega = omega;
MGParam.tol = tol;
MGParam.kMax = kMax;

%
% Call the solver.
[u,errVals,kStop] = multiGridSolver(u,f,MGParam,uExact);

figure(1)
clf

plot(x,u,'k-','LineWidth',1.5)
xlabel('$x$','Interpreter','latex');
ylabel('$u_L(x)$','Interpreter','latex');
title('Finite Element Approximation of the Poisson Equation');
axis([0.0,1.0,-1.0,2.0])
printstr = strcat('Solution_nL_',num2str(nL),'.pdf');
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
legend('$\left\|{\bf u}_L^{\rm E}-{\bf u}_L^k\right\|_L$', ...
  '$\left\|{\bf u}_L^k-{\bf u}_L^{k-1} \right\|_L$', ...
  'log-linear fit','Interpreter','latex');
text(1.5,128*tol,strcat( ...
  '$\gamma_{\rm comp} =\hspace{.1cm}$', ...
  num2str(rate,'%10.5e')),'FontSize',14, ...
  'Interpreter','latex')
text(1.5,32*tol,strcat('$m_1 =\hspace{.1cm}$', ...
  num2str(MGParam.m1)),'FontSize',14, ...
  'Interpreter','latex')
text(1.5,8*tol,strcat('$m_2 =\hspace{.1cm}$', ...
  num2str(MGParam.m2)),'FontSize',14, ...
  'Interpreter','latex')
text(1.5,2*tol,strcat('$p =\hspace{.1cm}$', ...
  num2str(MGParam.p)),'FontSize',14, ...
  'Interpreter','latex')
text(1.5,0.5*tol,strcat('$\omega =\hspace{.1cm}$', ...
  num2str(MGParam.omega)),'FontSize',14, ...
  'Interpreter','latex')
printstr = strcat('Err_nL_',num2str(nL), ...
  '_m1_',num2str(MGParam.m1), '_omega_', ...
  num2str(MGParam.omega),'.pdf');
exportgraphics(gca, printstr)
hold off