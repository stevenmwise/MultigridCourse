tol = 1.0e-09;
nL = 4096;
L = 12;
xLower = 0.0;
xUpper = 3.2;
hL = (xUpper-xLower)./nL;

pCycle = 1;
m1 = 2;
m2 = 2;
omega = 1.0;
kMax = 100;
%
% Check mesh coarsening.
if nL < 0
  disp('nL must be a positive integer.')
  return
end
nl = nL;
for k = L-1:-1:0
  if mod(nl,2) ~= 0
    disp('Coarsening error: L is too large.')
    return
  end
  nl = nl/2;
end

MGParam.nL = nL;
MGParam.xLower = xLower;
MGParam.xUpper = xUpper;
MGParam.L = L;
MGParam.pCycle = pCycle;
MGParam.m1 = m1;
MGParam.m2 = m2;
MGParam.omega = omega;
MGParam.kMax = kMax;
MGParam.tol = tol;

SEED = 1234;
rng(SEED);
u = rand(nL+2,1)-0.5;
uExact = zeros(nL+2,1);
f = zeros(nL,1);
    
px = xUpper-xLower;

for i = 1:nL
  x = (i-0.5)*hL+xLower;
  uExact(i+1) = exp(cos(2.0*pi*x/px));
end

%
% Manufacture the forcing vector.
uExact = applyBCs(uExact);
f = FDOperator(uExact,hL);
%
% Get the approximate solution.
[u,errVals,kStop] = multiGridSolver(u,f,hL,MGParam,uExact);
uPlot(1:nL) = u(2:nL+1);
%
% Plot solution details.
x =linspace(xLower+hL/2,xUpper-hL/2,nL);

figure(1)
clf

plot(x,uPlot,'k-');
title('Numerical Solution u');
xlabel('x');
ylabel('u(x)');
axis equal tight;

figure(2)
clf

semilogy(errVals(1:kStop,3),'bo','LineWidth',1.5)
hold on
semilogy(errVals(1:kStop,2),'rs','LineWidth',1.5)
hold on
semilogy(errVals(1:kStop,1),'kd','LineWidth',1.5)
hold on

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
  '$\left\|{\bf r}_L^k\right\|_L$', ...
  'log-linear fit','Interpreter','latex');
text(1.5,128*MGParam.tol,strcat( ...
  '$\gamma_{\rm comp} =\hspace{.1cm}$', ...
  num2str(rate,'%10.5e')),'FontSize',14, ...
  'Interpreter','latex')
text(1.5,32*MGParam.tol,strcat('$m_1 =\hspace{.1cm}$', ...
  num2str(MGParam.m1)),'FontSize',14, ...
  'Interpreter','latex')
text(1.5,8*MGParam.tol,strcat('$m_2 =\hspace{.1cm}$', ...
  num2str(MGParam.m2)),'FontSize',14, ...
  'Interpreter','latex')
text(1.5,2*MGParam.tol,strcat('$p =\hspace{.1cm}$', ...
  num2str(MGParam.pCycle)),'FontSize',14, ...
  'Interpreter','latex')
printstr = strcat('Err_nL_',num2str(MGParam.nL), ...
  '_m1_',num2str(MGParam.m1),'.pdf');
exportgraphics(gca, printstr)
hold off