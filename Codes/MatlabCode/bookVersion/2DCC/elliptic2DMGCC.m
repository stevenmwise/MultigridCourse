tol = 1.0e-11;
nL = [2048,2048];
L = 11;
xLower = [0.0,0.0];
xUpper = [3.2,3.2];
hhL = (xUpper-xLower)./nL;

pCycle = 1;
m1 = 2;
m2 = 2;
omega = 1.0;
kMax = 100;

%
% Check mesh spacing.
if abs(hhL(1)-hhL(2))>1.0e-10
  print('Grid size error: hL.')
  return
else
  hL = hhL(1);
end
%
% Check mesh coarsening.
%
if any(nL < 0)
  disp('nL entries must be positive integers.')
  return
end
nl = nL;
for k = L-1:-1:0
  if any(mod(nl,2) ~= 0)
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
u = rand(nL(1)+2,nL(2)+2)-0.5;
uExact = zeros(nL(1)+2,nL(2)+2);
f = zeros(nL(1),nL(2));
    
px = xUpper-xLower;

for i = 1:nL(1)
  for j = 1:nL(2)
    x = (i-0.5)*hL+xLower(1);
    y = (j-0.5)*hL+xLower(2);
    uExact(i+1,j+1) = exp(cos(2.0*pi*x/px(1))*cos(2.0*pi*y/px(2)));
  end
end
%
% Manufacture the forcing vector.
uExact = applyBCs(uExact);
f = FDOperator(uExact,hL);
%
% Get the approximate solution.
[u,errVals,kStop] = multiGridSolver(u,f,hL,MGParam,uExact);
uPlot(1:nL(1),1:nL(1)) = u(2:nL(1)+1,2:nL(2)+1);
%
% Plot solution details.
[x, y] = meshgrid(linspace(xLower(1)+hL/2,xUpper(1)-hL/2,nL(1)), ...
                  linspace(xLower(2)+hL/2,xUpper(2)-hL/2,nL(2)));

figure(1)
clf

surf(x,y,uPlot,'EdgeColor','none');
title('Numerical Solution u');
xlabel('x');
ylabel('y');
zlabel('u(x,y)');
colorbar;
view(2); % View from top
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
printstr = strcat('Err_nL_',num2str(MGParam.nL(1)), ...
  '_m1_',num2str(MGParam.m1),'.pdf');
exportgraphics(gca, printstr)
hold off