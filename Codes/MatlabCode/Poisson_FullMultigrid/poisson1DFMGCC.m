clear
clc

tol    = 1.0e-09;
nL     = 2^5;      
L      = 3;        
xLower = 0.0;
xUpper = 1.0;
hL     = (xUpper - xLower)./nL;

pCycle = 2;
m1     = 3;
m2     = 3;
omega  = 2/3;
kMax   = 100;

if nL < 1
  disp('nL must be a positive integer.')
  return
end

nl = nL;
for kk = L-1 : -1 : 0
  if mod(nl,2) ~= 0
    disp('Coarsening error: L is too large or nL not divisible by 2^L.')
    return
  end
  nl = nl/2;
end


MGParam.nL     = nL;
MGParam.xLower = xLower;
MGParam.xUpper = xUpper;
MGParam.L      = L;
MGParam.pCycle = pCycle;
MGParam.m1     = m1;
MGParam.m2     = m2;
MGParam.omega  = omega;
MGParam.kMax   = kMax;
MGParam.tol    = tol;


MGParam.levelHistory = [];
MGParam.colorHistory = [];

MGParam.levelHistory(end+1) = 0;    

uExact = zeros(nL+2,1);
f      = zeros(nL,1);
px     = xUpper - xLower;

for i = 1:nL
  x = (i - 0.5)*hL + xLower;
  uExact(i+1) = sin(2.0*pi*x/px)^3;   % your chosen exact function
end
uExact = applyBCs(uExact);
f      = FDOperator(uExact,hL);


[u, MGParam, errVals, kStop] = fullMultiGridSolver(MGParam, uExact);


uPlot = u(2:nL+1);
x     = linspace(xLower + hL/2, xUpper - hL/2, nL);

figure(1)
clf
plot(x, uPlot, 'k-');
title('Numerical Solution u');
xlabel('x');
ylabel('u(x)');
axis tight;

figure(2)
clf
semilogy(errVals(1:kStop,3), 'bo', 'LineWidth',1.5)
hold on
semilogy(errVals(1:kStop,2), 'rs', 'LineWidth',1.5)
semilogy(errVals(1:kStop,1), 'kd', 'LineWidth',1.5)

rate = 1.0;
if kStop >= 4
  kv = [kStop-3 : kStop];
  le = log(errVals(kv,3));
  p1 = polyfit(kv, le, 1);
  rate = exp(p1(1));
  p1k = polyval(p1, [1:kStop]);
  semilogy([1:kStop], exp(p1k), '-k')
end

xlabel('$k$','Interpreter','latex');
title('Multigrid Iteration Errors','Interpreter','latex');
legend('$\|u^E - u^k\|_L$', ...
       '$\|u^k - u^{k-1}\|_L$', ...
       '$\|r^k\|_L$', ...
       'log-linear fit', ...
       'Interpreter','latex');

text(1.5,128*tol, ...
  ['$\gamma_{\mathrm{comp}} =$ ', num2str(rate,'%10.5e')], ...
  'FontSize',14,'Interpreter','latex')
text(1.5,32*tol, ['$m_1 =$ ', num2str(MGParam.m1)], ...
  'FontSize',14,'Interpreter','latex')
text(1.5,8*tol, ['$m_2 =$ ', num2str(MGParam.m2)], ...
  'FontSize',14,'Interpreter','latex')
text(1.5,2*tol, ['$p =$ ', num2str(MGParam.pCycle)], ...
  'FontSize',14,'Interpreter','latex')
text(1.5,0.5*tol, ['$\omega =$ ', num2str(MGParam.omega,'%10.3f')], ...
  'FontSize',14,'Interpreter','latex')

hold off

if length(MGParam.levelHistory) <= 2000
  figure(3)
  clf
  LH = MGParam.levelHistory;    % short alias
  CH = MGParam.colorHistory;    % same length-1 as LH, typically

  if length(LH) < 2
    text(0.5,0.5,'No or insufficient levels recorded','HorizontalAlignment','center')
    title('Path of Levels Visited During Full Multigrid')
    xlabel('Step #')
    ylabel('Grid Level')
    axis([0 1 0 1])
    return
  end

  hold on
  for i = 1 : (length(LH)-1)
    xvals = [i, i+1];
    yvals = [LH(i), LH(i+1)];
  
    if CH(i) == 1
      % Outer FMG step => red
      plot(xvals, yvals, '-or','LineWidth',1.5,'MarkerSize',5);
    else
      % V-cycle step => blue
      plot(xvals, yvals, '-ob','LineWidth',1.5,'MarkerSize',5);
    end
  end

  title('Path of Levels Visited During Full Multigrid')
  xlabel('Step #')
  ylabel('Grid Level')
  xlim([0, length(LH)+1]);
  ylim([-0.5, max(LH)+0.5])
  grid on
  hold off
else
    fprintf("Level traversal history is not plotted since the number of iterations is large.\n")
end
