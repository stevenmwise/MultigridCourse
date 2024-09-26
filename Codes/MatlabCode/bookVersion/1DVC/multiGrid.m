function [u,rate] = multiGrid(f,MGParam,kMax,tol,uExact)

L = MGParam.L;
nL = MGParam.nL;

nl = nL;
for k = L-1:-1:0
  if mod(nl+1,2) ~= 0
    disp('Coarsening error: L is too large.')
    u = zeros(nL,1);
    rate = -1;
    return
  end
  nl = (nl+1)/2-1;
end
if nl < 1
  disp('Coarsening error: n0 >= 1 required.')
  u = zeros(nL,1);
  rate = -1;
  return
end

SEED = 1234;
rng(SEED);
u = rand(nL,1)-0.5;
rate = 1.0;

hL = 1/(nL+1);
x = linspace(hL,1.0-hL,nL)';

corr(1:kMax) = 0.0;
errv(1:kMax) = 0.0;
k = 0;
err = tol+1.0;

while ((k<kMax) & (err > tol))
  
  uo = u;
  u = MGOperator(f,u,L,MGParam);

  err = norm(uExact-u,2)
  corr(k+1) = norm(u-uo,2);
  errv(k+1) = err;
  k = k+1
  
end

clf

semilogy(errv(1:k),'bo','LineWidth',1.5)
hold on
semilogy(corr(1:k),'rs','LineWidth',1.5)
hold on

kv = [k-3:k];
le = log(errv(k-3:k));
p1 = polyfit(kv,le,1);
rate = exp(p1(1));
p1k = polyval(p1,[1:k]);
semilogy([1:k],exp(p1k),'-k')

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

end