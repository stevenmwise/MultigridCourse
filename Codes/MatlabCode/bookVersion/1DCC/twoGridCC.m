function [u,rate] = twoGridCC(f,n1,m1,omega,kMax,tol,uExact)

SEED = 1234;
rng(SEED);
u = rand(n1+2,1)-0.5;
u(   1) = -u(   2);
u(n1+2) = -u(n1+1);
res = zeros(n1,1);
cgc = zeros(n1,1);
rate = 1.0;

h1 = 1/n1;

corr(1:kMax) = 0.0;
errv(1:kMax) = 0.0;
k = 0;
err = 10.0;

while ((k<kMax) & (err > tol))
  
  u(1:n1+2) = smoothCC(f(1:n1),u(1:n1+2),m1,omega);
  res(1:n1) = getResidualCC(f(1:n1),u(1:n1+2));
  cgc(1:n1) = coarseGridCorrectionCC(res(1:n1));
  u(2:n1+1) = u(2:n1+1)+cgc(1:n1);
  u(   1) = -u(   2);
  u(n1+2) = -u(n1+1);

  err = norm(uExact-u(2:n1+1),2)
  corr(k+1) = norm(cgc,2);
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
title('Two-grid Iteration Errors', ...
  'Interpreter','latex');
legend('$\|{\bf u}_1^{\rm E}-{\bf u}_1^k\|_1$', ...
  '$\|{\bf q}_1^{(1)} \|_1$','log-linear fit', ...
  'Interpreter','latex');
text(1.5,10*tol,strcat( ...
  '$\gamma_{\rm comp} =\hspace{.1cm}$', ...
  num2str(rate,'%10.5e')),'FontSize',14, ...
  'Interpreter','latex')
text(1.5,2*tol,strcat('$m_1 =\hspace{.1cm}$', ...
  num2str(m1)),'FontSize',14, ...
  'Interpreter','latex')
text(1.5,0.5*tol,strcat('$\omega =\hspace{.1cm}$', ...
  num2str(omega)),'FontSize',14, ...
  'Interpreter','latex')
printstr = strcat('Err_n1_',num2str(n1), ...
  '_m1_',num2str(m1), '_omega_', num2str(omega),'.pdf');
exportgraphics(gca, printstr)
hold off

end