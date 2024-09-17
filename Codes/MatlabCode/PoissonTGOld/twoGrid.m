function [u,rate] = twoGrid(f,n1,m1,omega,kMax,tol,uExact)

SEED = 1234;
rng(SEED);
u = rand(n1,1);
res = zeros(n1,1);
cgc = zeros(n1,1);
rate = 1.0;

h1 = 1/(n1+1);
x = linspace(h1,1.0-h1,n1)';

corr(1:kMax) = 0.0;
errv(1:kMax) = 0.0;
k = 0;
err = 10.0;

while ((k<kMax) & (err > tol))
    
    clf
    
    u = smooth(f,u,m1,omega);
    
    subplot(2,1,1)
    plot(x,uExact-u)
    xlabel('$x$','Interpreter','latex');
    ylabel('$u_1^{\rm E}-u_1^{(1)}$','Interpreter','latex');
    title('Error After Smoothing','Interpreter','latex');
    
    res = getResidual(f,u);
    cgc = coarseGridCorrection(res);
    u = u+cgc;
    
    subplot(2,1,2)
    plot(x,uExact-u)
    xlabel('$x$','Interpreter','latex');
    ylabel('$u_1^{\rm E}-u_1^{(2)}$','Interpreter','latex');
    title('Error After Coarse Grid Correction', ...
        'Interpreter','latex');
    
    printstr = strcat('Err_k_',num2str(k+1), ...
        '_n1_',num2str(n1), ...
        '_m1_',num2str(m1), ...
        '_omega_', num2str(omega),'.pdf');
    exportgraphics(gcf, printstr)

    err = norm(uExact-u,2)
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
text(1.5,10*tol,strcat('$\gamma_{\rm comp} =\hspace{.1cm}$', ...
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

function uSmooth = smooth(f,u,m1,omega)

n1 = length(f);
h1 = 1/(n1+1);
uSmooth = zeros(n1,1);
omegaPrime = 1.0-omega;

for j = 1:m1

    uSmooth(1) = (h1*f(1)+u(2))/2.0;
    uSmooth(2:n1-1)= (h1*f(2:n1-1)+u(1:n1-2)+u(3:n1))/2.0;
    uSmooth(n1) = (h1*f(n1)+u(n1-1))/2.0;

    uSmooth(1:n1) = omega*uSmooth(1:n1)+omegaPrime*u(1:n1);

    u = uSmooth;
    
end

end

function [res] = getResidual(f,u)

n1 = length(f);
h1 = 1.0/(n1+1);

res(1) = f(1)-(2.0*u(1)-u(2))/h1;

res(2:n1-1) = f(2:n1-1)-(-u(1:n1-2)+2*u(2:n1-1)-u(3:n1))/h1;

res(n1) = f(n1)-(-u(n1-1)+2.0*u(n1))/h1; 

end

function [cgc] = coarseGridCorrection(res)

n1 = length(res);
cgc = zeros(n1,1);
n0 = (n1+1)/2-1;
h0 = 1/(n0+1);
A0 = gallery('tridiag',n0,-1.0/h0,2.0/h0,-1.0/h0);

cgc = prolongation(A0\restriction(res));

end

function [vRest] = restriction(v)

n1 = length(v);
n0 = (n1+1)/2-1;
vRest = zeros(n0,1);

vRest(1:n0) = 0.5*(v(1:2:2*n0-1)+v(3:2:2*n0+1))+v(2:2:2*n0);

end

function [vProl] = prolongation(v)

n0 = length(v);
n1 = 2*(n0+1)-1;
vProl = zeros(n1,1);

vProl(1 ) = 0.5*v(1 );
vProl(n1) = 0.5*v(n0);

vProl(2:2:2*n0  ) = v(1:n0);
vProl(3:2:2*n0-1) = (v(1:n0-1)+v(2:n0))/2.0;

end