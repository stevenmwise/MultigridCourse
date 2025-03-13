clear

n = 3;
kMax = 50000;
tol = 1e-12;
maxLevel = 6;
iterParam.kMax = kMax;
iterParam.tol = tol;
%
% Initialize w.
SEED = 1234;
rng(SEED);
w = rand(n,1)-0.5;

errVals = zeros(kMax,1);
kStop = zeros(maxLevel,1);

for level = 1:maxLevel
%
% Increase grid size.
  n = 2*(n+1)-1;
  h = 1/(n+1);
  h2 = h*h;
%
% Store variables in a parameter structure.
  iterParam.n = n;
  iterParam.h = h;
  iterParam.h2 = h2;
%
% Construct the right-hand-side vector, f.
  x = linspace(h,1.0-h,n)';
  A = gallery('tridiag',n,-1.0/h2,2.0/h2,-1.0/h2);
  wExact = exp(sin(3.0*pi*x))-1.0;
  f = A*wExact;
%
% Move approximation from lower level.
  w = prolongation(w);
  figure(2*level-1)
  clf
  plot(x,w)
  hold on
%
% Jacobi iteration loop.
  [w,errVals(:,level),kStop(level)] = GaussSeidel(f,w,wExact,iterParam);  
%
  plot(x,w)
%
% Plot Jacobi errors.
  figure(2*level)
  clf
  semilogy(errVals(1:kStop(level),level),'b-','LineWidth',1.5)
  hold on
end 

% Make a nice figure.
%xlabel('$k$','Interpreter','latex');
%title('Jacobi-Cascade Errors','Interpreter','latex');
%printstr = strcat('Err_Cascade_',num2str(n),'.pdf');
%exportgraphics(gca, printstr)

hold off
