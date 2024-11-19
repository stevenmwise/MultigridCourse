function [u,errVals,kStop] = twoGridSolver(u,f,n1,m1,omega,kMax, ...
  tol,uExact)

errVals = zeros(kMax,3);

res = zeros(n1,1);
cgc = zeros(n1,1);
rate = 1.0;

h1 = 1/(n1+1);
x = linspace(h1,1.0-h1,n1)';

k = 0;
err = tol+1.0;

while ((k<kMax) & (err > tol))
  
  u = smooth(f,u,m1,omega);
  res = getResidual(f,u);
  cgc = coarseGridCorrection(res);
  u = u+cgc;

  err = norm(uExact-u,2)
  errVals(k+1,2) = norm(cgc,2);
  errVals(k+1,3) = err;
  k = k+1
  
end

kStop = k;

end