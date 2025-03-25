function [u,errVals,kStop] = twoGridSolver(u,f,n1,m1,...
    omega,kMax,tol,uExact)

errVals = zeros(kMax,3);
res = zeros(n1,1);
cgc = zeros(n1,1);

h1 = 1/(n1+1);
x = linspace(h1,1.0-h1,n1)';

k = 0;
err = tol+1.0;

while ((k < kMax) && (err > tol))
    
  % Pre-smoothing (using damped Jacobi by default).
  u = smoothJacobiDamped(f,u,m1,omega);

  % Alternative pre-smoothing options:
  % u = smoothJacobi(f,u,m1);
  % u = smoothGaussSeidel(f,u,m1,'Forward');

  % Compute and apply coarse grid correction.
  res = getResidual(f,u);
  cgc = coarseGridCorrection(res);
  u = u+cgc;

  % Compute errors.
  err = norm(uExact-u,2);
  errVals(k+1,2) = norm(cgc,2);
  errVals(k+1,3) = err;

  k = k+1;
  
  % Print iteration information.
  fprintf('Two-grid iteration %d:\n',k);
  fprintf('   residual  = %.6e\n',err);
  fprintf('   correction= %.6e\n',errVals(k,2));
  fprintf('   true error= %.6e\n',errVals(k,3));
end

kStop = k;

end