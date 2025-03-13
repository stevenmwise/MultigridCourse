function [u,errVals,kStop] = multiGridSolver(u,f, ...
    MGParam,uExact)

L = MGParam.L;
nL = MGParam.nL;
kMax = MGParam.kMax;
tol = MGParam.tol;
errVals = zeros(kMax,3);

hL = 1/(nL+1);
x = linspace(hL,1.0-hL,nL)';

k = 0;
err = tol+1.0;

while ((k < kMax) && (err > tol))
  uo = u;
  u = MGOperator(L,f,u,MGParam);

  err = norm(uExact-u,2);
  errVals(k+1,2) = norm(u-uo,2);
  errVals(k+1,3) = err;

  k = k+1;
  
  % Print iteration information.
  fprintf('Multigrid iteration %d:\n',k);
  fprintf('   residual  = %.6e\n',err);
  fprintf('   correction= %.6e\n',errVals(k,2));
  fprintf('   true error= %.6e\n',errVals(k,3));
end

kStop = k;

end