function [u,errVals,kStop] = multiGridSolver(u0,f, ...
    hf,MGParam,uExact,D)

n = length(f);
k = 0;
errorVal = MGParam.tol+1.0;
errVals = zeros(MGParam.kMax,3);
uo = u0;

while ((k < MGParam.kMax) && (errorVal > MGParam.tol))
  
  u = MGOperator(MGParam.L,f,uo,hf,MGParam,D);
  
  % Compute various error measurements.
  resNorm = normScaledL2(getResidual(u,f,hf,MGParam,D));
  corrNorm = normScaledL2(uo(2:n+1)-u(2:n+1));
  errNorm = normScaledL2(uExact(2:n+1)-u(2:n+1));
  
  k = k+1;
  errVals(k,1) = resNorm;
  errVals(k,2) = corrNorm;
  errVals(k,3) = errNorm;
  
  fprintf('Multigrid iteration %d:\n', k);
  fprintf('   residual  = %.6e\n', resNorm);
  fprintf('   correction= %.6e\n', corrNorm);
  fprintf('   true error= %.6e\n', errNorm);
  
  errorVal = errNorm;
  uo = u;
end

kStop = k;
end
