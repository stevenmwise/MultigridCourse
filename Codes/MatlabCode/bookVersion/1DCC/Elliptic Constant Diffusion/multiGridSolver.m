function [u,errVals,kStop] = multiGridSolver(uo,f,hf, ...
    MGParam,uExact)

nf = length(f);
k = 0;
error = MGParam.tol + 1.0;
errVals = zeros(MGParam.kMax,3);

while ((k < MGParam.kMax) && (error > MGParam.tol))
  
  u = MGOperator(MGParam.L,f,uo,hf,MGParam);
  
  % Compute various error measurements.
  resNow = normScaledL2(getResidual(u,f,hf));
  corrNow = normScaledL2(uo(2:nf+1) - u(2:nf+1));
  exactErr = normScaledL2(uExact(2:nf+1) - u(2:nf+1));
  
  k = k+1;
  errVals(k,1) = resNow;
  errVals(k,2) = corrNow;
  errVals(k,3) = exactErr;
  
  fprintf('Multigrid iteration %d:\n', k);
  fprintf('   residual  = %.6e\n', resNow);
  fprintf('   correction= %.6e\n', corrNow);
  fprintf('   true error= %.6e\n', exactErr);
  
  error = exactErr;
  uo = u;
end

kStop = k;
end
