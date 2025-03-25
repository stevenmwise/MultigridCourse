function [u,errVals,kStop] = multiGridSolver(u0,f,hf,...
    MGParam,uExact)

nf = length(f);
k = 0;
errorVal = MGParam.tol+1.0;
errVals = zeros(MGParam.kMax,3);
uo = u0;

while ((k < MGParam.kMax) && (errorVal > MGParam.tol))
  
  u = MGOperator(MGParam.L,f,uo,hf,MGParam);
  
  % Compute various error measurements.
  errVals(k+1,1) = normScaledL2(getResidual(u,f,hf));
  errVals(k+1,2) = normScaledL2(uo(2:nf+1)-u(2:nf+1));
  errVals(k+1,3) = normScaledL2(uExact(2:nf+1)-u(2:nf+1));
  
  fprintf('Multigrid iteration %d:\n', k+1);
  fprintf('   residual  = %.6e\n', errVals(k+1,1));
  fprintf('   correction= %.6e\n', errVals(k+1,2));
  fprintf('   true error= %.6e\n', errVals(k+1,3));
  
  errorVal = errVals(k+1,3);
  k = k+1;
  uo = u;
end

kStop = k;

end
