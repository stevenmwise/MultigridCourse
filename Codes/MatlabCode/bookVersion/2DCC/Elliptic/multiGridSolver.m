function [u,errVals,kStop] = multiGridSolver(uo,f,hf,...
    MGParam,uExact)

nf = size(f);
k = 0;
error = MGParam.tol+1.0;
errVals = zeros(MGParam.kMax,3);

while ((k < MGParam.kMax) && (error > MGParam.tol))
  
  u = MGOperator(MGParam.L,f,uo,hf,MGParam);
  
  % Various error measurements. If the exact solution is not
  % available, use the residual (errVals(k+1,1)) or...
  % the correction (errVals(k+1,2)) as an error estimate.
  errVals(k+1,1) = normScaledL2(getResidual(u,f,hf,MGParam));
  errVals(k+1,2) = normScaledL2(uo(2:nf(1)+1,2:nf(2)+1)-...
    u(2:nf(1)+1,2:nf(2)+1));
  errVals(k+1,3) = normScaledL2(uExact(2:nf(1)+1,2:nf(2)+1)-...
    u(2:nf(1)+1,2:nf(2)+1));
  
  fprintf('Multigrid iteration %d:\n',k+1);
  fprintf('   residual  = %.6e\n',errVals(k+1,1));
  fprintf('   correction= %.6e\n',errVals(k+1,2));
  fprintf('   true error= %.6e\n',errVals(k+1,3));
  
  error = errVals(k+1,3);
  k = k+1;
  uo = u;
end

kStop = k;

end
