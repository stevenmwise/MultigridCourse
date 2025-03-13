function [u,errVals,kStop] = multiGridSolver(u0,f,hf, ...
    MGParam,uExact,D_EW,D_NS)

nf = size(f);
k = 0;
errorVal = MGParam.tol+1.0;
errVals = zeros(MGParam.kMax,3);

uo = u0;

while ((k < MGParam.kMax) && (errorVal > MGParam.tol))

  % Apply one multigrid cycle.
  u = MGOperator(MGParam.L,f,uo,hf,MGParam,D_EW,D_NS);

  % Evaluate various norms.
  resNow = normScaledL2(getResidual(u,f,hf,MGParam,D_EW,D_NS));
  corrNow = normScaledL2(uo(2:nf(1)+1,2:nf(2)+1)- ...
      u(2:nf(1)+1,2:nf(2)+1));
  exactErr = normScaledL2(uExact(2:nf(1)+1,2:nf(2)+1)- ...
      u(2:nf(1)+1,2:nf(2)+1));

  k = k+1;
  errVals(k,1) = resNow;
  errVals(k,2) = corrNow;
  errVals(k,3) = exactErr;

  fprintf('Multigrid iteration %d:\n', k);
  fprintf('   residual  = %.6e\n', resNow);
  fprintf('   correction= %.6e\n', corrNow);
  fprintf('   true error= %.6e\n', exactErr);

  % Stopping criterion (artificial, for demonstration).
  errorVal = exactErr;  
  uo = u;
end

kStop = k;
end
