function [u,errVals,kStop] = multiGridSolver(uo,f,hf,MGParam,uExact)

  nf = length(f);
  k = 0;
  error = MGParam.tol+1.0;
  errVals = zeros(MGParam.kMax,3);
  
  while ((k<MGParam.kMax) & (error > MGParam.tol))

    u = MGOperator(uo,f,hf,MGParam.L,MGParam);

%
% Various error measurements. If the exact solution is not
% available, use the residual, errVals(k+1,1), or the correction,
% errVals(k+1,2), as an error estimate:
    errVals(k+1,1) = normScaledL2(getResidual(u,f,hf));
    errVals(k+1,2) = normScaledL2(    uo(2:nf+1)-u(2:nf+1));
    errVals(k+1,3) = normScaledL2(uExact(2:nf+1)-u(2:nf+1));
%
% Print errors to the screen.
    fprintf('Multigrid iteration %d: Residual   l2 = %.6e\n', ...
      k+1,errVals(k+1,1));
    fprintf('                       Correction l2 = %.6e\n', ...
          errVals(k+1,2));
    fprintf('                       Error      l2 = %.6e\n', ...
          errVals(k+1,3));
%
% Using the true error for stopping. This is not practical.
    error = errVals(k+1,3);
    k = k+1;
    uo = u;
  end

  kStop = k;

end