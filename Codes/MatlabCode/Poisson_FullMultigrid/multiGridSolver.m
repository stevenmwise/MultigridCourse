function [u, MGParam, errVals, kStop] = multiGridSolver(uo, f, hf, MGParam, uExact)
%
% MULTIGRIDSOLVER  Repeated V-cycles on the current grid (spacing hf),
% with recursion down to level=0 if MGParam.L>0.
%

  nf      = length(f);
  k       = 0;
  error   = MGParam.tol + 1.0;
  errVals = zeros(MGParam.kMax,3);

  fprintf('\n=== Starting V-cycle on Level %d ===\n', MGParam.L);

  while (k < MGParam.kMax) && (error > MGParam.tol)

    % Perform one V-cycle
    [u, MGParam] = MGOperator(uo, f, hf, MGParam.L, MGParam);

    % Evaluate various error measurements
    r = getResidual(u, f, hf);
    errVals(k+1,1) = normScaledL2(r);
    errVals(k+1,2) = normScaledL2(uo(2:nf+1) - u(2:nf+1));
    errVals(k+1,3) = normScaledL2(uExact(2:nf+1) - u(2:nf+1));

    % Detailed, informative output
    fprintf('Level %d, Iteration %d:\n', MGParam.L, k+1);
    fprintf('    Residual   l2 = %.6e\n', errVals(k+1,1));
    fprintf('    Correction l2 = %.6e\n', errVals(k+1,2));
    fprintf('    Error      l2 = %.6e\n', errVals(k+1,3));

    % Update stopping criteria
    error = errVals(k+1,3);
    k     = k + 1;
    uo    = u;
  end

  fprintf('=== Finished V-cycle on Level %d after %d iterations ===\n\n', MGParam.L, k);

  kStop = k;
end
