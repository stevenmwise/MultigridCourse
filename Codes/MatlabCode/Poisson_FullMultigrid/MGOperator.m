function [u, MGParam] = MGOperator(u, f, hf, level, MGParam)
%
% MGOPERATOR  A single V-cycle step at the given 'level'. 
%
  nfPlusGhostLayers = length(u);
  nf = nfPlusGhostLayers - 2;

  % Pre-smoothing
  u = smooth(u, f, hf, MGParam.m1, MGParam.omega);

  if level > 0

    hc = 2*hf;
    nc = nf/2;

    % Restrict the residual
    resFine = getResidual(u, f, hf);
    cGr = restriction(resFine);
    cGc = zeros(nc+2,1);

    % Recursively solve on coarse level
    for sig = 1:MGParam.pCycle

      % Going down to (level-1) => color = 0 (blue)
      MGParam.levelHistory(end+1) = level - 1;
      MGParam.colorHistory(end+1) = 0;

      [cGc, MGParam] = MGOperator(cGc, cGr, hc, level-1, MGParam);

      % Coming back up to 'level' => color = 0 (blue)
      MGParam.levelHistory(end+1) = level;
      MGParam.colorHistory(end+1) = 0;

    end

    % Prolong coarse correction
    u(2:nf+1) = u(2:nf+1) + prolongation(cGc(2:nc+1));

  end

  % Post-smoothing
  u = smooth(u, f, hf, MGParam.m2, MGParam.omega);

end
