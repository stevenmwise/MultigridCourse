function u = MGOperator(level,f,u,hf,MGParam)

nPlusGhostLayers = length(u);
nf = nPlusGhostLayers-2;

% Pre-smoothing using damped Jacobi.
u = smoothQJacDamped(f,u,hf,MGParam.m1,MGParam.omega);

% Alternative pre-smoothing options:
% u = smoothRichardson(f,u,hf,MGParam.m1,MGParam.omega);
% u = smoothQGSDamped(f,u,hf,...
%   MGParam.m1,MGParam.omega,'Forward');

if level > 0
  hc = 2*hf;
  nc = nf/2;
  
  % Compute fine grid residual.
  cGr = restriction(getResidual(u,f,hf));
  
  % Initialize coarse grid solution (with ghost cells).
  cGc = zeros(nc+2, 1);
  
  for sig = 1:MGParam.pCycle
    cGc = MGOperator(level-1,cGr,cGc,hc,MGParam);
  end
  
  % Prolongate the coarse grid correction and update the...
  % fine grid solution.
  u(2:nf+1) = u(2:nf+1)+prolongation(cGc(2:nc+1));
  
  % Post-smoothing using damped Jacobi.
   u = smoothQJacDamped(f,u,hf,MGParam.m2,MGParam.omega);

  % Alternative post-smoothing options:
  % u = smoothRichardson(f,u,hf,MGParam.m2,MGParam.omega);
  % u = smoothQGSDamped(f,u,hf,...
  %  MGParam.m2,MGParam.omega,'Backward');
end

end
