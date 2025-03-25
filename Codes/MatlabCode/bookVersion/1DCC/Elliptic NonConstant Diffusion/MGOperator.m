function u = MGOperator(level,f,u,hf,MGParam,D)

n = length(u)-2;

% Pre-smoothing.
u = smoothQJacDamped(f,u,hf,D,MGParam.m1,MGParam.omega);

% Alternative pre-smoothing options:
% u = smoothRichardson(f,u,hf,...
%   D,MGParam.m1,MGParam.omega,MGParam);
% u = smoothQGSDamped(f,u,hf,...
%   D,MGParam.m1,MGParam.omega,'Forward');

if level > 0
  hc = 2*hf;
  nc = n/2;
  
  % Compute fine grid residual.
  fineRes = getResidual(u,f,hf,MGParam,D);
  
  % Restrict residual to the coarse grid.
  cGr = restriction(fineRes);
  
  % Restrict the diffusion coefficient to the coarse grid.
  cD = restrictDiff(D);
  
  % Initialize coarse grid solution (with ghost cells).
  cGc = zeros(nc+2,1);
  
  % Solve on the coarse level using pCycle cycles.
  for sig = 1:MGParam.pCycle
    cGc = MGOperator(level-1,cGr,cGc,hc,MGParam,cD);
  end
  
  % Prolongate the coarse grid correction and update the...
  % fine grid solution.
  u(2:n+1) = u(2:n+1)+prolongation(cGc(2:nc+1));
  
  % Post-smoothing.
  u = smoothQJacDamped(f,u,hf,D,MGParam.m2,MGParam.omega);
  
  % Alternative post-smoothing options:
  % u = smoothRichardson(f,u,hf,...
  %   D,MGParam.m2,MGParam.omega,MGParam);
  % u = smoothQGSDamped(f,u,hf,D,MGParam.m2,...
  %   MGParam.omega,'Backward');
end

end
