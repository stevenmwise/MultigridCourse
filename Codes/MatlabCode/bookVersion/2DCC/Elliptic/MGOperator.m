function u = MGOperator(level,f,u,hf,MGParam)

% Get the number of interior cells from the size of u.
nfPlusGhost = size(u);

% For 2D, nf is a 2-element vector.
nf = nfPlusGhost-2;

% Pre-smoothing.
u = smoothQJacDamped(f,u,hf,MGParam.m1,MGParam.omega);

% Alternative smoothers (uncomment one of the following):
% u = smoothQGSDamped(f,u,hf,MGParam.m1,...
%   MGParam.omega,'Forward');

if level > 0
  hc = 2*hf;

  % nc is a 2-element vector [n1/2, n2/2].
  nc = nf/2;
    
  % Restrict the residual from the fine grid to the coarse grid.
  cGr = restriction(getResidual(u,f,hf,MGParam));
  
  % Initialize the coarse grid solution (with ghost cells).
  cGc = zeros(nc(1)+2,nc(2)+2);
  
  % Recursively solve on the coarse grid using p-cycles.
  for sig = 1:MGParam.pCycle
    cGc = MGOperator(level-1,cGr,cGc,hc,MGParam);
  end
  
  % Prolongate the coarse grid correction and add it to the 
  % fine grid solution.
  u(2:nf(1)+1,2:nf(2)+1) = u(2:nf(1)+1,2:nf(2)+1)+...
      prolongation(cGc(2:nc(1)+1,2:nc(2)+1));
  
  % Post-smoothing.
  u = smoothQJacDamped(f,u,hf,MGParam.m2,MGParam.omega);
  
  % Alternative smoothers (uncomment one of the following):
  % u = smoothQGSDamped(f,u,hf,...
  %   MGParam.m2,MGParam.omega,'Backward');
end

end
