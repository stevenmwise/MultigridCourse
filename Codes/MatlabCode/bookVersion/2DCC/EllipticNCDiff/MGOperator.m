function u = MGOperator(level,f,DeW,DnS,u,hf,MGParam)

rC = MGParam.rC;

% Pre-smoothing
u = smoothQJacDamped(f,DeW,DnS,rC,u,hf,MGParam.m1, ...
  MGParam.omega);

% Alternative pre-smoothing option:
% u = smoothQGSDamped(f,DeW,DnS,rC,u,hf, ...
%   MGParam.m1,MGParam.omega,'fwd');

if level > 0
  hc = 2*hf;
  nfPlusGhost = size(u);
  nf = nfPlusGhost-2;

% Coarse grid interior dimensions are [nc(1),nc(2)].
  nc = nf/2;

% Compute the residual and restric to the coarse grid.
  cGr = restriction(getResidual(u,f,DeW,DnS,hf,MGParam));
    
% Restrict the diffusion coefficients to the coarse grid.
  [cDeW,cDnS] = restrictDiff(DeW,DnS);
    
% Approximate on the coarse grid using recursive MG.
  cGc = zeros(nc(1)+2,nc(2)+2);
  for s = 1:MGParam.pCycle
    cGc = MGOperator(level-1,cGr,cDeW,cDnS,cGc,hc,MGParam);
  end

% Prolongate the coarse grid correction and update the
% fine grid approximation.
  u(2:nf(1)+1,2:nf(2)+1) = u(2:nf(1)+1,2:nf(2)+1) ...
    + prolongation(cGc(2:nc(1)+1,2:nc(2)+1));
    
% Post-smoothing.
  u = smoothQJacDamped(f,DeW,DnS,rC,u,hf,MGParam.m2, ...
    MGParam.omega);

% Alternative post-smoothing option:
%   u = smoothQGSDamped(f,DeW,DnS,rC,u,hf, ...
%     MGParam.m2,MGParam.omega,'bwd');
end

end