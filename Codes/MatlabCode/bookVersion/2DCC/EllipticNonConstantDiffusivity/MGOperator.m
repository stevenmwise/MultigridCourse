function u = MGOperator(level,f,u,hf,MGParam,D_EW,D_NS)

% Pre-smoothing
u = smoothJacobiDamped(f,u,hf,D_EW,D_NS, ...
    MGParam.m1,MGParam.omega);

% Alternative pre-smoothing options:
% u = smoothJacobi(f,u,hf,D_EW,D_NS,MGParam.m1);
% u = smoothGaussSeidel(f,u,hf,D_EW,D_NS, ...
%     MGParam.m1,'Forward');

if level > 0
  hc = 2*hf;
  nfPlusGhost = size(u);
  nf = nfPlusGhost-2;

  % Coarse grid interior dimensions is [nc(1),nc(2)].
  nc = nf/2; 
    
  % Compute the residual on the fine grid.
  fineRes = getResidual(u,f,hf,MGParam,D_EW,D_NS);
    
  % Restrict the residual to the coarse grid.
  cGr = restriction(fineRes);
    
  % Restrict the diffusion coefficients to the coarse grid.
  [cD_EW,cD_NS] = restrictDiff(D_EW,D_NS);
    
  % Solve on the coarse grid using p-cycles.
  cGc = zeros(nc(1)+2,nc(2)+2);
  for s = 1:MGParam.pCycle
    cGc = MGOperator(level-1,cGr,cGc,hc,MGParam,cD_EW,cD_NS);
  end
    
  % Prolongate the coarse grid correction and update the ...
  % fine grid solution.
  u(2:nf(1)+1,2:nf(2)+1) = u(2:nf(1)+1,2:nf(2)+1)+ ...
      prolongation(cGc(2:nc(1)+1,2:nc(2)+1));
    
  % Post-smoothing.
  u = smoothJacobiDamped(f,u,hf,D_EW,D_NS, ...
      MGParam.m2,MGParam.omega);
    
  % Alternative post-smoothing options:
  % u = smoothJacobi(f,u,hf,D_EW,D_NS,MGParam.m2);
  % u = smoothGaussSeidel(f,u,hf,D_EW,D_NS, ...
  %     MGParam.m2,'Backward');
end
end
