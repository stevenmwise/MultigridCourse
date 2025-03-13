function u = MGOperator(level,f,u,hf,MGParam)

[nFull,~] = size(u);

n = nFull-2;

% Pre-smoothing.
u = smoothJacobiDamped(f,u,hf,MGParam.m1,MGParam.omega);

% Alternative pre-smoothing options:
% u = smoothJacobi(f,u,hf,MGParam.m1);
% u = smoothGaussSeidel(f,u,hf,MGParam.m1, 'Forward');

if level > 0
  hc = 2*hf;
  
  % Assumes n is divisible by 2.
  nc = n/2;  
  
  % Restrict the residual from the fine grid to the coarse grid.
  cGr = restriction(getResidual(u,f,hf,MGParam));
  
  % Initialize coarse grid solution (with ghost cells).
  cGc = zeros(nc+2,nc+2);
  
  % Solve on the coarse grid using p-cycle.
  for sig = 1:MGParam.pCycle
    cGc = MGOperator(level-1,cGr,cGc,hc,MGParam);
  end
  
  % Prolongate the coarse grid correction and update the ...
  % fine grid solution.
  u(2:n+1,2:n+1) = u(2:n+1,2:n+1)+ ...
    prolongation(cGc(2:nc+1,2:nc+1));
  
  % Post-smoothing.
  u = smoothJacobiDamped(f,u,hf,MGParam.m2,MGParam.omega);

  % Alternative post-smoothing options:
  % u = smoothJacobi(f,u,hf,MGParam.m2);
  % u = smoothGaussSeidel(f,u,hf,MGParam.m2,'Backward');
end

end
