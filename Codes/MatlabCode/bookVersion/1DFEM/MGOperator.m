function u = MGOperator(level,f,u,MGParam)

nl = length(f);
hl = 1/(nl+1);

if level == 0
    
  % Direct solve at the coarsest level.
  Al = gallery('tridiag',nl,-1.0/hl,2.0/hl,-1.0/hl);
  u = Al\f;
else

  % Pre-smoothing.
  u = smoothJacobiDamped(f,u,MGParam.m1,MGParam.omega);

  % Alternative pre-smoothing options:
  % u = smoothJacobi(f,u,MGParam.m1);
  % u = smoothGaussSeidel(f,u,MGParam.m1,'Forward');

  % Get residual and restrict to coarse grid.
  cGr = restriction(getResidual(f,u));
  nlm1 = (nl+1)/2-1;
  cGc = zeros(nlm1,1);

  % Compute the coarse-grid correction with recursive MG call.
  for sig = 1:MGParam.p
    cGc = MGOperator(level-1,cGr,cGc,MGParam);
  end

  % Correct the solution using coarse-grid correction.
  u = u+prolongation(cGc);

  % Post-smoothing.
  u = smoothJacobiDamped(f,u,MGParam.m2,MGParam.omega);

  % Alternative post-smoothing options:
  % u = smoothJacobi(f,u,MGParam.m2);
  % u = smoothGaussSeidel(f,u,MGParam.m2,'Backward');
end

end