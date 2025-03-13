function [u] = MGOperator(g,u,level,MGParam)

nl = length(g);
hl = 1/(nl+1);

if level == 0
  Al = gallery('tridiag',nl,-1.0/hl,2.0/hl,-1.0/hl);
  % Direct solve at the coarsest level.
  u = Al\g;
else
  % Pre-smoothing.
  u = smooth(g,u,MGParam.m1,MGParam.omega); 
  cGr = restriction(getResidual(g,u));
  nlm1 = (nl+1)/2-1;
  cGc = zeros(nlm1,1);
  % Compute the coarse-grid correction with recursive MG call.
  for sig = 1:MGParam.p                     
    cGc = MGOperator(cGr,cGc,level-1,MGParam);
  end
  % Correct the solution using coarse-grid correction.
  u = u+prolongation(cGc);
  % Post-smoothing.
  u = smooth(g,u,MGParam.m2,MGParam.omega); 
end

end