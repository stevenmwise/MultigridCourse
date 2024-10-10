function [u] = MGOperator(g,u,level,MGParam)

nl = length(g);
hl = 1/(nl+1);

if level == 0
  Al = gallery('tridiag',nl,-1.0/hl,2.0/hl,-1.0/hl);
  u = Al\g;
else
  u = smooth(g,u,MGParam.m1,MGParam.omega);
  res = restriction(getResidual(g,u));
  nlm1 = (nl+1)/2-1;
  cGc = zeros(nlm1,1);
  for sig = 1:MGParam.p
    cGc = MGOperator(res,cGc,level-1,MGParam);
  end
  u = u+prolongation(cGc);
  u = smooth(g,u,MGParam.m2,MGParam.omega);
end

end