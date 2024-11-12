function [u,errVals,kStop] = multiGridSolver(f,MGParam,uExact)

L = MGParam.L;
nL = MGParam.nL;
kMax = MGParam.kMax;
tol = MGParam.tol;
errVals = zeros(kMax,3);

nl = nL;
for k = L-1:-1:0
  if mod(nl+1,2) ~= 0
    disp('Coarsening error: L is too large.')
    u = zeros(nL,1);
    rate = -1;
    return
  end
  nl = (nl+1)/2-1;
end

SEED = 1234;
rng(SEED);
u = rand(nL,1)-0.5;
rate = 1.0;

hL = 1/(nL+1);
x = linspace(hL,1.0-hL,nL)';

k = 0;
err = tol+1.0;

while ((k<kMax) & (err > tol))
  
  uo = u;
  u = MGOperator(f,u,L,MGParam);

  err = norm(uExact-u,2)
  errVals(k+1,2) = norm(u-uo,2);
  errVals(k+1,3) = err;
  k = k+1
  
end

kStop = k;

end