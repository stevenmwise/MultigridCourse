function [w,errVals,kStop] = SOR(f,w,wExact,iterParam)

n = iterParam.n;
h = iterParam.h;
h2 = iterParam.h2;
kMax = iterParam.kMax;
tol = iterParam.tol;
sig = cos(pi*h);
omega = 2/(1+sqrt(1-sig*sig))
omOmega = 1.0-omega;

errVals = zeros(kMax,1);

k = 0;
err = tol+1.0;
%
% SOR loop.
while ((k < kMax) && (err > tol))
  
  tmp = (h2*f(1)+w(2))/2.0;
  w(1) = omega*tmp+omOmega*w(1);
  for i = 2:n-1
    tmp = (h2*f(i)+w(i+1)+w(i-1))/2.0;
    w(i) = omega*tmp+omOmega*w(i);
  end
  tmp = (h2*f(n)+w(n-1))/2.0;
  w(n) = omega*tmp+omOmega*w(n);
  
  err = norm(wExact-w,2);
  k = k+1;
  errVals(k) = err;
end
kStop = k;

end