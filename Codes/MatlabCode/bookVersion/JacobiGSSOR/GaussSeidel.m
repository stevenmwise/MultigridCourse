function [w,errVals,kStop] = GaussSeidel(f,w,wExact,iterParam)

n = iterParam.n;
h2 = iterParam.h2;
kMax = iterParam.kMax;
tol = iterParam.tol;

errVals = zeros(kMax,1);

k = 0;
err = tol+1.0;
%
% Gauss-Seidel loop.
while ((k < kMax) && (err > tol))
  
  w(1) = (h2*f(1)+w(2))/2.0; 
  for i = 2:n-1
    w(i) = (h2*f(i)+w(i+1)+w(i-1))/2.0;
  end
  w(n) = (h2*f(n)+w(n-1))/2.0; 
  
  err = norm(wExact-w,2);
  k = k+1;
  errVals(k) = err;
end
kStop = k;

end