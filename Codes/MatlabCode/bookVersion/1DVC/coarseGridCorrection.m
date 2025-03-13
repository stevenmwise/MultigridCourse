function cgc = coarseGridCorrection(res)

n1 = length(res);
cgc = zeros(n1,1);
n0 = (n1+1)/2-1;
h0 = 1/(n0+1);

% Assemble coarse grid matrix.
A0 = gallery('tridiag',n0,-1.0/h0,2.0/h0,-1.0/h0);

% Solve coarse grid problem and interpolate.
cgc = prolongation(A0\restriction(res));

end