function [cgc] = coarseGridCorrection(res)

n1 = length(res);
cgc = zeros(n1,1);
n0 = (n1+1)/2-1;
h0 = 1/(n0+1);
A0 = gallery('tridiag',n0,-1.0/h0,2.0/h0,-1.0/h0);

cgc = prolongation(A0\restriction(res));

end