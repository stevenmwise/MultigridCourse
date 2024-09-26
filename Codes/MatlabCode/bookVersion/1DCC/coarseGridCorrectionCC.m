function [cgc] = coarseGridCorrectionCC(res)

n1 = length(res);
cgc = zeros(n1,1);
n0 = n1/2;
h0 = 1/n0;
A0 = gallery('tridiag',n0,-1.0,2.0,-1.0);
A0( 1, 1) = 3.0;
A0(n0,n0) = 3.0;

cgc = prolongationCC(A0\restrictionCC(res));

end