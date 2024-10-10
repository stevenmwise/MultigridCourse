function operatorResult = FDOperator(u,hf)
  nfPlusGhostLayers = length(u);
  nf = nfPlusGhostLayers-2;
  hf2 = hf*hf;
%
% -dicreteLaplacian(u) + u: 
%
  operatorResult = (-     u(3:nf+2) ...
                    -     u(1:nf  ) ...
                    + 2.0*u(2:nf+1))/hf2 ...
                    +     u(2:nf+1);
end