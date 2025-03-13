function operatorResult = FDOperator(u,hf)

nPlusGhostLayers = length(u);

nf = nPlusGhostLayers-2;
hf2 = hf*hf;

% Discrete -Laplacian(u).
operatorResult = (-u(3:nf+2)-u(1:nf)+2.0*u(2:nf+1))/hf2;

end