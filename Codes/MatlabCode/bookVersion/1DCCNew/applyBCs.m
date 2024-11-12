function u = applyBCs(u)

  nfPlusGhostLayers = length(u);
  nf = nfPlusGhostLayers-2;
%
% Homogeneous Dirichlet BC:
%  u(   1) = -u(   2);
%  u(nf+2) = -u(nf+1);
%
% Homogeneous Neumann BC:
  u(   1) = u(  2);
  u(nf+2) = u(nf+1);

end