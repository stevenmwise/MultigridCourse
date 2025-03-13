function u = applyBCs(u)
  nfPlusGhostLayers = size(u);
  nf = nfPlusGhostLayers-2;
%
% Homogeneous Dirichlet BC:
%
%  u(2:nf(1)+1,        1) = -u(2:nf(1)+1,        2);
%  u(2:nf(1)+1,  nf(2)+2) = -u(2:nf(1)+1,  nf(2)+1);
%  u(        1,1:nf(2)+2) = -u(        2,1:nf(2)+2);
%  u(  nf(1)+2,1:nf(2)+2) = -u(  nf(1)+1,1:nf(2)+2);
%
% Homogeneous Neumann BC:
%
  u(2:nf(1)+1,        1) = u(2:nf(1)+1,        2);
  u(2:nf(1)+1,  nf(2)+2) = u(2:nf(1)+1,  nf(2)+1);
  u(        1,1:nf(2)+2) = u(        2,1:nf(2)+2);
  u(  nf(1)+2,1:nf(2)+2) = u(  nf(1)+1,1:nf(2)+2);
%
end