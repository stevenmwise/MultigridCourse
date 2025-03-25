function operatorResult = FDOperator(u,DeW,DnS,hf,MGParam)

nfPlusGhostLayers = size(u);
nf = nfPlusGhostLayers-2;

nx = nf(1);
ny = nf(2);
hf2 = hf*hf;
%
% Vectorized Code:
%
% - div(D grad(u)) + u 
%
% Compute east-west and north-south fluxes:
feW(1:nx+1,1:ny) ...
  = DeW(1:nx+1,1:ny).*(u(2:nx+2,2:ny+1)-u(1:nx+1,2:ny+1));
fnS(1:nx,1:ny+1) ...
  = DnS(1:nx,1:ny+1).*(u(2:nx+1,2:ny+2)-u(2:nx+1,1:ny+1));

% Compute the divergence of the fluxes and add reaction term:
operatorResult(1:nx,1:ny) ...
  = -(feW(2:nx+1,1:ny)-feW(1:nx,1:ny) ...
   +  fnS(1:nx,2:ny+1)-fnS(1:nx,1:ny))/hf2 ...
   + u(2:nx+1,2:ny+1);

end