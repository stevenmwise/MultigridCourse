function operatorResult = FDOperator(u,hf,MGParam,D_EW,D_NS)

nfPlusGhostLayers = size(u);
nf = nfPlusGhostLayers-2;

nx = nf(1);
ny = nf(2);
hf2 = hf*hf;

operatorResult = zeros(nx,ny);

for i = 1:nx
  for j = 1:ny
    iU = i+1;
    jU = j+1;
        
    termEW = -(D_EW(i+1,j)*(u(iU+1,jU)-u(iU,jU))- ...
        D_EW(i,j)*(u(iU,jU)-u(iU-1,jU)));
    termNS = -(D_NS(i,j+1)*(u(iU,jU+1) - u(iU,jU))- ...
        D_NS(i,j)*(u(iU,jU)-u(iU,jU-1)));
        
    operatorResult(i,j) = (termEW+termNS)/hf2+u(iU,jU);
    
  end
end

end
