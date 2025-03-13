function u = smoothJacobi(f,u,hf,D_EW,D_NS,m)

[nfPlusGhostX,nfPlusGhostY] = size(u);

% Number of interior cells in x.
nx = nfPlusGhostX-2; 

% Number of interior cells in y.
ny = nfPlusGhostY-2;  

hf2 = hf^2;

% Vectorized precomputation of diagonal coefficients.
diagCoeff = (D_EW(2:nx+1,1:ny)+D_EW(1:nx,1:ny)+ ...
    D_NS(1:nx,2:ny+1)+D_NS(1:nx,1:ny))/hf2+1.0;

u = applyBCs(u);
for sweep = 1:m
  uOld = u;
  
  % Loop over interior points
  for j = 1:ny
    for i = 1:nx
      iU = i+1;
      jU = j+1;

      % Neighbors (from the old iteration).
      uE = uOld(iU+1,jU);  
      uW = uOld(iU-1,jU);
      uN = uOld(iU,jU+1);
      uS = uOld(iU,jU-1);
            
      % Diffusion coefficients for cell (i,j).
      De = D_EW(i+1,j);  
      Dw = D_EW(i,j);  
      Dn = D_NS(i,j+1); 
      Ds = D_NS(i,j); 
            
      % Use precomputed diagonal coefficient.
      rhs = f(i,j)+(De*uE+Dw*uW+Dn*uN+Ds*uS)/hf2;
      u(iU,jU) = rhs/diagCoeff(i,j);
    end
  end
  u = applyBCs(u);
end

end