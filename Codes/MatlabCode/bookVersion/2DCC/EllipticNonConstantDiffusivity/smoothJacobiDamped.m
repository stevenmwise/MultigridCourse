function u = smoothJacobiDamped(f,u,hf,D_EW,D_NS,m,omega)

[nfPlusGhostX,nfPlusGhostY] = size(u);

% Number of interior cells in x.
nx = nfPlusGhostX-2;

% Number of interior cells in y.
ny = nfPlusGhostY-2;

hf2 = hf*hf;

% Vectorized precomputation of diagonal coefficients.
diagCoeff = (D_EW(2:nx+1,1:ny)+D_EW(1:nx,1:ny)+...
  D_NS(1:nx,2:ny+1)+D_NS(1:nx,1:ny))/hf2+1.0;

for sweep = 1:m
  uOld = u;
  
  % Loop over interior points
  for j = 1:ny
    for i = 1:nx
      iU = i+1;
      jU = j+1;

      % Old iteration neighbor values:
      uE = uOld(iU+1,jU);
      uW = uOld(iU-1,jU);
      uN = uOld(iU,jU+1);
      uS = uOld(iU,jU-1);

      % Diffusion coefficients:
      De = D_EW(i+1,j);
      Dw = D_EW(i,j);
      Dn = D_NS(i,j+1);
      Ds = D_NS(i,j);

      % Right-hand side.
      rhs = f(i,j)+(De*uE+Dw*uW+Dn*uN+Ds*uS)/hf2;

      % Pure Jacobi value with precomputed diagonal.
      z = rhs/diagCoeff(i,j);

      % Damped Jacobi update.
      u(iU,jU) = (1.0-omega)*u(iU,jU)+omega*z;
    end
  end

  % Re-apply boundary conditions each sweep.
  u = applyBCs(u);
end

end