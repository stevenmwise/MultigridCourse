function u = smoothGaussSeidel(f,u,hf,D_EW,D_NS,m,sweepDirection)

[nfPlusGhostX,nfPlusGhostY] = size(u);

% Number of interior cells in x.
nx = nfPlusGhostX-2;

% Number of interior cells in y.
ny = nfPlusGhostY-2;

hf2 = hf*hf;

% Vectorized precomputation of diagonal coefficients.
diagCoeff = (D_EW(2:nx+1,1:ny)+D_EW(1:nx,1:ny)+ ...
    D_NS(1:nx,2:ny+1)+D_NS(1:nx,1:ny))/hf2+1.0;

if strcmp(sweepDirection,'Forward')
  xSweepIndex = [1:nx];
  ySweepIndex = [1:ny];
elseif strcmp(sweepDirection,'Backward')
  xSweepIndex = [nx:-1:1];
  ySweepIndex = [ny:-1:1];
else
  xSweepIndex = [1:nx];
  ySweepIndex = [1:ny];
end

for sweep = 1:m

  % Loop over interior points
  for i = xSweepIndex
    for j = ySweepIndex
      iU = i+1;
      jU = j+1;
            
      % Fix the coefficient calculations.
      De = D_EW(i+1,j);
      Dw = D_EW(i,j);
      Dn = D_NS(i,j+1);
      Ds = D_NS(i,j);
            
      % Update to match the corrected operator form
      rhs = f(i,j)+(De*u(iU+1,jU)+Dw*u(iU-1,jU)+ ...
          Dn*u(iU,jU+1)+Ds*u(iU,jU-1))/hf2;
            
      u(iU,jU) = rhs/diagCoeff(i,j);
    end
  end
  u = applyBCs(u);
end

end