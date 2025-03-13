function u = smoothJacobiDamped(f,u,h,D,numSweeps,omega)

n = length(u)-2;
h2 = h*h;

for sweep = 1:numSweeps

  % Copy current iterate.
  uOld = u;
  
  for i = 1:n
    iU = i+1;
    
    % Compute the diagonal coefficient.
    diagC = (D(i)+D(i+1))/h2+1.0;
    
    % Compute the right-hand side using uOld.
    rhs = f(i)+(D(i+1)*uOld(iU+1)+D(i)*uOld(iU-1))/h2;
    
    % Compute pure Jacobi value.
    z = rhs/diagC;
    
    % Damped Jacobi update.
    u(iU) = omega*z+(1.0-omega)*uOld(iU);
  end
  
  u = applyBCs(u);
end

end
