function u = smoothJacobi(f,u,h,D,m)

n = length(u)-2;
h2 = h*h;

for sweep = 1:m

  % Copy current iterate.
  uOld = u;
  
  for i = 1:n
    iU = i+1;
    
    % Compute the diagonal coefficient.
    diagC = (D(i)+D(i+1))/h2+1.0;
    
    % Compute the right-hand side using uOld.
    rhs = f(i)+(D(i+1)*uOld(iU+1)+D(i)*uOld(iU-1))/h2;
    u(iU) = rhs/diagC;
  end
  
  u = applyBCs(u);
end
end
