function u = smoothQJacDamped(f,u,hf,m,omega)

[n1Full,n2Full] = size(u);
n1 = n1Full-2;
n2 = n2Full-2;
hf2 = hf^2;

u = applyBCs(u);
for sweep = 1:m
  uOld = u;
  for i = 2:n1+1
    for j = 2:n2+1
        
      % Old iteration neighbor values.
      uE = uOld(i+1,j);
      uW = uOld(i-1,j);
      uN = uOld(i,j+1);
      uS = uOld(i,j-1);

      % Diagonal entry = 4/h^2 + 1.
      diagC = 4.0/hf2+1.0;

      % Right-hand side.
      rhs = f(i-1,j-1)+(uE+uW+uN+uS)/hf2;

      % Pure Jacobi value.
      z = rhs/diagC;

      % Damped Jacobi update.
      u(i,j) = (1-omega)*uOld(i,j)+omega*z;
    end
  end
  u = applyBCs(u);
end

end