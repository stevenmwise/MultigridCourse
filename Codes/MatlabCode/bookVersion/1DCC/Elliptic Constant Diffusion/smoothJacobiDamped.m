function u = smoothJacobiDamped(f,u,h,m,omega)

nPlusGhostLayers = length(u);
n = nPlusGhostLayers-2;
h2 = h*h;

u = applyBCs(u);
for sweep = 1:m
  uOld = u;
  for i = 2:n+1
      
    % Old iteration neighbor values.
    uE = uOld(i+1);
    uW = uOld(i-1);

    % Diagonal entry = 2/h^2 + 1.
    diagC = 2.0/h2+1.0;

    % Right-hand side.
    rhs = f(i-1)+(uE+uW)/h2;

    % Pure Jacobi value.
    z = rhs/diagC;

    % Damped Jacobi update.
    u(i) = (1.0-omega)*uOld(i)+omega*z;
  end
  u = applyBCs(u);
end

end