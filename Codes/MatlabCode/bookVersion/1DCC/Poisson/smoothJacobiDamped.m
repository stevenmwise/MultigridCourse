function u = smoothJacobiDamped(f,u,h,m,omega)

nPlusGhostLayers = length(u);
n = nPlusGhostLayers-2;
h2 = h*h;

u = applyBCs(u);

for sweep = 1:m

  % Temporary variable to store the updated values.
  z = zeros(nPlusGhostLayers,1);
  
  for i = 2:n+1
      
    % Compute the right-hand side.
    rhs = f(i-1)+(u(i+1)+u(i-1))/h2;
    
    % Compute the diagonal element.
    diagC = 2.0/h2;
    
    z(i) = rhs/diagC;
  end
  
  u(2:n+1) = omega*z(2:n+1)+(1.0-omega)*u(2:n+1);
  
  u = applyBCs(u);
end

end
