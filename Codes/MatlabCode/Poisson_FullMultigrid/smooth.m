function u = smooth(u,f,h,m,omega)

  nPlusGhostLayers = length(u);
  n = nPlusGhostLayers-2;
  h2 = h*h;

  omegaPrime = 1.0 - omega;

  u = applyBCs(u);

  for k = 1:m
    % Temporary variable to store the updated values.
    uNew = zeros(nPlusGhostLayers,1); 

    for i = 2:n+1
      uNew(i) = (h2*f(i-1)+u(i+1)+u(i-1))/(2.0); 
    end

    u(2:n+1) = omega * uNew(2:n+1) + omegaPrime * u(2:n+1);

    u = applyBCs(u);

  end
end
