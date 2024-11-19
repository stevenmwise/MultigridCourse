function u = smooth(u,f,h,m,omega)

  nPlusGhostLayers = length(u);
  n = nPlusGhostLayers-2;
  h2 = h*h;
    
  u = applyBCs(u);

  for k = 1:m
    for i = 2:n+1
      u(i) = (h2*f(i-1)+u(i+1)+u(i-1))/(2.0+h2); 
    end

    u = applyBCs(u);

  end
end