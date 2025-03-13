function u = smoothJacobi(f,u,h,m)

nPlusGhostLayers = length(u);
n = nPlusGhostLayers-2;
h2 = h*h;
u = applyBCs(u);
for sweep = 1:m
  uNew = zeros(nPlusGhostLayers,1);
  
  % Loop over interior points.
  for i = 2:n+1
    uNew(i) = (h2*f(i-1)+u(i+1)+u(i-1))/(2.0+h2);
  end
  
  u(2:n+1) = uNew(2:n+1);
  u = applyBCs(u);
end
end
