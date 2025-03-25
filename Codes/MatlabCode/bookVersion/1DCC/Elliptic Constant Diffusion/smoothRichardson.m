function u = smoothRichardson(f,u,h,m,omega)

nPlusGhostLayers = length(u);
n = nPlusGhostLayers-2;
h2 = h*h;

u = applyBCs(u);

for sweep = 1:m
  r = getResidual(u,f,h);
  u(2:n+1) = u(2:n+1)+omega*(h2/2)*r;
  u = applyBCs(u);
end

end