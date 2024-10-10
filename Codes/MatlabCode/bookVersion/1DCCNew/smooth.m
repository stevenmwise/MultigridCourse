function u = smooth(u,f,hf,m,omega)

  nfPlusGhostLayers = size(u);
  nf = nfPlusGhostLayers-2;
  hf2 = hf*hf;
    
  u = applyBCs(u);

  for k = 1:m
    for i = 2:nf(1)+1
      u(i) = (hf2*f(i-1)+u(i+1)+u(i-1))/(2.0+hf2); 
    end

    u = applyBCs(u);

  end
end