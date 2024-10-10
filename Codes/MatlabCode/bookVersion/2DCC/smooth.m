function u = smooth(u,f,hf,m,omega)

  nfPlusGhostLayers = size(u);
  nf = nfPlusGhostLayers-2;
  hf2 = hf*hf;
    
  u = applyBCs(u);

  for k = 1:m
    for j = 2:nf(2)+1
      for i = 2:nf(1)+1
        u(i,j) = (hf2*f(i-1,j-1) ...
               + u(i+1,j)+u(i-1,j) ...
               + u(i,j+1)+u(i,j-1))/(4.0+hf2);
      end
    end

    u = applyBCs(u);

  end
end