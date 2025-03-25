function u = smoothQJacDamped(f,DeW,DnS,rC,u,hf,m,omega)

nfPlusGhostLayers = size(u);
nf = nfPlusGhostLayers-2;

z = zeros(nf(1),nf(2));

hf2 = hf*hf;
omegaPrime = 1.0-omega;

u = applyBCs(u);

for sweep = 1:m
  for j = 1:nf(2)
    for i = 1:nf(1)

      z(i,j) = (hf2*f(i,j)+u(i+2,j+1)*DeW(i+1,j  ) ...
             +             u(i  ,j+1)*DeW(i  ,j  ) ...
             +             u(i+1,j+2)*DnS(i  ,j+1) ...
             +             u(i+1,j  )*DnS(i  ,j  )) ...
             / (DeW(i+1,j)+DeW(i,j)+DnS(i,j+1)+DnS(i,j) ...
             + rC*hf2);

    end
  end

  u(2:nf(1)+1,2:nf(2)+1) = omega*z(1:nf(1),1:nf(2)) ...
    + omegaPrime*u(2:nf(1)+1,2:nf(2)+1);
  u = applyBCs(u);

end

end