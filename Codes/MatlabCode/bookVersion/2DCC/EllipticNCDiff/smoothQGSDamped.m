function u = smoothQGSDamped(f,DeW,DnS,rC,u,hf,m,...
  omega,sweepDirection)

nfPlusGhostLayers = size(u);
nf = nfPlusGhostLayers-2;

hf2 = hf*hf;
omegaPrime = 1.0-omega;

switch sweepDirection
  case 'fwd'
    xSweepIndex = [1:nf(1)];
    ySweepIndex = [1:nf(2)];
  otherwise
    xSweepIndex = [nf(1):-1:1];
    ySweepIndex = [nf(2):-1:1];
end

u = applyBCs(u);

for sweep = 1:m
  for j = ySweepIndex
    for i = xSweepIndex

      tmp = (hf2*f(i,j)+u(i+2,j+1)*DeW(i+1,j  ) ...
          +             u(i  ,j+1)*DeW(i  ,j  ) ...
          +             u(i+1,j+2)*DnS(i  ,j+1) ...
          +             u(i+1,j  )*DnS(i  ,j  )) ...
          / (DeW(i+1,j)+DeW(i,j)+DnS(i,j+1)+DnS(i,j) ...
          + rC*hf2);
      u(i+1,j+1) = omega*tmp+omegaPrime*u(i+1,j+1);

    end
  end
  u = applyBCs(u);

end

end