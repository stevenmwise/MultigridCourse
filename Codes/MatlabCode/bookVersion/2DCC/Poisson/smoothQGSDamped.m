function u = smoothQGSDamped(f,u,hf,m,omega,sweepDirection)

[nFull,~] = size(u);
n = nFull-2;
hf2 = hf*hf;
omegaPrime = 1.0-omega;

switch sweepDirection
  case 'Forward'
    xSweepIndex=[2:n+1];
    ySweepIndex=[2:n+1];
  case 'Backward'
    xSweepIndex=[n+1:-1:2];
    ySweepIndex=[n+1:-1:2];
  otherwise
    xSweepIndex=[2:n+1];
    ySweepIndex=[2:n+1];
end

for sweep = 1:m
  for i = xSweepIndex
    for j = ySweepIndex

      z = (hf2*f(i-1,j-1)+u(i+1,j)+u(i-1,j)...
          +u(i,j+1)+u(i,j-1))/4.0;

      u(i,j) = omega*z+omegaPrime*u(i,j);

    end
  end
  u = applyBCs(u);
end

end