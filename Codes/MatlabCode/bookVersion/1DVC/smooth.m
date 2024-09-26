function [u] = smooth(f,u,m,omega)

if m == 0
  return
end

nl = length(f);
hl = 1/(nl+1);
uSmooth = zeros(nl,1);
omegaPrime = 1.0-omega;

for j = 1:m

  uSmooth(1) = (hl*f(1)+u(2))/2.0;
  uSmooth(2:nl-1)= (hl*f(2:nl-1)+u(1:nl-2)+u(3:nl))/2.0;
  uSmooth(nl) = (hl*f(nl)+u(nl-1))/2.0;

  uSmooth(1:nl) = omega*uSmooth(1:nl)+omegaPrime*u(1:nl);

  u = uSmooth;
  
end

end