function [u] = smooth(f,u,m,omega)

if m == 0
  return
end

nl = length(f);
hl = 1.0/nl;
uSmooth = zeros(nl+2,1);
omegaPrime = 1.0-omega;

% Enforce Dirichlet boundary conditions.
u(   1) = -u(   2);
u(nl+2) = -u(nl+1);

for j = 1:m

  uSmooth(2:nl+1)= (hl*hl*f(1:nl)+u(1:nl)+u(3:nl+2))/2.0;
  uSmooth(2:nl+1) = omega*uSmooth(2:nl+1)+omegaPrime*u(2:nl+1);

  u(2:nl+1) = uSmooth(2:nl+1);

% Enforce Dirichlet boundary conditions.
  u(   1) = -u(   2);
  u(nl+2) = -u(nl+1);
  
end

end