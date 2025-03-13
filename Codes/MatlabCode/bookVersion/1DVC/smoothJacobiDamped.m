function u = smooth(f,u,m,omega)

if m == 0
  return
end

nl = length(f);
hl = 1/(nl+1);
z = zeros(nl,1);

for j = 1:m

  % Compute smoothing step.
  z(1) = (hl*f(1)+u(2))/2.0;
  z(2:nl-1) = (hl*f(2:nl-1)+u(1:nl-2)+u(3:nl))/2.0;
  z(nl) = (hl*f(nl)+u(nl-1))/2.0;

  % Apply damping.
  z(1:nl) = omega*z(1:nl)+(1.0-omega)*u(1:nl);

  u = z;
end

end