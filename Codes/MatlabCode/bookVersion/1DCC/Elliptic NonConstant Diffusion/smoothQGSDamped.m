function u = smoothQGSDamped(f,u,h,D,m,omega,sweepDirection)

n = length(u)-2;
h2 = h*h;
omegaPrime = 1-omega;

if strcmp(sweepDirection,'Forward')
  sweepIndex=[1:n];
elseif strcmp(sweepDirection,'Backward')
  sweepIndex=[n:-1:1];
else
  sweepIndex=[1:n];
end

for sweep = 1:m
  for i = sweepIndex
    iU = i+1;
    
    % Compute the diagonal coefficient.
    diagC = D(i)+D(i+1)+h2;
    
    % Compute the right-hand side.
    rhs = h2*f(i)+D(i+1)*u(iU+1)+D(i)*u(iU-1);
    
    % Update using Gauss-Seidel.
    u(iU) = omega*rhs/diagC+omegaPrime*u(iU);
  end
  
  u = applyBCs(u);
end
end
