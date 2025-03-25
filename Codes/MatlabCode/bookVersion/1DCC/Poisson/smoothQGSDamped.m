function u = smoothQGSDamped(f,u,h,m,omega,sweepDirection)

nPlusGhostLayers = length(u);
n = nPlusGhostLayers-2;
h2 = h*h;
omegaPrime = 1-omega;

if strcmp(sweepDirection,'Forward')
  sweepIndex=[2:n+1];
elseif strcmp(sweepDirection,'Backward')
  sweepIndex=[n+1:-1:2];
else
  sweepIndex=[2:n+1];
end

u = applyBCs(u);

for sweep = 1:m
  for i = sweepIndex

    % Compute the right-hand side.
    rhs = f(i-1)+(u(i+1)+u(i-1))/h2;
    
    % Compute the diagonal element.
    diagC = 2.0/h2;
    
    % Update immediately (Gauss-Seidel).
    u(i) = omega*rhs/diagC+omegaPrime*u(i);
  end
  
  u = applyBCs(u);
end

end