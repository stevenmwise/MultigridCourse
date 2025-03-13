function u = smoothGaussSeidel(f,u,h,m,sweepDirection)

nPlusGhostLayers = length(u);
n = nPlusGhostLayers-2;
h2 = h*h;

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
    u(i) = rhs/diagC;
  end
  
  u = applyBCs(u);
end

end