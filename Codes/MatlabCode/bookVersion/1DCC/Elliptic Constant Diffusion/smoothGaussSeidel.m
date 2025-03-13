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

for sweep = 1:m
  for i = sweepIndex
      
    % Neighbors (using newest values available).
    uE = u(i+1);
    uW = u(i-1);

    % Diagonal entry = 2/h^2 + 1.
    diagC = 2.0/h2+1.0;

    % Right-hand side.
    rhs = f(i-1)+(uE+uW)/h2;

    u(i) = rhs/diagC;
  end
  u = applyBCs(u);
end

end