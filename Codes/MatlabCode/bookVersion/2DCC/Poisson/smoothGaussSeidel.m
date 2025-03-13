function u = smoothGaussSeidel(f,u,hf,m,sweepDirection)

[nFull,~] = size(u);
n = nFull-2;
hf2 = hf^2;

if strcmp(sweepDirection,'Forward')
  xSweepIndex=[2:n+1];
  ySweepIndex=[2:n+1];
elseif strcmp(sweepDirection,'Backward')
  xSweepIndex=[n+1:-1:2];
  ySweepIndex=[n+1:-1:2];
else
  xSweepIndex=[2:n+1];
  ySweepIndex=[2:n+1];
end

for sweep = 1:m
  for i = xSweepIndex
    for j = ySweepIndex
        
      % Neighbors (using newest values available).
      uE = u(i+1,j);
      uW = u(i-1,j);
      uN = u(i,j+1);
      uS = u(i,j-1);

      % Diagonal entry = 4/h^2.
      diagC = 4.0/hf2;

      % Right-hand side.
      rhs = f(i-1,j-1)+(uE+uW+uN+uS)/hf2;

      u(i,j) = rhs/diagC;
    end
  end
  u = applyBCs(u);
end

end