function u = QGSDamped(f,u,hf,m,omega,sweepDirection)

[n1Full,n2Full] = size(u);
n1 = n1Full-2;
n2 = n2Full-2;
hf2 = hf^2;
omegaPrime = 1-omega;

if strcmp(sweepDirection,'Forward')
  xSweepIndex=[2:n1+1];
  ySweepIndex=[2:n2+1];
elseif strcmp(sweepDirection,'Backward')
  xSweepIndex=[n1+1:-1:2];
  ySweepIndex=[n2+1:-1:2];
else
  xSweepIndex=[2:n1+1];
  ySweepIndex=[2:n2+1];
end

for sweep = 1:m
  for i = xSweepIndex
    for j = ySweepIndex
        
      % Neighbors (using newest values available).
      uE = u(i+1,j);
      uW = u(i-1,j);
      uN = u(i,j+1);
      uS = u(i,j-1);

      % Diagonal entry = 4/h^2 + 1.
      diagC = 4.0/hf2+1.0;

      % Right-hand side.
      rhs = f(i-1,j-1)+(uE+uW+uN+uS)/hf2;

      u(i,j) = omega*rhs/diagC+omegaPrime*u(i,j);
    end
  end
  u = applyBCs(u);
end

end