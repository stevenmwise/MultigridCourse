function operatorResult = FDOperator(u,hf,MGParam,D)
n = length(u)-2;

h2 = hf^2;
operatorResult = zeros(n,1);

for i = 1:n
  iU = i+1;
  
  % Compute the flux difference.
  operatorResult(i) = (-(D(i+1)*(u(iU+1)-u(iU))-...
    D(i)*(u(iU)-u(iU-1)))/h2)+u(iU);
end

end
