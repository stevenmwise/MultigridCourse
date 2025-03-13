function normResult = normScaledL2(u)
  nf = size(u);
  normResult = sqrt(sum(sum(u.*u))/(nf(1)*nf(2)));
end