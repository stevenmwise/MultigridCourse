function normResult = normScaledL2(u)
  nf = length(u);
  normResult = sqrt(sum(u.*u)/nf);
end