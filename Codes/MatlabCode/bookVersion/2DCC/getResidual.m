function residualResult = getResidual(u,f,hf)
  u = applyBCs(u);
  residualResult = f-FDOperator(u,hf);
end