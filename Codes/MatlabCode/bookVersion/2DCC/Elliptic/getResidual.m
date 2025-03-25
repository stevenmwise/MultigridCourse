function residualResult = getResidual(u,f,hf,MGParam)

u = applyBCs(u);
residualResult = f-FDOperator(u,hf,MGParam);

end
