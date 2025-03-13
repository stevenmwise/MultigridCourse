function residualResult = getResidual(u,f,hf,MGParam,D)

u = applyBCs(u);
residualResult = f-FDOperator(u,hf,MGParam,D);
end
