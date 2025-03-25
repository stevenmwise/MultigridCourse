function residualResult = getResidual(u,f,DeW,DnS,hf,MGParam)

% Apply the boundary conditions
u = applyBCs(u);

residualResult = f-FDOperator(u,DeW,DnS,hf,MGParam);

end
