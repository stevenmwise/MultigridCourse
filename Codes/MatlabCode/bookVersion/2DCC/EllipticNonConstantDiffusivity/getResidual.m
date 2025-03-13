function residualResult = getResidual(u,f,hf,MGParam,D_EW,D_NS)

% Apply the boundary conditions
u = applyBCs(u);

% Residual = f - ( -div(D grad(u)) + u )
Au = FDOperator(u,hf,MGParam,D_EW,D_NS);
residualResult = f-Au;

end
