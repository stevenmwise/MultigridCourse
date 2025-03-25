%% Initialization and Parameter Setting

% Clear command line window and workspace.
clear; clc; close all;

% Tolerance of residual to stop multigrid iterations.
tol = 1.0e-10;

% Number of multigrid levels.
L = 5;

% Number of interior cells in x,y at finest level.
nL = [2^L-1,2^L-1];

% Define the domain of the equation.
xLower = [0.0,0.0];
xUpper = [1.0,1.0];

% Define the exact solution and source term.
syms x y;
px = [xUpper(1)-xLower(1),xUpper(2)-xLower(2)];
u_exact = exp(sin(4.0*pi*x/px(1))*sin(4.0*pi*y/px(2)))-1;
f = -diff(diff(u_exact,x),x)-diff(diff(u_exact,y),y);

% Convert to MATLAB functions.
f_func = matlabFunction(f,'Vars',{[x,y]});
u_exact_func = matlabFunction(u_exact,'Vars',{[x,y]});

%% Setting up the Mesh and Boundary Conditions

% Initial coarse mesh.
[node,elem] =...
  squaremesh([xLower(1),xUpper(1),xLower(2),xUpper(2)],px(1));

% Refine the mesh.
for i = 1:L+1
  [node,elem] = uniformrefine(node,elem);
end

% Set Dirichlet boundary conditions on all boundaries.
bdFlag = setboundary(node,elem,'Dirichlet','all');

%% Setting up the Problem and Solver Parameters

% Problem definition.

% Homogeneous Dirichlet Condition.
pde.g_D = 0;

% Forcing function.
pde.f = f_func;

% Exact solution for error calculations.
pde.exactu = u_exact_func;

% Multigrid solver parameters.
option.solver                 = 'mg';
option.mgoption.solver        = 'VCYCLE';
option.mgoption.smoothingstep = 3;
option.mgoption.printlevel    = 2;
option.tol                    = tol;
option.maxIt                  = 100;

%% Calling the iFEM Solver

[soln,eqn,info] = Poisson(node,elem,bdFlag,pde,option);

%% Plotting the Numerical Solution

figure(1)
showresult(node,elem,soln.u);
title('Numerical Solution of the Poisson Equation');
xlabel('x');
ylabel('y');
zlabel('u');
colorbar;

%% Error Analysis and Convergence Plots

% Extract error values.
errVals = info.error;
kStop = size(errVals,1);

figure(2)
clf

% Plot available error quantities.
semilogy(errVals(1:kStop,1),'rs','LineWidth',1.5)
hold on
semilogy(errVals(1:kStop,2),'kd','LineWidth',1.5)

% Log-linear fit for convergence rate.
rate = 1.0;
if kStop>=4
  kv = kStop-3:kStop;
  le = log(errVals(kStop-3:kStop,2));
  p1 = polyfit(kv,le,1);
  rate = exp(p1(1));
  p1k = polyval(p1,1:kStop);
  semilogy(1:kStop,exp(p1k),'-k')
end

% Labels and title.
xlabel('$k$','Interpreter','latex');
title('Multigrid Iteration Errors','Interpreter','latex');
legend('$\left\|{\bf u}_L^k-{\bf u}_L^{k-1} \right\|_L$',...
  '$\left\|{\bf r}_L^k\right\|_L$',...
  'log-linear fit','Interpreter','latex');

% Add annotations for parameters.
text(1.5,2048*option.tol,strcat(...
  '$\gamma_{\rm comp} =\hspace{.1cm}$',...
  num2str(rate,'%10.5e')),'FontSize',14,'Interpreter','latex');
text(1.5,512*option.tol,strcat('$m_1 =\hspace{.1cm}$',...
  num2str(option.mgoption.smoothingstep)),'FontSize',14,...
  'Interpreter','latex');
text(1.5,128*option.tol,strcat('$m_2 =\hspace{.1cm}$',...
  num2str(option.mgoption.smoothingstep)),'FontSize',14,...
  'Interpreter','latex');

% Save figure as PDF.
printstr = strcat('Err_nL_',num2str(nL(1)),'_m1_',...
  num2str(option.mgoption.smoothingstep),'.pdf');
exportgraphics(gca,printstr)
hold off

%% L2 Error Computation

% Note about error types.
fprintf(['Note: All the errors listed above are residual ',...
  'errors from the iterative solver.\n']);

L2_error = getL2error(node,elem,pde.exactu,soln.u);
fprintf('\nFinal error in L2 sense = %g\n',L2_error);