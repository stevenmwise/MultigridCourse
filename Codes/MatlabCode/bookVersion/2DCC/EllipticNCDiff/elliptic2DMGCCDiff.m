%
%   elliptic2DMGCCDiff.m
%
% Calling script for solving the 2D elliptic equation
%
%   -div(D grad u) + rC*u = f
%
% on a rectangular domain, where D is a non-constant
% diffusion coefficient, and rC is a constant. The boundary 
% conditions may be taken to be (i) homogeneous Dirichlet,
% (ii) homogeneous Neumann, or (iii) periodic.
% 
% Discretization is accomplished using the cell-centered 
% finite difference method. The cell-centered multigrid method
% is used as the solver.
%
%% Initialization and Parameter

% Clear command line window and workspace.
clear; clc;

% Tolerance of residual to stop multigrid iterations.
tol = 1.0e-10;

% Number of interior cells in x,y at finest level.
nL = [512,512]; 

% Number of multigrid levels.
L = 9;            

% Define the domain of the equation.
xLower = [0.0,0.0];
xUpper = [3.2,3.2];

hLvec = (xUpper-xLower)./nL;

% Parameters for the multigrid solver.
pCycle = 1;
m1     = 3;
m2     = 3;
omega  = 0.8;
kMax   = 100;
rC     = 2.0;

%% Setting up the multigrid solver

% Check mesh spacing.
if abs(hLvec(1)-hLvec(2)) > 1.0e-14
  error('Grid spacing mismatch in x vs y');
end
hL = hLvec(1);

% Check that nL is divisible by 2^(L-1).
nl = nL;
for k = L-1:-1:0
  if any(mod(nl,2) ~= 0)
    error('nL not divisible by 2^L.  Coarsening error.');
  end
  nl = nl/2;
end

MGParam.nL     = nL;
MGParam.xLower = xLower;
MGParam.xUpper = xUpper;
MGParam.L      = L;
MGParam.pCycle = pCycle;
MGParam.m1     = m1;
MGParam.m2     = m2;
MGParam.omega  = omega;
MGParam.kMax   = kMax;
MGParam.tol    = tol;
MGParam.rC     = rC;


% Allocate solution arrays (with ghost layers).
u      = zeros(nL(1)+2,nL(2)+2);
uExact = zeros(nL(1)+2,nL(2)+2);
f      = zeros(nL(1)  ,nL(2)  );
DeW    = zeros(nL(1)+1,nL(2)  );
DnS    = zeros(nL(1)  ,nL(2)+1);

% Build an exact solution and define f accordingly, 
% so we can test the true algebraic error.
for i = 1:nL(1)
  x = (i-0.5)*hL+xLower(1);
  xArg = 2.0*pi*x/xUpper(1);
  for j = 1:nL(2)
    y = (j-0.5)*hL+xLower(2);
    yArg = 2.0*pi*y/xUpper(2);
    uExact(i+1,j+1) = exp(cos(xArg)*cos(yArg));
  end
end

% We next build the finest-grid diffusion arrays DeW, DnS:
for i = 1:nL(1)+1
  xEdge = (i-1.0)*hL+xLower(1);
  xArg = 2.0*pi*x/xUpper(1);
  for j = 1:nL(2)
    yEdge = (j-0.5)*hL+xLower(2);
    yArg = 2.0*pi*y/xUpper(2);

% Evaluate D at the east-west edge point:
    DeW(i,j) = 1.0+0.5*sin(xArg)*cos(yArg);
  end
end

for i = 1:nL(1)
  xEdge = (i-0.5)*hL+xLower(1);
  xArg = 2.0*pi*x/xUpper(1);
  for j = 1:nL(2)+1
    yEdge = (j-1.0)*hL+xLower(2);
    yArg = 2.0*pi*y/xUpper(2);
    
% Evaluate D at the north-south edge point:
    DnS(i,j) = 1.0+0.5*sin(xArg)*cos(yArg);
  end
end

% Apply the BCs to uExact, then compute f by calling 
% FDOperator.
uExact = applyBCs(uExact);
f = FDOperator(uExact,DeW,DnS,hL,MGParam);

% Initialize our approximate solution (e.g. random):
%rng(1234);
%u = rand(size(u))-0.5;

%% Call the Multigrid solver.

tic;
[u,errVals,kStop] = multiGridSolver(u,f,DeW,DnS,hL,MGParam, ...
  uExact);
toc;

%% Plotting the numerical solution and error plots.

uPlot(1:nL(1),1:nL(2)) = u(2:nL(1)+1,2:nL(2)+1);

% Plot solution details.
xCC = linspace(xLower(1)+hL/2,xUpper(1)-hL/2,nL(1));
yCC = linspace(xLower(2)+hL/2,xUpper(2)-hL/2,nL(2));
[x,y] = meshgrid(xCC,yCC);

figure(1)
clf
surf(x,y,uPlot,'EdgeColor','none');
title('Numerical Solution u');
xlabel('x');
ylabel('y');
zlabel('u(x,y)');
colorbar;
view(2);
axis equal tight;

figure(2)
clf
semilogy(errVals(1:kStop,3),'bo','LineWidth',1.5)
hold on
semilogy(errVals(1:kStop,2),'rs','LineWidth',1.5)
hold on
semilogy(errVals(1:kStop,1),'kd','LineWidth',1.5)
hold on

rate = 1.0;
if kStop >= 4
  kv = kStop-3:kStop;
  le = log(errVals(kStop-3:kStop,3));
  p1 = polyfit(kv,le,1);
  rate = exp(p1(1));
  p1k = polyval(p1,1:kStop);
  semilogy(1:kStop,exp(p1k),'-k')
else
  semilogy(errVals(1:kStop,1),'.','LineWidth',1.5)
end

xlabel('$k$','Interpreter','latex');
title('Multigrid Iteration Errors',...
  'Interpreter','latex');
legend('$\left\|{\bf u}_L^{\rm E}-{\bf u}_L^k\right\|_L$',...
  '$\left\|{\bf u}_L^k-{\bf u}_L^{k-1} \right\|_L$',...
  '$\left\|{\bf r}_L^k\right\|_L$',...
  'log-linear fit','Interpreter','latex');
text(1.5,128*MGParam.tol,strcat(...
  '$\gamma_{\rm comp} =\hspace{.1cm}$',...
  num2str(rate,'%10.5e')),'FontSize',14,...
  'Interpreter','latex')
text(1.5,32*MGParam.tol,strcat('$m_1 =\hspace{.1cm}$',...
  num2str(MGParam.m1)),'FontSize',14,...
  'Interpreter','latex')
text(1.5,8*MGParam.tol,strcat('$m_2 =\hspace{.1cm}$',...
  num2str(MGParam.m2)),'FontSize',14,...
  'Interpreter','latex')
text(1.5,2*MGParam.tol,strcat('$p =\hspace{.1cm}$',...
  num2str(MGParam.pCycle)),'FontSize',14,...
  'Interpreter','latex')
printstr = strcat('Err_nL_',num2str(MGParam.nL(1)),...
  '_m1_',num2str(MGParam.m1),'.pdf');
exportgraphics(gca, printstr)
hold off