function [soln,eqn,info] = Elliptic(node,elem,bdFlag,pde,option)
%% ELLIPTIC Solve -div(d*grad(u)) + u = f: P1 linear element.
%
%   u = Elliptic(node,elem,bdFlag,pde) produces the linear finite element
%   approximation of the elliptic equation
%
%       - div(d * grad(u)) + u = f   in \Omega,
%
%   with boundary conditions:
%       Dirichlet boundary:    u = g_D on \Gamma_D,
%       Neumann boundary:     d*grad(u)*n = g_N on \Gamma_N,
%       Robin boundary:   g_R*u + d*grad(u)*n = g_N on \Gamma_R.
%
%   The mesh is given by node and elem, and bdFlag marks boundary edges
%   (see meshdoc, bddoc in iFEM). The data structure pde can contain
%   function handles for f, g_D, g_N, g_R, and d. The optional structure
%   'option' determines the solver, multigrid settings, etc.
%
%   Syntax:
%       [soln,eqn,info] = Elliptic(node,elem,bdFlag,pde,option)
%
%   The output structure soln contains:
%       - soln.u:  the nodal solution
%       - soln.Du: the piecewise gradients of the solution
%
%   The output structure eqn contains:
%       - eqn.A:        the modified stiffness matrix, including boundary
%                       conditions and the mass matrix
%       - eqn.b:        the modified load vector
%       - eqn.Lap:      the basic stiffness matrix from the diffusion term
%       - eqn.freeNode: the list of free nodes (not Dirichlet-fixed)
%
%   The output structure info contains:
%       - info.assembleTime: time (sec) to assemble matrices and vector
%       - info.solverTime:   time (sec) to solve the system
%       - info.itStep:       iteration count of mg/amg solver
%       - info.err:          residual norm
%       - info.flag:         solver status
%           0 = converged in fewer than maxIt steps
%           1 = did not converge by maxIt
%           2 = direct solver used
%           3 = no solver called
%
%   Boundary data pde.g_D, pde.g_N, pde.g_R can be scalars/vectors or
%   function handles. Similarly for pde.d (the diffusion coefficient)
%   and pde.f (the right-hand side).
%
%   By default:
%       - Direct solver is used for small problems (Ndof <= 2e3).
%       - Multigrid solver is used for larger ones (Ndof > 2e3).
%
%   See also Poisson, mg, amg, gradbasis, findboundary, quadpts.

%% -------------------- PART 0: SETUP --------------------
if ~exist('bdFlag','var'), bdFlag = []; end
if ~exist('option','var'), option = []; end

% Basic sizes
N  = size(node,1);   % number of vertices
NT = size(elem,1);   % number of triangles
Ndof = N;            % P1 element: Ndof = number of vertices

time = cputime;      % record assembling start time

%% -------------------- PART 1: DIFFUSION COEFFICIENT --------------------
if ~isfield(pde,'d'), pde.d = []; end
if ~isfield(option,'dquadorder'), option.dquadorder = 1; end

% K will be the elementwise average of the diffusion coefficient d
% if d is not empty.
K = []; 
if ~isempty(pde.d)
    if isnumeric(pde.d)
        % pde.d is a scalar or vector matching the size of elem
        K = pde.d;
    else
        % pde.d is a function handle
        [lambda,weight] = quadpts(option.dquadorder);
        nQuad = size(lambda,1);
        K = zeros(NT,1);
        for p = 1:nQuad
            pxy = lambda(p,1)*node(elem(:,1),:) ...
                + lambda(p,2)*node(elem(:,2),:) ...
                + lambda(p,3)*node(elem(:,3),:);
            K = K + weight(p)*pde.d(pxy);
        end
    end
end

%% -------------------- PART 2: GRADIENT BASIS --------------------
% Compute geometric quantities and gradient of local basis
[Dphi,area] = gradbasis(node,elem);
% Dphi(:,1,i) and Dphi(:,2,i) are partial derivatives of the i-th local basis
% function on each element.

%% -------------------- PART 3: ASSEMBLE STIFFNESS (DIFFUSION) MATRIX A ---
A = sparse(Ndof,Ndof);
for i = 1:3
    for j = i:3
        % Stiffness integrand: (d * grad phi_i) . (grad phi_j)
        Aij = (Dphi(:,1,i).*Dphi(:,1,j) + Dphi(:,2,i).*Dphi(:,2,j)).*area;
        if ~isempty(K), Aij = K.*Aij; end
        if (i == j)
            A = A + sparse(elem(:,i),elem(:,j),Aij,Ndof,Ndof);
        else
            A = A + sparse([elem(:,i); elem(:,j)], ...
                           [elem(:,j); elem(:,i)], ...
                           [Aij; Aij], Ndof, Ndof);
        end
    end
end

%% -------------------- PART 4: ASSEMBLE MASS MATRIX M FOR +u -------------
% The PDE has "+u" term, so we need:
%    M_{ij} = \int_{\Omega} phi_i * phi_j dx.
M = sparse(Ndof,Ndof);

% By default, we can integrate mass terms at order 2 or 3 to be safe
if ~isfield(option,'mquadorder'), option.mquadorder = 2; end
[lambda_m,weight_m] = quadpts(option.mquadorder);
nQuadM = size(lambda_m,1);

for i = 1:3
    for j = i:3
        % M_{ij} on each element
        Mij = zeros(NT,1);
        for p = 1:nQuadM
            pxy = lambda_m(p,1)*node(elem(:,1),:) ...
                + lambda_m(p,2)*node(elem(:,2),:) ...
                + lambda_m(p,3)*node(elem(:,3),:);
            % phi_i and phi_j are just barycentric coords for P1
            phi_i_val = lambda_m(p,i);
            phi_j_val = lambda_m(p,j);
            Mij = Mij + weight_m(p) * phi_i_val * phi_j_val;
        end
        Mij = Mij .* area;  % multiply by element area
        if (i == j)
            M = M + sparse(elem(:,i),elem(:,j),Mij,Ndof,Ndof);
        else
            M = M + sparse([elem(:,i); elem(:,j)], ...
                           [elem(:,j); elem(:,i)], ...
                           [Mij; Mij], Ndof, Ndof);
        end
    end
end

% Combine diffusion + mass
A0 = A;    % keep a copy of the original diffusion matrix
A  = A + M;

%% -------------------- PART 5: ASSEMBLE LOAD VECTOR b --------------------
b = zeros(Ndof,1);
if ~isfield(option,'fquadorder'), option.fquadorder = 3; end

% If pde.f is absent or zero, we skip. Otherwise, we handle numeric or function forms:
if ~isfield(pde,'f'), pde.f = []; end
if isreal(pde.f) && all(pde.f == 0)
    pde.f = [];
end

if isreal(pde.f) && ~isempty(pde.f)
    % pde.f might be a scalar, a length-NT vector, or length-N vector
    if length(pde.f) == NT
        % piecewise constant on each element
        bt = pde.f .* area / 3; 
        b  = accumarray(elem(:), [bt; bt; bt], [Ndof,1]);
    elseif length(pde.f) == N
        % piecewise linear if we have nodal data
        bt = zeros(NT,3);
        bt(:,1) = area .* (2*pde.f(elem(:,1)) + pde.f(elem(:,2)) + pde.f(elem(:,3)))/12;
        bt(:,2) = area .* (2*pde.f(elem(:,2)) + pde.f(elem(:,3)) + pde.f(elem(:,1)))/12;
        bt(:,3) = area .* (2*pde.f(elem(:,3)) + pde.f(elem(:,1)) + pde.f(elem(:,2)))/12;
        b       = accumarray(elem(:), bt(:), [Ndof,1]);
    elseif length(pde.f) == 1
        % constant scalar f
        bt = pde.f * area/3;
        b  = accumarray(elem(:), [bt; bt; bt], [Ndof,1]);
    else
        error('Dimension mismatch: pde.f is numeric but not of length N or NT.');
    end
elseif isa(pde.f,'function_handle')
    % f is a function, integrate with chosen quadrature
    [lambda_f,weight_f] = quadpts(option.fquadorder);
    nQuadF = size(lambda_f,1);
    bt = zeros(NT,3);
    for p = 1:nQuadF
        pxy = lambda_f(p,1)*node(elem(:,1),:) ...
            + lambda_f(p,2)*node(elem(:,2),:) ...
            + lambda_f(p,3)*node(elem(:,3),:);
        fp = pde.f(pxy);
        for i = 1:3
            bt(:,i) = bt(:,i) + weight_f(p)*lambda_f(p,i)*fp;
        end
    end
    bt = bt .* repmat(area,1,3);
    b  = accumarray(elem(:), bt(:), [Ndof,1]);
end

%% -------------------- PART 6: APPLY BOUNDARY CONDITIONS -----------------
[AD,b,u,freeNode,isPureNeumann] = getbd(A,b);

%% Record assembly time
assembleTime = cputime - time;
if ~isfield(option,'printlevel'), option.printlevel = 1; end
if option.printlevel >= 2
    fprintf('Time to assemble matrix equation = %4.2g s\n',assembleTime);
end

%% -------------------- PART 7: SOLVE THE LINEAR SYSTEM ------------------
if isempty(freeNode)
    % If no free nodes, problem is trivial (all Dirichlet?)
    soln = u; eqn = []; info = [];
    return;
end

% Decide solver if not explicitly given
Ndof = length(u);  %#ok<NASGU>
if ~isfield(option,'solver') || isempty(option.solver)
    if length(u) <= 2e3
        option.solver = 'direct';
    else
        option.solver = 'mg';
    end
end

solver = option.solver;
timeSolve = cputime;

switch solver
    case 'direct'
        u(freeNode) = AD(freeNode,freeNode)\b(freeNode);
        residual = norm(b - AD*u);
        info = struct('solverTime',cputime - timeSolve,...
                      'itStep',0,'err',residual,'flag',2,'stopErr',residual);

    case 'none'
        % Just return assembled system, no solution
        info = struct('solverTime',[],'itStep',0,'err',[],'flag',3,'stopErr',[]);

    case 'mg'
        if ~isfield(option,'mgoption'), option.mgoption = []; end
        if ~isfield(option.mgoption,'x0')
            option.mgoption.x0 = u;
        end
        if ~isfield(option.mgoption,'solver')
            option.mgoption.solver = 'CG';
        end
        [u,info] = mg(AD,b,elem,option.mgoption);

    case 'amg'
        if ~isfield(option,'amgoption'), option.amgoption = []; end
        if ~isfield(option.amgoption,'x0')
            option.amgoption.x0 = u;
        end
        if ~isfield(option.amgoption,'solver')
            option.amgoption.solver = 'CG';
        end
        [u(freeNode),info] = amg(AD(freeNode,freeNode), b(freeNode), option.amgoption);

    otherwise
        error('Unknown solver type.');
end

% Handle pure Neumann (if relevant)
if isPureNeumann
    % Force mean(u) = 0
    patchArea = accumarray(elem(:), [area; area; area]/3, [N,1]);
    uc = sum(u .* patchArea) / sum(area);
    u  = u - uc;
end

%% -------------------- PART 8: GRADIENT POST-PROCESS ---------------------
% Compute the gradient Du on each element
dudx =  u(elem(:,1)).*Dphi(:,1,1) + ...
        u(elem(:,2)).*Dphi(:,1,2) + ...
        u(elem(:,3)).*Dphi(:,1,3);
dudy =  u(elem(:,1)).*Dphi(:,2,1) + ...
        u(elem(:,2)).*Dphi(:,2,2) + ...
        u(elem(:,3)).*Dphi(:,2,3);
Du = [dudx, dudy];

%% -------------------- PART 9: OUTPUT STRUCTURES -------------------------
if nargout == 1
    soln = u;
else
    soln = struct('u',u,'Du',Du);
    eqn  = struct('A',AD,'b',b,'freeNode',freeNode,'Lap',A0);
    info.assembleTime = assembleTime;
    info.solverTime   = info.solverTime  ;  % from the solve
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      S U B F U N C T I O N  getbd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [AD,b,u,freeNode,isPureNeumann] = getbd(A,b)
    %% Enforce Dirichlet/Neumann/Robin boundary conditions.
    %
    %  This function modifies A and b to incorporate boundary data:
    %    - Dirichlet:   u = g_D on \Gamma_D
    %    - Neumann:     d grad(u)*n = g_N on \Gamma_N
    %    - Robin:       g_R*u + d grad(u)*n = g_N on \Gamma_R
    %
    %  In the pure Neumann case, the matrix is singular. We fix the mean zero
    %  condition by adjusting b and adding a tiny diagonal to the matrix at
    %  node 1 (or simply leaving it singular and removing the rank-1 null).

    u = zeros(Ndof,1);

    if ~isfield(pde,'g_D'), pde.g_D = []; end
    if ~isfield(pde,'g_N'), pde.g_N = []; end
    if ~isfield(pde,'g_R'), pde.g_R = []; end

    % Identify boundary edges and boundary nodes
    if ~isempty(bdFlag)
        [fixedNode,bdEdge,isBdNode] = findboundary(elem,bdFlag);
        freeNode = find(~isBdNode);
    else
        fixedNode = [];
        freeNode  = (1:Ndof)';
    end

    % Handle Robin boundary (bdFlag==3)
    isRobin  = (bdFlag(:) == 3);
    allEdges = [elem(:,[2,3]);
                elem(:,[3,1]);
                elem(:,[1,2])];
    Robin = allEdges(isRobin,:);
    if ~isempty(Robin) && ~isempty(pde.g_R) && ~(isnumeric(pde.g_R)&&(all(pde.g_R==0)))
        ve = node(Robin(:,1),:) - node(Robin(:,2),:);
        edgeLength = sqrt(sum(ve.^2,2));
        mid        = 0.5*(node(Robin(:,1),:) + node(Robin(:,2),:));
        % Add \int g_R phi_i phi_j ds
        ii = [Robin(:,1); Robin(:,1); Robin(:,2); Robin(:,2)];
        jj = [Robin(:,1); Robin(:,2); Robin(:,1); Robin(:,2)];
        temp  = pde.g_R(mid).*edgeLength;
        ss    = [temp/3; temp/6; temp/6; temp/3];
        A     = A + sparse(ii,jj,ss,Ndof,Ndof);
    end

    % Determine if pure Neumann:
    % It happens if there are no fixedNodes and no Robin conditions that fix
    % the solution, i.e. the matrix is singular. We fix the rank deficiency.
    isPureNeumann = false;
    if isempty(fixedNode) && isempty(Robin)
        isPureNeumann = true;
        % Keep A but typically fix dof #1 or add a small shift
        A(1,1) = A(1,1) + 1e-8;
    end

    % Neumann boundary
    Neumann = [];
    if ~isempty(bdFlag)
        % bdFlag marks edges with 1,2,3,...; we assume 2 => Neumann, 3 => Robin
        Neumann = bdEdge(bdFlag(bdEdge(:,end))==2,1:2);  %#ok<*NASGU>
    else
        % If no bdFlag is present but pde.g_N or pde.g_R is non-empty,
        % treat the entire boundary as Neumann or Robin
        if ~isempty(pde.g_N) || ~isempty(pde.g_R)
            [~,Neumann] = findboundary(elem);
        end
    end
    % If pde.g_N is not empty or zero, add the boundary integral
    if ~isempty(Neumann) && ~isempty(pde.g_N)
        if ~isfield(option,'gNquadorder'), option.gNquadorder = 2; end
        [lambdaN,weightN] = quadpts1(option.gNquadorder);
        nQuadN = length(weightN);
        el     = sqrt(sum((node(Neumann(:,1),:) - node(Neumann(:,2),:)).^2,2));
        ge     = zeros(size(Neumann,1),2);
        for p = 1:nQuadN
            pxy  = lambdaN(p,1)*node(Neumann(:,1),:) ...
                 + lambdaN(p,2)*node(Neumann(:,2),:);
            gNp  = pde.g_N(pxy);
            for ig = 1:2
                ge(:,ig) = ge(:,ig) + weightN(p)*lambdaN(p,ig)*gNp;
            end
        end
        ge = ge.*repmat(el,1,2);
        b  = b + accumarray(Neumann(:), ge(:), [Ndof,1]);
    end

    % Dirichlet boundary
    if ~isempty(fixedNode)
        % Possibly non-homogeneous g_D
        if ~isempty(pde.g_D) && ~(isnumeric(pde.g_D) && all(pde.g_D==0))
            if isnumeric(pde.g_D)
                u(fixedNode) = pde.g_D(fixedNode);
            else
                u(fixedNode) = pde.g_D(node(fixedNode,:));
            end
        end
        % Impose in matrix
        bdidx        = zeros(Ndof,1);
        bdidx(fixedNode) = 1;
        Tbd          = spdiags(bdidx,0,Ndof,Ndof);
        T            = spdiags(1-bdidx,0,Ndof,Ndof);
        AD           = T*A*T + Tbd;
        b            = b - A*u;        % subtract A*u from RHS
        b(fixedNode) = u(fixedNode);   % then set these values
    else
        % pure Neumann or Robin
        AD = A;
    end

    if isPureNeumann
        % Enforce compatibility:
        % The integral of b over domain must be zero for existence of solution
        b = b - mean(b);
    end
    end % getbd
end
