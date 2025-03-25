function [u,flag,relres,iter,resvec] = myPCG(A,b,tol,maxit,M,u0)

% Initial residual.
u = u0;
r = b-A*u;

% Apply preconditioner.
z = M(r);
p = z;

% Compute initial residual norm.
resnorm = norm(r);
normb = norm(b);
if normb == 0
  normb = 1;
end

resvec = zeros(maxit+1,1);
resvec(1) = resnorm;

% Check initial convergence.
if resnorm/normb<tol
  flag = 0;
  relres = resnorm/normb;
  iter = 0;
  resvec = resvec(1);
  return;
end

% Main PCG loop.
flag = 1;
for iter = 1:maxit
  Ap = A*p;

  % Compute step size alpha.
  rz = r'*z;
  alpha = rz/(p'*Ap);

  % Update solution.
  u = u+alpha*p;

  % Update residual.
  r_new = r-alpha*Ap;
  resnorm = norm(r_new);
  resvec(iter+1) = resnorm;

  % Check convergence.
  if resnorm/normb<tol
    flag = 0;
    break;
  end

  % Apply preconditioner.
  z_new = M(r_new);

  % Compute beta.
  beta = (r_new'*z_new)/rz;

  % Update search direction.
  p = z_new+beta*p;

  % Prepare for next iteration.
  r = r_new;
  z = z_new;
end

relres = resnorm/normb;
resvec = resvec(1:iter+1);

if flag~=0
  warning('myPCG did not converge.');
end

end