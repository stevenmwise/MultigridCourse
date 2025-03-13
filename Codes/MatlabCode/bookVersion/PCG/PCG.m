%
% THIS NEEDS TO BE CHECKED FOR CORRECTNESS AND BEAUTIFIED
%
function [x, iStop, err] = PCG( A, x0, f, xexact, maxit, tol, Prec )
% The preconditioned conjugate gradient method to approximate 
% the solution to
  err = zeros(maxit,1);
  x = x0;
  r = f - A*x;
  z = Prec(r);
  p = z;
  for its = 1:maxit
     Ap = A*p;
     denom = p'*Ap;
     num = z'*p;
     alpha = num/denom;
     x = x + alpha*p;
     r = r - alpha*Ap;
     z = Prec(r);
     denom = num;
     num = z'*p;
     beta = num/denom;
     p = z + beta*p;
     eee = norm( xexact - x );
     err(its) = eee;
     %fprintf( 'iteration = %d \t \t residual = %g \n', its, eee );
     if eee < tol
       iStop = its;
       return;
     end
  end
  iStop = maxit;
end
