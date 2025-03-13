%% WARNING! THIS CODE NEEDS TO BE BEAUTIFIED!!

theta = 2.000116726443145e-01;
rot = [cos(theta),sin(theta);-sin(theta),cos(theta)];
A = [1,0;0,10];
A = rot'*A*rot;
f = [1;1];
uexact = A\f

phi = @(x,y,A,f)  0.5*A(1,1)*x.^2 +0.5*A(2,2)*y.^2 + A(1,2)*x.*y - f(1)*x - f(2)*y;

[X,Y] = meshgrid(-3:0.01:3);
EE = phi(X,Y,A,f);
contour(X,Y,EE);
grid on
hold on

fg = figure(1);

plot(uexact(1),uexact(2),'ro')
hold on

u = [2;1];
for i=1:2000
   r = f - A*u;
   p = prec(r);
   Ap = A*p;
   t = (r'*p)/(p'*Ap);
   unew = u + t*p;
   plot( [u(1),unew(1)],[u(2),unew(2)], 'rx-' )
   hold on

   nerror = norm( unew - uexact );
   fprintf(' %d %g \n ', i, nerror );

   if nerror < 1e-6
       fprintf('Solution found after %d iterations! \n ', i)
       break
   end
   u = unew;
end

%% Octave output - make it matlab
% print( fg, 'SteepestDescent.pdf', '-dpdflatex' );

hold off

function p = prec( r )
    p = r;
end
