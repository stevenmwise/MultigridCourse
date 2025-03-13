clear

n = 63;
h = 1/(n+1);
h2 = h*h;
kMax = 30000;
tol = 1e-07;

errVals = zeros(kMax,2);
kStop = zeros(2,1);
%
x = linspace(h,1.0-h,n)';
A = gallery('tridiag',n,-1.0/h2,2.0/h2,-1.0/h2);
wExact = exp(sin(3.0*pi*x))-1.0;
f = A*wExact;
%
% Initialize w
SEED = 1234;
rng(SEED);
w = rand(n,1)-0.5;

% We run the PCG with identity preconditioner
[~, kStop(1,1), errVals(:,1) ] = PCG( A, w, f, wExact, kMax, tol, @PrecIdentity);
%
% We create the figure for the identity preconditioner
semilogy( 1:kStop(1,1), errVals(1:kStop(1,1),1), 'r-', 'LineWidth', 1.5 )
grid on
hold on
%
% We run the PCG with Jacobi preconditioner
%
PrecRichardson = @(z) PrecRichardsonDriver( 4.0*h, z );
[~, kStop(2,1), errVals(:,2) ] = PCG( A, w, f, wExact, kMax, tol, PrecRichardson );
%
% We create the figure for the Jacobi preconditioner
semilogy( 1:kStop(2,1), errVals(1:kStop(2,1),2), 'b-', 'LineWidth', 1.5 )
grid on
hold on
% Make the figure nice
xlabel('$k$', 'Interpreter','latex');
title('PCG Iteration Errors', 'Interpreter','latex');
legend( 'Identity Preconditioner', 'Richradson Preconditioner', 'Interpreter', 'latex');

printstr = 'PCGErrors.pdf';
exportgraphics( gca, printstr );

hold off


function p = PrecIdentity( z )
    p = z;
end

function p = PrecRichardsonDriver( h, z )
    p = h*z;
end
