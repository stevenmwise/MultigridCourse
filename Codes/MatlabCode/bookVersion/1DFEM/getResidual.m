function res = getResidual(f,u)

nl = length(f);
hl = 1.0/(nl+1);
res = zeros(nl,1);

% Compute residual at boundaries.
res(1) = f(1)-(2.0*u(1)-u(2))/hl;
res(nl) = f(nl)-(-u(nl-1)+2.0*u(nl))/hl;

% Compute residual at interior points.
res(2:nl-1) = f(2:nl-1)-(-u(1:nl-2)+2*u(2:nl-1)-u(3:nl))/hl;

end