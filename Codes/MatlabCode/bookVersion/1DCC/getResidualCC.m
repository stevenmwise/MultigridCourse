function [res] = getResidualCC(f,u)

nl = length(f);
hl = 1.0/nl;
res = zeros(nl,1);

u(   1) = -u(   2);
u(nl+2) = -u(nl+1);

res(1:nl) = hl*hl*f(1:nl)-(-u(1:nl)+2*u(2:nl+1)-u(3:nl+2));

end