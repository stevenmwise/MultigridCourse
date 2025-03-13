function restrictionResult = restriction(u)

nf = length(u);
nc = nf/2;
restrictionResult = 0.5*(u(2:2:nf)+u(1:2:nf-1));

end