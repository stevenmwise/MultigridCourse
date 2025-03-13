function cD = restrictDiff(D)

% Restrict the diffusion coefficient by injection.
cD = D(1:2:end);

end
