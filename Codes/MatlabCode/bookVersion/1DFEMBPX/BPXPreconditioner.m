function z = BPXPreconditioner(r,MGParam)

% Getting number of levels.
L = MGParam.L;

% Phase 1: Restriction.
residuals = cell(L+1,1);
residuals{L+1} = r;
for l = L:-1:1
  residuals{l} = restriction(residuals{l+1});
end

% Phase 2: Scaling.
scaled = cell(L+1,1);
for l = 1:(L+1)
  nl = length(residuals{l});
  hl = 1/(nl+1);
  scaled{l} = hl*residuals{l};
end

% Phase 3: Accumulation using Horner-like scheme.
zTemp = scaled{1};
for l = 2:(L+1)
  zTemp = scaled{l}+prolongation(zTemp);
end
z = zTemp;
end