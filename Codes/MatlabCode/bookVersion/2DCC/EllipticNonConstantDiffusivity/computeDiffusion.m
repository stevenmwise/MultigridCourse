function [D_EW,D_NS] = computeDiffusion(nf,hf,xLower)

nx = nf(1);
ny = nf(2);

D_EW = zeros(nx+1,ny);
D_NS = zeros(nx,ny+1);

for iE = 1:(nx+1)
  for j = 1:ny

    % This edge is conceptually between cell iE-1 and ...
    % cell iE (for x). % The center of that edge is:
    xEdge = (iE-0.5)*hf+xLower(1);
    yEdge = (j -0.5)*hf+xLower(2);

    % Evaluate D at that midpoint.
    D_EW(iE,j) = 1.0+0.7*sin(pi*xEdge)*cos(pi*yEdge);
  end
end

for i = 1:nx
  for jE = 1:(ny+1)

    % This edge is conceptually between cell jE-1 and ...
    % cell jE (for y). The center of that edge is:
    xEdge = (i  -0.5)*hf + xLower(1);
    yEdge = (jE -0.5)*hf + xLower(2);
    
    % Evaluate D at that midpoint.
    D_NS(i,jE) = 1.0+0.7*sin(pi*xEdge)*cos(pi*yEdge);
  end
end

end
