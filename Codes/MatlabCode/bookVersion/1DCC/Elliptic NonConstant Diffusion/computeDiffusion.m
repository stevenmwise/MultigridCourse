function D = computeDiffusion(nf,hf,xLower)
D = zeros(nf+1,1);

for i = 1:nf+1
  xEdge = (i-1)*hf+xLower;
  
  % Example variable coefficient.
  D(i) = 1.0+0.3*sin(2*pi*xEdge);
end

end