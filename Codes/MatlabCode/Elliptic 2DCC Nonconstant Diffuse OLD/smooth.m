function u = smooth(u,f,hf,MGParam)

nfPlusGhostLayers = size(u);
nf = nfPlusGhostLayers-2;
hf2 = hf*hf;

m = MGParam.m1;
xLower = MGParam.xLower;

dEW = zeros(nf(1)+1,nf(2));
dNS = zeros(nf(1),nf(2)+1);
cCC = zeros(nf(1),nf(2));
%
% EW edges
for i = 1:nf(1)+1
  x = (i-1)*hf+xLower(1);
  for j = 1:nf(2)
    y = (j-0.5)*hf+xLower(2);
    dEW(i,j) = 2.0*x*x+y*y+1.0;
  end
end
%
% NS edges
for i = 1:nf(1)
  x = (i-0.5)*hf+xLower(1);
  for j = 1:nf(2)+1
    y = (j-1)*hf+xLower(2);
    dNS(i,j) = 2.0*x*x+y*y+1.0;
  end
end
%
% Cell centers
for i = 1:nf(1)
  x = (i-0.5)*hf+xLower(1);
  for j = 1:nf(2)
    y = (j-0.5)*hf+xLower(2);
    cCC(i,j) = exp(x)+exp(y);
  end
end

u = applyBCs(u);

for k = 1:m
  for j = 2:nf(2)+1
    for i = 2:nf(1)+1
      u(i,j) = (hf2*f(i-1,j-1) ...
        + dEW(i,j-1)*u(i+1,j)+dEW(i-1,j-1)*u(i-1,j) ...
        + dNS(i-1,j)*u(i,j+1)+dNS(i-1,j-1)*u(i,j-1)) ...
        / (dEW(i,j-1)+dEW(i-1,j-1)+dNS(i-1,j)+dNS(i-1,j-1) ...
        + hf2*cCC(i-1,j-1));  
    end
  end

  u = applyBCs(u);

end

end