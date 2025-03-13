function operatorResult = FDOperator(u,hf,MGParam)
  
nfPlusGhostLayers = size(u);
nf = nfPlusGhostLayers-2;
hf2 = hf*hf;
%
% Vectorized...
%
% -dicreteLaplacian(u) + u: 
%
%  operatorResult = (-     u(3:nf(1)+2,2:nf(2)+1) ...
%                    -     u(1:nf(1)  ,2:nf(2)+1) ...
%                    -     u(2:nf(1)+1,3:nf(2)+2) ...
%                    -     u(2:nf(1)+1,1:nf(2)  ) ...
%                    + 4.0*u(2:nf(1)+1,2:nf(2)+1))/hf2 ...
%                    +     u(2:nf(1)+1,2:nf(2)+1);
%
% Non-Vectorized...
%
operatorResult = zeros(nf(1),nf(2));
% 
xLower = MGParam.xLower;
%
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
%
for i = 2:nf(1)+1
  for j = 2:nf(2)+1
    operatorResult(i-1,j-1) ...
      = -( dEW(i  ,j-1)*(u(i+1,j)-u(i,j)) ...
          -dEW(i-1,j-1)*(u(i,j)-u(i-1,j)))/hf2 ...
        -( dNS(i-1,  j)*(u(i,j+1)-u(i,j)) ...
          -dNS(i-1,j-1)*(u(i,j)-u(i,j-1)))/hf2 ...
        + u(i,j)*cCC(i-1,j-1);
  end
end

end