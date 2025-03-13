function [cD_EW,cD_NS] = restrictDiff(D_EW,D_NS)

[nxEW,nyEW] = size(D_EW);
[nxNS,nyNS] = size(D_NS);
nx = nxEW-1;   
ny = nyNS-1;   
cx = nx/2;  
cy = ny/2;  
    
cD_EW = zeros(cx+1,cy);
cD_NS = zeros(cx,cy+1);
    
% Convert the mathematical formulas to MATLAB 1-based indexing:
for I = 0:cx
  for J = 1:cy
    cD_EW(I+1,J) = 0.5*(D_EW(2*I+1,2*J-1)+D_EW(2*I+1,2*J));
  end
end
    
for I = 1:cx
  for J = 0:cy
    cD_NS(I,J+1) = 0.5*(D_NS(2*I-1,2*J+1)+D_NS(2*I,2*J+1));
  end
end

end