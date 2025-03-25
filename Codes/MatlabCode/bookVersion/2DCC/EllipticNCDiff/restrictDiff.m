function [cDeW,cDnS] = restrictDiff(DeW,DnS)

[nxEW,ny] = size(DeW);
[nx,nyNS] = size(DnS);

ncx = nx/2;  
ncy = ny/2;

cDeW(1:ncx+1,1:ncy) ...
  = 0.5*(DeW(1:2:nx+1,1:2:ny-1)+DeW(1:2:nx+1,2:2:ny));

cDnS(1:ncx,1:ncy+1) ...
  = 0.5*(DnS(1:2:nx-1,1:2:ny+1)+DnS(2:2:nx,1:2:ny+1));

end