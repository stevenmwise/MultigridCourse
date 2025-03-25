function operatorResult = FDOperator(u,hf,MGParam)
  
nfPlusGhostLayers = size(u);
nf = nfPlusGhostLayers-2;
hf2 = hf*hf;

operatorResult = (-     u(3:nf(1)+2,2:nf(2)+1)...
                  -     u(1:nf(1)  ,2:nf(2)+1)...
                  -     u(2:nf(1)+1,3:nf(2)+2)...
                  -     u(2:nf(1)+1,1:nf(2)  )...
                  + 4.0*u(2:nf(1)+1,2:nf(2)+1))/hf2;


end