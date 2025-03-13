function presult = prolongation(u)
  nc = size(u);
  nf = nc*2;
  presult = zeros(nf(1), nf(2));
    
  presult(2:2:nf(1)  ,2:2:nf(2)  ) = u;
  presult(1:2:nf(1)-1,2:2:nf(2)  ) = u;
  presult(2:2:nf(1)  ,1:2:nf(2)-1) = u;
  presult(1:2:nf(1)-1,1:2:nf(2)-1) = u;
end