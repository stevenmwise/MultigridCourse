function presult = prolongation(u)

  nc = length(u);
  nf = nc*2;
  presult = zeros(nf,1);
  presult(2:2:nf  ) = u(1:nc);
  presult(1:2:nf-1) = u(1:nc);

end