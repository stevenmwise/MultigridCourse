function presult = prolongation(u)

nc = size(u);
nf = nc*2;
presult = zeros(nf);

presult(2:2:nf(1)  ,2:2:nf(2)  ) = u(1:nc(1),1:nc(2));
presult(1:2:nf(1)-1,2:2:nf(2)  ) = u(1:nc(1),1:nc(2));
presult(2:2:nf(1)  ,1:2:nf(2)-1) = u(1:nc(1),1:nc(2));
presult(1:2:nf(1)-1,1:2:nf(2)-1) = u(1:nc(1),1:nc(2));

end
