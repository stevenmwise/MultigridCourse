function [vProl] = prolongationCC(v)

nc = length(v);
nf = 2*nc;
vProl = zeros(nf,1);

vProl(1:2:nf-1) = v(1:nc);
vProl(2:2:nf  ) = v(1:nc);

end