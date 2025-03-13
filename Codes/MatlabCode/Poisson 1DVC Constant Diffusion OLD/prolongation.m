function [vProl] = prolongation(v)

nc = length(v);
nf = 2*(nc+1)-1;
vProl = zeros(nf,1);

vProl( 1) = v( 1)/2.0;
vProl(nf) = v(nc)/2.0;

vProl(2:2:2*nc  ) =  v(1:nc);
vProl(3:2:2*nc-1) = (v(1:nc-1)+v(2:nc))/2.0;

end