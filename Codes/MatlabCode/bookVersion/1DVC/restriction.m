function [vRest] = restriction(v)

nf = length(v);
nc = (nf+1)/2-1;
vRest = zeros(nc,1);

vRest(1:nc) = 0.5*(v(1:2:2*nc-1)+v(3:2:2*nc+1))+v(2:2:2*nc);

end