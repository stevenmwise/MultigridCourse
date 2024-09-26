function [vRest] = restrictionCC(v)

nf = length(v);
nc = nf/2;
vRest = zeros(nc,1);

vRest(1:nc) = v(1:2:nf-1)+v(2:2:nf);

end