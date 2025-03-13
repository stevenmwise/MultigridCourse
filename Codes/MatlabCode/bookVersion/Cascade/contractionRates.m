function [rate] = contractionRates(errVals,kStop)

kStop(:) = kStop(:)-mod(kStop(:),4);
rate = zeros(3,1);
%
% Estimate the rates of contraction:
if min(kStop) >= 4 
  for i = 1:3
    kv = [kStop(i)/2,kStop(i)/2+kStop(i)/4,kStop(i)];
    le = log(errVals(kv,i));
    p1 = polyfit(kv,le,1);
    rate(i) = exp(p1(1));
  end
end

end