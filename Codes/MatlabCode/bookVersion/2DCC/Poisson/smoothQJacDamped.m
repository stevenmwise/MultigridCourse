function u = smoothQJacDamped(f,u,hf,m,omega)

[nFull,~] = size(u);
n = nFull-2;
hf2 = hf*hf;
omegaPrime = 1.0-omega;

z = zeros(n,n);

u = applyBCs(u);
for sweep = 1:m
  for i = 2:n+1
    for j = 2:n+1
        
       z(i-1,j-1) = (hf2*f(i-1,j-1)+u(i+1,j)+u(i-1,j)...
           +u(i,j+1)+u(i,j-1))/4.0;
    end
  end
  u(2:n+1,2:n+1) = omega*z(1:n,1:n)+omegaPrime*u(2:n+1,2:n+1);
  u = applyBCs(u);
end

end