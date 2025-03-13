%% Grid with Boundary Function
function f = grid_with_boundary(k,H1,H2)
    n = 2^k;
    
    % Grid:
    f = zeros(n+1, n+1);
    
    % Boundary Conditions:
    % Bottom
    f(1, 1:(n/2 +1)) = linspace(13, 5, n/2 +1);
    f(1, (n/2+1):(3*n/4)) = 5;
    f(1, (3*n/4+1):(n+1)) = linspace(5, 13, n/4+1);
    
    % Top
    f(n+1, 1:n+1) = 21;
    
    % Left
    f(1:(3*n/8+1), 1) = linspace(13, H1, 3*n/8+1);
    f((n/2+1):n+1, 1) = linspace(H1, 21, n/2+1);
    
    % Right
    f(1:(n/2+1), n+1) = linspace(13, H2, n/2+1);
    f((5*n/8+1):n+1, n+1) = linspace(H2, 21, 3*n/8+1);

    % Heaters:
    f = add_heaters(f,H1,H2);
end