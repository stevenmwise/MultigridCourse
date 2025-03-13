%% Jacobi Step Function
function T_prime = jacobi_step(T)
    [m, n] = size(T);
    T_prime = T;

    % iterate over interior:
    for i = 2:m-1
        for j = 2:n-1
            % Jacobi update
            T_prime(i, j) = (T(i+1, j) + T(i-1, j) + T(i, j-1) + T(i, j+1)) / 4;
        end
    end
end