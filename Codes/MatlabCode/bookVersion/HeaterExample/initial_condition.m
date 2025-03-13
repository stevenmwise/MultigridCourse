function T_new = initial_condition(T)
    [m, n] = size(T);
    T_new = T;

    for i = 2:m-1
        for j = 2:n-1
            % weighted average of boundaries:
            T_new(i, j) = (j * T(i, n) + (n + 1 - j) * T(i, 1) + ...
                (m + 1 - i) * T(1, j) + i * T(m, j)) / (m + n + 2);
        end
    end
end