%% Uplevel Function
function T_new = uplevel(T, H1, H2)
    n = size(T, 1);
    k = floor(log2(n - 1));
    T_new = grid_with_boundary(k + 1, H1, H2);
    T_new = add_heaters(T_new, H1, H2);

    new_n = size(T_new, 1);

    for i = 2:new_n-1
        for j = 2:new_n-1
            T_new(i, j) = (T(floor((i+1)/2), floor((j+1)/2)) + ...
                T(ceil((i+1)/2), floor((j+1)/2)) + ...
                T(floor((i+1)/2), ceil((j+1)/2)) + ...
                T(ceil((i+1)/2), ceil((j+1)/2))) / 4;
        end
    end
end
