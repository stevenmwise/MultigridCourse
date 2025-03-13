%% Multilevel Function
function [T, ks, iterations, level_solutions] = multilevel(T, k_final, H1, H2)
    n = size(T, 1);
    k_start = floor(log2(n - 1));
    
    ks = [];
    iterations = [];
    level_solutions = {};
    
    for k = k_start:k_final
        ks(end+1) = k;
        
        results = simulation(T, 1e-5, 5000, H1, H2);
        iterations(end+1) = length(results);
        
        T = results{end};
        level_solutions{end+1} = T;  % Store the solution at this level
        
        if k < k_final
            T = uplevel(T, H1, H2);
            T = add_heaters(T, H1, H2);
        end
    end
end
