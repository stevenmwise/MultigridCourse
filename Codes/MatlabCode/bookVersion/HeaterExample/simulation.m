%% Simulation Function
function results = simulation(T, epsilon, num_steps,H1,H2)
    results = {T};

    for i = 1:num_steps
        T_prime = jacobi_step(T);
        T_prime = add_heaters(T_prime,H1,H2);

        results{end+1} = T_prime;
        
        if max(abs(T_prime(:) - T(:))) < epsilon
            return;
        end

        T = T_prime;
    end
end
