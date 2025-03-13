%% Add Heaters Function
function T = add_heaters(T,H1,H2)
    n = size(T, 1) - 1;
    
    T((3*n/8+1):(n/2+1), 1:(n/8+1)) = H1;
    T((n/2+1):(6*n/8+1), (7*n/8+1):n+1) = H2;
end
