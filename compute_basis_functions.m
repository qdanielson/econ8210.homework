function basis_functions = compute_basis_functions(grid, x)
    % Ensure x is within the bounds of the grid
    if x < grid(1) || x > grid(end)
        error('x is outside the bounds of the grid.');
    end

    % Initialize the basis function vector
    basis_functions = zeros(1, length(grid));

    % Define the basis function for each segment
    for i = 1:length(grid)
        if i == 1  % First basis function
            if x >= grid(i) && x <= grid(i + 1)
                basis_functions(i) = (grid(i + 1) - x) / (grid(i + 1) - grid(i));
            else
                basis_functions(i) = 0;
            end
        elseif i == length(grid)  % Last basis function
            if x >= grid(i - 1) && x <= grid(i)
                basis_functions(i) = (x - grid(i - 1)) / (grid(i) - grid(i - 1));
            else
                basis_functions(i) = 0;
            end
        else  % Middle basis functions
            if x >= grid(i - 1) && x < grid(i)
                basis_functions(i) = (x - grid(i - 1)) / (grid(i) - grid(i - 1));
            elseif x >= grid(i) && x <= grid(i + 1)
                basis_functions(i) = (grid(i + 1) - x) / (grid(i + 1) - grid(i));
            else
                basis_functions(i) = 0;
            end
        end
    end
end
