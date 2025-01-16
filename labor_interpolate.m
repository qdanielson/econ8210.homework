function [elements] = labor_interpolate(k,z_index,k_grid,z_grid)
    n_k = length(k_grid);
    n_z = length(z_grid);

    % Create a 3-element vector initialized with zeros
    z_basis = zeros(1, 3);

    % Set the value at the random index to 1
    z_basis(z_index) = 1;

    k_basis = compute_basis_functions(k_grid,k);

    elements = reshape(k_basis' * z_basis, n_k*n_z,1);
end