function [res] = residual_chebyshev(coefficients, T, chebyshev_function, k_nodes, z_grid, P, params)
    % Unpack params
    delta = params.delta;
    beta = params.beta;
    n_k = length(k_nodes);
    n_z = length(z_grid);
    coefficients = reshape(coefficients, n_k, n_z);
    % Anonymous functions
    U_l = @(l) -l;                    % Marginal Utility of labor
    U_c = @(c) 1 / c;                 % Marginal Utility of consumption
    inv_U_c = @(Uc) 1 / Uc;           % Inverse marginal utility of consumption
    Y = @(z, k, l, params) exp(z) * k^params.alpha_k * l^params.alpha_l;
    Y_l = @(z, k, l, params) params.alpha_l * exp(z) * k^params.alpha_k * l^(params.alpha_l - 1);
    Y_k = @(z, k, l, params) params.alpha_k * exp(z) * k^(params.alpha_k - 1) * l^params.alpha_l;
    find_c = @(l, z, k, params) inv_U_c(-U_l(l) / Y_l(z, k, l, params));
    find_i = @(l, z, k, params) Y(z, k, l, params) - find_c(l, z, k, params);
    find_k1 = @(l, z, k, params) (1 - delta) * k + find_i(l, z, k, params);
    
    % Initialize outputs
    residuals = zeros(n_k, n_z);
    
    % Double loop over state grid
    for i = 1:n_k
        for j = 1:n_z
            z = z_grid(j);
            k = k_nodes(i);
            l = T(1,:)*coefficients(:,1);
            % Define r for euler eq. + copntrol equations
            r = Y_k(z, k, l, params);
            c = find_c(l, z, k, params);
            k1 = find_k1(l, z, k, params);

            % Evaluate next-period labor for each z'
            l1 = zeros(1, n_z);
            for zp = 1:n_z
                l1(zp) = funeval(coefficients(:, zp), chebyshev_function, k1);
            end
            % solve for next period marginal utility
            % Calculate marginal utility in next period
            marginal_value = 0;
            for zp = 1:n_z
                marginal_value = marginal_value + P(1, zp) * U_c(find_c(l1(zp), z_grid(zp), k1, params));
            end
            residuals(i, j) = (U_c(c) - beta * (1 - delta + r) * marginal_value)^2;
        end
    end
    res = sum(residuals, 'all')
end