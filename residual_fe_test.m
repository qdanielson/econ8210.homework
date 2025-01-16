function [res] = residual_fe_test(weights, k, z_index, k_grid, z_grid, P, params)
    % Unpack params
    delta = params.delta;
    beta = params.beta;
    n_k = length(k_grid);
    n_z = length(z_grid);
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

    % Calculate elements
    elements = labor_interpolate(k, z_index, k_grid, z_grid);
    l = approx_labor(elements, weights);
    z = z_grid(z_index);

    % Define r for Euler equation and control equations
    r = Y_k(z, k, l, params);
    c = find_c(l, z, k, params);
    k1 = find_k1(l, z, k, params); % Compute k1
    k1 = max(k_grid(1), min(k1, k_grid(end))); % Restrict k1 to grid bounds
            
    % Calculate next-period marginal utility
    marginal_value = 0;
    for zp = 1:n_z
        l1 = approx_labor(labor_interpolate(k1, zp, k_grid, z_grid), weights);
        marginal_value = marginal_value + P(z_index, zp) * U_c(find_c(l1, z_grid(zp), k1, params));
    end

    % Residual for this point
    res = U_c(c) - beta * (1 - delta + r) * marginal_value;
end