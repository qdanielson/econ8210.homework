function [res, c] = steady_state(x, params)
    % Set values, z is hard-coded
    z = 0;
    k = x(1);
    l = x(2);
    % Unpack params
    delta = params.delta;
    beta = params.beta;
    % define anonymous temporary functions
    % Utility functions
    U_c = @(c) 1 / c;
    
    % Production function
    Y = @(z, k, l, params) exp(z) * (k^params.alpha_k) * (l^params.alpha_l);
    Y_k = @(z, k, l, params) params.alpha_k * exp(z) * (k^(params.alpha_k - 1)) * (l^params.alpha_l);
    Y_l = @(z, k, l, params) params.alpha_l * exp(z) * (k^params.alpha_k) * (l^(params.alpha_l - 1));

    % Calculate Consumption
    c = Y(z,k,l,params) - delta * k;
    % Set optimality conditions
    res(1) = 1 - beta * (Y_k(z,k,l,params) + 1 - delta);
    res(2) = l - (Y_l(z,k,l,params) * U_c(c));
end