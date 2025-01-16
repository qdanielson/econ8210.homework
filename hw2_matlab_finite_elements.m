% Set working directory
cd('C:/Users/qdani/OneDrive/Documents/grad_school/coursework/ECON8210/homework/hw2/');
rng(2025);

% Define parameter struct
params.delta = 0.1;
params.alpha_k = 0.33;
params.alpha_l = 0.67;
params.rho = 0.95;
params.sigma = 0.007;
params.beta = 0.97;


% Steady-state calculation
x0 = [3, 0.9]';

% Use an anonymous function to pass `params` into `steady_state`
ss_values = fsolve(@(x) steady_state(x, params), x0);

% Compute steady-state consumption
[~, css] = steady_state(ss_values, params);

% Extract steady-state values for k and l
k_ss = ss_values(1);
l_ss = ss_values(2);

% Discretize distribution
[z_grid,P] = discretizeAR1_Tauchen(0,params.rho,params.sigma,3,1);

% Set grids
k_min = k_ss - 2;
k_max = k_ss + 4;

z_min = min(z_grid);
z_max = max(z_grid);

% Set no. of elements on capital, productivity
n_k = 6;
n_z = 3;

k_grid = linspace(k_min,k_max,n_k);
z_grid = linspace(z_min,z_max,n_z);

first_weights = l_ss * ones(n_k*n_z,1);

% Number of samples
n = 100000;

% Generate 100 random k values in the range [k_min, k_max]
k_points = k_min + (k_max - k_min) * rand(1, n);

% Generate 100 random integers between 1 and 3
z_points = randi([1, 3], 1, n);

% Define the residual function
solve_fe = @(weights) residual_fe(weights, k_points, z_points, n, k_grid, z_grid, P, params);

% Set up options for lsqnonlin
options = optimoptions('lsqnonlin', ...
    'Display', 'iter', ...          % Display iteration progress
    'MaxIterations', 1000, ...      % Increase maximum iterations
    'FunctionTolerance', 1e-6, ...  % Adjust tolerances as needed
    'StepTolerance', 1e-6);

% Solve for optimal weights
optimal_weights = lsqnonlin(solve_fe, first_weights, [], [], options);

% See final error
solve_fe(first_weights)
solve_fe(optimal_weights)

% Anonymous functions
    U_l = @(l) -l;                    % Marginal Utility of labor
    U_c = @(c) 1 / c;                 % Marginal Utility of consumption
    inv_U_c = @(Uc) 1 / Uc;           % Inverse marginal utility of consumption
    Y = @(z, k, l, params) exp(z) * k^params.alpha_k * l^params.alpha_l;
    Y_l = @(z, k, l, params) params.alpha_l * exp(z) * k^params.alpha_k * l^(params.alpha_l - 1);
    Y_k = @(z, k, l, params) params.alpha_k * exp(z) * k^(params.alpha_k - 1) * l^params.alpha_l;
    find_c = @(l, z, k, params) inv_U_c(-U_l(l) / Y_l(z, k, l, params));
    find_i = @(l, z, k, params) Y(z, k, l, params) - find_c(l, z, k, params);
    find_k1 = @(l, z, k, params) (1 - params.delta) * k + find_i(l, z, k, params);
    
k_grid_test = linspace(k_min,k_max,100);
n_k_grid_test = length(k_grid_test);
l1 = zeros(n_k_grid_test, n_z);
c = zeros(n_k_grid_test, n_z);
k1 = zeros(n_k_grid_test, n_z);
for g = 1:n_k_grid_test
    for zp = 1:n_z
        % Calculate elements
        elements = labor_interpolate(k_grid_test(g), zp, k_grid, z_grid);
        l1(g,zp) = approx_labor(elements, optimal_weights);
        c(g,zp) = find_c(l1(g,zp),z_grid(zp),k_grid_test(g),params);
        k1(g,zp) = find_k1(l1(g,zp),z_grid(zp),k_grid_test(g),params);
    end
end



% Assuming c and k1 are matrices of size (100, n_z)
% k_grid_test is a vector of size (100, 1)
% z_grid is a vector of size (1, n_z)

% Plot consumption policy function
figure;
hold on;
for zp = 1:n_z
    plot(k_grid_test, c(:, zp), 'DisplayName', ['z\_grid[', num2str(zp), '] = ', num2str(z_grid(zp))]);
end
title('Consumption Policy Function');
xlabel('Capital (k)');
ylabel('Consumption (c)');
legend show;
grid on;
hold off;

% Plot capital policy function with 45-degree line
figure;
hold on;

% Plot the policy function lines
for zp = 1:n_z
    plot(k_grid_test, k1(:, zp), 'DisplayName', ['z\_grid[', num2str(zp), '] = ', num2str(z_grid(zp))]);
end

% Plot the 45-degree line
plot(k_grid_test, k_grid_test, '--k', 'DisplayName', '45-degree line');

title('Capital Policy Function');
xlabel('Capital (k)');
ylabel('Next-period Capital (k1)');
legend show;
grid on;
hold off;

% Initialize residuals matrix
residuals = zeros(n_k_grid_test, n_z);

% Compute residuals for each value of k and z
for g = 1:n_k_grid_test
    for zp = 1:n_z
        residuals(g, zp) = residual_fe_test(optimal_weights, k_grid_test(g), zp, k_grid, z_grid, P, params);
    end
end

% Plot the residuals
figure;
hold on;
for zp = 1:n_z
    plot(k_grid_test, residuals(:, zp), 'DisplayName', ['z\_grid[', num2str(zp), '] = ', num2str(z_grid(zp))]);
end
title('Residuals for test k Grid at Different z');
xlabel('Capital (k)');
ylabel('Residual (res)');
legend show;
grid on;
hold off;

