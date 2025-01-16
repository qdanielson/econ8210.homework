% Plotting Dynare Policy fn, calculating residuals

%% Setup
% View endogenous variables (in declared order)
M_.endo_names
% Should be: y, c, k, l, z (1, 2, 3, 4, 5)

% View DR order:
oo_.dr.order_var
% Should be  1, 4, 3, 5, 2 (y, l, k, z, c)

% View state vars
M_.state_var
% Should be 5, 3 -- in our vector we will follow DR order (3, 5) = (k, z)

%% Define functions, pull steady-state values
% Extract steady states and coefficients
k_ss = oo_.dr.ys(3);
z_ss = oo_.dr.ys(5);

% Policy function for c
policy_c = @(k, z) oo_.dr.ys(2) ...
    + oo_.dr.ghx(5,:) * [k - k_ss; z - z_ss] ...
    + oo_.dr.ghxx(5,:) * kron([k - k_ss; z - z_ss], [k - k_ss; z - z_ss]) ...
    + oo_.dr.ghxxx(5,:) * kron([k - k_ss; z - z_ss], kron([k - k_ss; z - z_ss], [k - k_ss; z - z_ss]));

% Policy function for k
policy_k = @(k, z) oo_.dr.ys(3) ...
    + oo_.dr.ghx(3,:) * [k - k_ss; z - z_ss] ...
    + oo_.dr.ghxx(3,:) * kron([k - k_ss; z - z_ss], [k - k_ss; z - z_ss]) ...
    + oo_.dr.ghxxx(3,:) * kron([k - k_ss; z - z_ss], kron([k - k_ss; z - z_ss], [k - k_ss; z - z_ss]));

% Policy function for l
policy_l = @(k, z) oo_.dr.ys(4) ...
    + oo_.dr.ghx(2,:) * [k - k_ss; z - z_ss] ...
    + oo_.dr.ghxx(2,:) * kron([k - k_ss; z - z_ss], [k - k_ss; z - z_ss]) ...
    + oo_.dr.ghxxx(2,:) * kron([k - k_ss; z - z_ss], kron([k - k_ss; z - z_ss], [k - k_ss; z - z_ss]));

%% Plotting Part 1 -- calculate data

% Define ranges around steady states
k_range = linspace(k_ss - 2, k_ss + 4, 500);
z_range = linspace(-0.025, 0.025, 3);

[K, Z] = meshgrid(k_range, z_range);

% Initialize matrices
C = zeros(size(K));
K_next = zeros(size(K));
L = zeros(size(K));

% Compute for each grid point
for i = 1:numel(K)
    C(i) = policy_c(K(i), Z(i));
    K_next(i) = policy_k(K(i), Z(i));
    L(i) = policy_l(K(i), Z(i));
end

%% Plotting part 2 -- Make plots

% Plot for c
figure;
surf(K, Z, C);
xlabel('k');
ylabel('z');
zlabel('c');
title('Policy Function for c');
shading interp;

% Plot for k'
figure;
figure;
surf(K, Z, K_next);
xlabel('k');
ylabel('z');
zlabel('k''');
title('Policy Function for k''');
shading interp;

% Plot for l
figure;
surf(K, Z, L);
xlabel('k');
ylabel('z');
zlabel('l');
title('Policy Function for l');
shading interp;

%% Plotting part 3.A -- Euler Residuals

% Define parameter struct
params.delta = 0.1;
params.alpha_k = 0.33;
params.alpha_l = 0.67;
params.rho = 0.95;
params.sigma = 0.007;
params.beta = 0.97;

% reuse functions from previous parts:
U_l = @(l) -l;                    % Marginal Utility of labor
U_c = @(c) 1 / c;                 % Marginal Utility of consumption
inv_U_c = @(Uc) 1 / Uc;           % Inverse marginal utility of consumption
Y = @(z, k, l, params) exp(z) * k^params.alpha_k * l^params.alpha_l;
Y_l = @(z, k, l, params) params.alpha_l * exp(z) * k^params.alpha_k * l^(params.alpha_l - 1);
Y_k = @(z, k, l, params) params.alpha_k * exp(z) * k^(params.alpha_k - 1) * l^params.alpha_l;
find_c = @(l, z, k, params) inv_U_c(-U_l(l) / Y_l(z, k, l, params));
find_i = @(l, z, k, params) Y(z, k, l, params) - find_c(l, z, k, params);
find_k1 = @(l, z, k, params) (1 - params.delta) * k + find_i(l, z, k, params);

% restrict z = 0 for simplicity
C_resid = C(2,:);
K1_resid = K_next(2,:);
L_resid = L(2,:);

% Define the Euler Residual as an Anonymous Function
% This function computes the residual for each k
euler_residual_fn = @(c, k1, l) abs(U_c(c) - (1 + Y_k(0, k1, l, params) - params.delta) .* U_c(find_c(policy_l(k1,0), 0, k1, params)));

resid = zeros(500,1);

% Compute for each grid point
for i = 1:500
    resid(i) = euler_residual_fn(C_resid(i), K1_resid(i), L_resid(i));
end


%% Plotting part 3.B -- Euler Residuals (actual plots)

figure;
plot(k_range, resid, 'b-', 'LineWidth', 1.5);
xlabel('k');
ylabel('Euler Residual');
title('Euler Residuals vs. k (z = 0)');
grid on;
hold on;

% Add vertical line at steady-state k
xline(k_ss, '--r', sprintf('k_{ss} = %.2f', k_ss), 'LineWidth', 1.5);

hold off;