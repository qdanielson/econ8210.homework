/*
 * Example 1 from F. Collard (2001): "Stochastic simulations with DYNARE:
 * A practical guide" (see "guide.pdf" in the documentation directory).
 */

/*
 * Copyright © 2001-2010 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
 */


var y, c, k, l, z;
varexo e;

parameters delta, alpha_k, alpha_l, rho, sigma, beta;

delta = 0.1;
alpha_k = 0.33;
alpha_l = 0.67;
rho = 0.95;
sigma = 0.007;
beta = 0.97;

model;
c * l = alpha_l * exp(z) * (k(-1)^alpha_k) * (l^(alpha_l - 1));
1 = beta * (c / c(+1)) * (alpha_k * exp(z) * (k^(alpha_k - 1)) * (l^alpha_l) + 1 - delta);
y = exp(z)*(k(-1)^alpha_k)*(l^alpha_l);
k = y - c + (1 - delta)*k(-1);
z = rho*z(-1) + e;
end;

initval;
y = 1.4923;
c = 1.1161;
l = 0.9465;
k = 3.7612;
z = 0;
e = 0;
end;

shocks;
var e; stderr sigma;

end;

stoch_simul(order=3, pruning, irf=20, periods=0);

