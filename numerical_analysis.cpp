#include <Rcpp.h>
using namespace Rcpp;

// Define g, f, u, and integrand functions

// [[Rcpp::export]]
double g(double x, double rho) {  
  return exp(-rho * x);
}

// [[Rcpp::export]]
double f(double x, double lambda) {
  return 1 - exp(-lambda * x);
}

// [[Rcpp::export]]
double u(double x) {
  return -exp(-x);
}

// [[Rcpp::export]]
double integrand(double x, double rho, double lambda) { // This function combines the three other functions into a single integrand. I split them up for tractability.
  return g(x, rho) * u(f(x, lambda));
}

// Midpoint rule function
// [[Rcpp::export]]
double midpoint_rule(NumericVector t_seq, double rho, double lambda) {
  // Initialize variables
  double V_midpoint = 0.0;
  
  // Midpoint loop
  for (int i = 1; i < t_seq.size(); i++) {
    double mid_point = (t_seq[i] + t_seq[i-1]) / 2;
    double slice = (t_seq[i] - t_seq[i-1]) * integrand(mid_point, rho, lambda); // pass rho and lambda
    V_midpoint += slice;
  }
  
  return V_midpoint;
  }
  
  // Trapezoid rule function
  // [[Rcpp::export]]
  double trapezoid_rule(NumericVector t_seq, double rho, double lambda) {
    // Initialize variables
    double V_trapezoid = 0.0;
    
    // Midpoint loop
    for (int i = 1; i < t_seq.size(); i++) {
      double mid_point = (integrand(t_seq[i], rho, lambda) + integrand(t_seq[i - 1], rho, lambda))/2;
      double slice = (t_seq[i] - t_seq[i-1]) * mid_point; 
      V_trapezoid += slice;
    }
    
    return V_trapezoid;
}

// Simpson rule function
// [[Rcpp::export]]
double simpson_rule(NumericVector t_seq, double rho, double lambda) {
  // Initialize variables
  double V_simpson = 0.0;
  
  // Simpson loop
  for (int i = 1; i < t_seq.size(); i++) {
    double simpson = (integrand(t_seq[i], rho, lambda) + 4*integrand((t_seq[i] + t_seq[i - 1])/2, rho, lambda) + integrand(t_seq[i - 1], rho, lambda))/6; //note that the formula is DIFFERENT than that found in the 
    double slice = (t_seq[i] - t_seq[i-1]) * simpson; 
    V_simpson += slice;
  }
  
  return V_simpson;
}

// Objective Function
// [[Rcpp::export]]
double obj_function(double x, double y) {
  return 100 * (y - (x * x)) * (y - (x * x)) + (1 - x) * (1 - x);
}


// Question 3
// Derivative and Hessian Calculation
// [[Rcpp::export]]
List derivative_and_hessian(NumericVector x_vec, NumericVector y_vec, double epsilon) {
  
  // Declare matrices for storing values
  NumericMatrix V_mat(x_vec.size(), y_vec.size());
  NumericMatrix dV_array_x(x_vec.size(), y_vec.size()); // First derivative wrt x
  NumericMatrix dV_array_y(x_vec.size(), y_vec.size()); // First derivative wrt y
  NumericMatrix hessian_xx(x_vec.size(), y_vec.size()); // Second derivative wrt x (dx^2)
  NumericMatrix hessian_xy(x_vec.size(), y_vec.size()); // Mixed second derivative (dxdy)
  NumericMatrix hessian_yx(x_vec.size(), y_vec.size()); // Mixed second derivative (dydx)
  NumericMatrix hessian_yy(x_vec.size(), y_vec.size()); // Second derivative wrt y (dy^2)
  
  // Loop over x and y vectors
  for (int i = 0; i < x_vec.size(); i++) {    
    for (int j = 0; j < y_vec.size(); j++) {  
      
      double x = x_vec[i];
      double y = y_vec[j];
      
      // Calculate the value of the objective function at the point (x, y)
      V_mat(i, j) = obj_function(x, y);
      
      // First derivatives (numerical differentiation)
      dV_array_x(i, j) = (obj_function(x + epsilon, y) - V_mat(i, j)) / epsilon;
      dV_array_y(i, j) = (obj_function(x, y + epsilon) - V_mat(i, j)) / epsilon;
      
      // Second derivatives (numerical differentiation)
      double dx_ahead_for_dx2 = (obj_function(x + 2 * epsilon, y) - obj_function(x + epsilon, y)) / epsilon;
      double dy_ahead_for_dy2 = (obj_function(x, y + 2 * epsilon) - obj_function(x, y + epsilon)) / epsilon;
      double dy_at_dx = (obj_function(x + epsilon, y + epsilon) - obj_function(x, y + epsilon)) / epsilon;
      double dx_at_dy = (obj_function(x + epsilon, y + epsilon) - obj_function(x + epsilon, y)) / epsilon;
      
      // Update the Hessian (second derivatives)
      hessian_xx(i, j) = (dx_ahead_for_dx2 - dV_array_x(i, j)) / epsilon;
      hessian_yy(i, j) = (dy_ahead_for_dy2 - dV_array_y(i, j)) / epsilon;
      
      // Mixed second derivatives
      hessian_xy(i, j) = (dy_at_dx - dV_array_x(i, j)) / epsilon;
      hessian_yx(i, j) = (dx_at_dy - dV_array_y(i, j)) / epsilon;
    }
  }
  
  // Return the results as a list
  return List::create(Named("V_mat") = V_mat,
                      Named("dV_array_x") = dV_array_x,
                      Named("dV_array_y") = dV_array_y,
                      Named("hessian_xx") = hessian_xx,
                      Named("hessian_xy") = hessian_xy,
                      Named("hessian_yx") = hessian_yx,
                      Named("hessian_yy") = hessian_yy);
}

// Question 4

// Utility function u
double u3(double x, double alpha, double omega) { // named to differentiate from other u function
  return alpha * pow(x, 1 + omega) / (1 + omega);
}

// Grid search function
// [[Rcpp::export]]
NumericMatrix grid_search(NumericMatrix e_mat, NumericVector alpha, NumericMatrix omega, NumericVector lambda, NumericVector split, int index) {
  
  // Initialize the value matrix and candidates matrix
  NumericMatrix V(split.size(), split.size());
  NumericMatrix candidates(split.size(), 3);
  
  // sum up endowment by column (represents goods)
  
  double e = 0; // initialize as empty variable
  
  for (int k = 0; k < 3; k++) {
    e += e_mat(k, index);
  }
  
  // Loop over the split vector
  for (int i = 0; i < split.size(); i++) {
    for (int j = 0; j < split.size(); j++) {
      
      // Calculate agent shares of endowment good
      double x1 = split[i] * e;
      double x2 = (e - x1) * split[j];
      double x3 = e - x1 - x2;
      
      // Calculate value in state using the utility function u
      V(i, j) = lambda[0] * u3(x1, alpha[index], omega(index, 0)) +
                lambda[1] * u3(x2, alpha[index], omega(index, 1)) +
                lambda[2] * u3(x3, alpha[index], omega(index, 2));
    }
    
    // Find the index of the maximum value in the current row of V
    NumericVector row_values = V(i, _);
    int max_index = which_max(row_values);
    
    // Store the candidate values (split[i], split[max_index], max(V[i,]))
    candidates(i, 0) = split[i];
    candidates(i, 1) = split[max_index];
    candidates(i, 2) = max(row_values);
  }
  
  return candidates;
}