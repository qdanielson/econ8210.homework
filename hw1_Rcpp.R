## ECON8210
## Homework 1
## Written by Quinn Danielson, Doctoral Student Wharton Finance

#Clear environment
rm(list=ls())

start_all <- Sys.time()

#Load Rccp for fast looping
library(Rcpp)

#Set working directory
setwd("C:/Users/qdani/OneDrive/Documents/grad_school/coursework/ECON8210/homework/hw1/rcpp_test/")

#load Cpp functions
sourceCpp("numerical_analysis.cpp")

set.seed(2024)

## Question 1 - Numerical Integration

# Set parameters
t <- 100
rho <- .04
lambda <- .02

# Quadrature methods
# Set up output matrix, precisions

precision_list <- c(10,1,.1,.01,.001)

output_mat <- matrix(nrow = 4,ncol = length(precision_list))
time_mat <- matrix(nrow = 4,ncol = length(precision_list))

colnames(output_mat) <- precision_list
colnames(time_mat) <- precision_list

rownames(output_mat) <- c("Midpoint","Trapezoid","Simpson","Monte Carlo")
rownames(time_mat) <- c("Midpoint","Trapezoid","Simpson","Monte Carlo")

# Define R functions for MC analysis which doesn't use Rcpp

# Write functions for integral

gR <- function(x){exp(-rho*x)}
fR <- function(x){1 - exp(-lambda*x)}
uR <- function(x){-exp(-1*x)}

int <- function(x){gR(x)*uR(fR(x))}

# Begin loop over precisions

for (precision in precision_list) {
  
  t_seq <- seq(0,t,by=precision)
  
  mid_start <- Sys.time()
  
  V_midpoint <- midpoint_rule(t_seq,rho,lambda)
  
  mid_end <- Sys.time()
  
  # Trapezoid

  trap_start <- Sys.time()
  
  V_trapezoid <- trapezoid_rule(t_seq,rho,lambda)
  
  trap_end <- Sys.time()
  
  # Simpson Rule
  
  sim_start <- Sys.time()
  
  V_simpson <- trapezoid_rule(t_seq,rho,lambda)
  
  sim_end <- Sys.time()
  
  # Monte Carlo
  
  MC_start <- Sys.time()
  
  draws <- runif(100/precision,0,100)
  
  V_MC <- 100*mean(int(draws))
  
  MC_end <- Sys.time()
  
  outputV <- c(V_midpoint,V_trapezoid,V_simpson,V_MC)
  
  outputT <- c(mid_end - mid_start, trap_end - trap_start, sim_end - sim_start, MC_start - MC_end)
  
  output_mat[,which(precision == precision_list)] <- outputV
  time_mat[,which(precision == precision_list)] <- outputT
  
  rm(outputV)
  rm(outputT)
  
}

# Define x_vec, y_vec, and epsilon in R
precision <- .001
epsilon <- .000001

x_vec <- seq(-1, 3, by = precision)
y_vec <- seq(-1, 3, by = precision)

start <- Sys.time()

# Call the C++ function
result <- derivative_and_hessian(x_vec, y_vec, epsilon)

end <- Sys.time()

end - start

# Access results
V_mat <- result$V_mat
dV_array_x <- result$dV_array_x
dV_array_y <- result$dV_array_y
hessian_xx <- result$hessian_xx
hessian_xy <- result$hessian_xy
hessian_yx <- result$hessian_yx
hessian_yy <- result$hessian_yy

dV_array <- array(dim = c(dim(dV_array_x)[1], dim(dV_array_x)[2], 2)) # derivative array
dV_array[,,1] <- dV_array_x
dV_array[,,2] <- dV_array_y

HD_array <- array(dim = c(dim(hessian_xx)[1], dim(hessian_xx)[2], 2, 2)) # Hessian array
HD_array[,,1,1] <- hessian_xx
HD_array[,,1,2] <- hessian_xy
HD_array[,,2,1] <- hessian_yx
HD_array[,,2,2] <- hessian_yy

# Gradient descent

l <- length(x_vec)

alpha <- 1

xstart <- sample(1:length(x_vec),1)
ystart <- sample(1:length(y_vec),1)

point_mat <- matrix(nrow = 10000, ncol = 2)
value_mat <- matrix(nrow = 10000, ncol = 2)

point_mat[1,] <- c(xstart,ystart)

err <- sum(abs(point_mat[1,] - c(2001,2001)))

start <- Sys.time()

i <- 1

while (err > 1) {
  
  i <- i + 1
  
  point <- c(x_vec[point_mat[i-1,1]],y_vec[point_mat[i-1,2]])
  
  value_mat[i-1,] <- point
  
  d <- -(dV_array[point_mat[i-1,1],point_mat[i-1,2],]/sqrt(sum(dV_array[point_mat[i-1,1],point_mat[i-1,2],]^2)))
  
  movement_r <-  d*(10*alpha/i)
  movement_p <-  round(movement_r/precision,0)
  
  point_mat[i,] <- point_mat[i-1,] + movement_p
  
  # Make sure we are in-bounds
  point_mat[i,1] <- min(point_mat[i,1],l)
  point_mat[i,1] <- max(point_mat[i,1],1)
  
  point_mat[i,2] <- min(point_mat[i,2],l)
  point_mat[i,2] <- max(point_mat[i,2],1)
  
  err <- sum(abs(point_mat[i,] - c(2001,2001)))
}

end <- Sys.time()

t_grad <- end - start

l_grad <- i


# Conjugate Descent

alpha <- 1
beta <- 0
d_lag <- matrix(c(0,0),nrow = 2, ncol = 1)

xstart <- sample(1:length(x_vec),1)
ystart <- sample(1:length(y_vec),1)

point_mat <- matrix(nrow = 1000, ncol = 2)
value_mat <- matrix(nrow = 1000, ncol = 2)

point_mat[1,] <- c(xstart,ystart)

err <- sum(abs(point_mat[1,] - c(2001,2001)))

i <- 1

start <- Sys.time()

while (err > 1) {
  
  i <- i + 1
  
  point <- c(x_vec[point_mat[i-1,1]],y_vec[point_mat[i-1,2]])
  
  value_mat[i-1,] <- point
  
  d <- -(dV_array[point_mat[i-1,1],point_mat[i-1,2],]) + beta[[1]] * d_lag
  
  movement_r <-  d/sqrt(sum(d^2))*(alpha/i)
  movement_p <-  round(movement_r/precision,0)
  
  point_mat[i,] <- point_mat[i-1,] + movement_p
  
  # Make sure we are in-bounds
  point_mat[i,1] <- min(point_mat[i,1],l)
  point_mat[i,1] <- max(point_mat[i,1],1)
  
  point_mat[i,2] <- min(point_mat[i,2],l)
  point_mat[i,2] <- max(point_mat[i,2],1)
  
  d_lead <- -(dV_array[point_mat[i,1],point_mat[i,2],])
  
  beta <- (t(d_lead)%*%d_lead)/(t(d)%*%d)
  
  d_lag <- d
  
  rm(d)
  rm(d_lead)
  
  err <- sum(abs(point_mat[i,] - c(2001,2001)))
  
}

end <- Sys.time()

end - start

t_conj <- end - start

l_conj <- i

# Newton-Raphson

xstart <- sample(1:length(x_vec),1)
ystart <- sample(1:length(y_vec),1)

point_mat <- matrix(nrow = 1000, ncol = 2)
value_mat <- matrix(nrow = 1000, ncol = 2)

point_mat[1,] <- c(xstart,ystart)

err <- sum(abs(point_mat[1,] - c(2001,2001)))

start <- Sys.time()

i <- 1

while (err > 1) {
  
  i <- i + 1
  
  point <- c(x_vec[point_mat[i-1,1]],y_vec[point_mat[i-1,2]])
  
  value_mat[i-1,] <- point
  
  H <- HD_array[point_mat[i-1,1],point_mat[i-1,2],,]
  
  d <- solve(H) %*% dV_array[point_mat[i-1,1],point_mat[i-1,2],]
  
  movement_r <-  -d
  movement_p <-  round(movement_r/precision,0)
  
  point_mat[i,] <- point_mat[i-1,] + movement_p
  
  # Make sure we are in-bounds
  point_mat[i,1] <- min(point_mat[i,1],l)
  point_mat[i,1] <- max(point_mat[i,1],1)
  
  point_mat[i,2] <- min(point_mat[i,2],l)
  point_mat[i,2] <- max(point_mat[i,2],1)
  
  rm(d)
  
  err <- sum(abs(point_mat[i,] - c(2001,2001)))
  
}

end <- Sys.time()

end - start

t_nr <- end - start

l_nr <- i

# BFGS

xstart <- sample(1:length(x_vec),1)
ystart <- sample(1:length(y_vec),1)

point_mat <- matrix(nrow = 1000, ncol = 2)
value_mat <- matrix(nrow = 1000, ncol = 2)

point_mat[1,] <- c(xstart,ystart)

err <- sum(abs(point_mat[1,] - c(2001,2001)))

H <- HD_array[point_mat[1,1],point_mat[1,2],,]
H_inv <- solve(H)

start <- Sys.time()

i <- 1

while (err > 1) {
  
  i <- i + 1
  
  point <- c(x_vec[point_mat[i-1,1]],y_vec[point_mat[i-1,2]])
  
  value_mat[i-1,] <- point
  
  d <- H_inv %*% dV_array[point_mat[i-1,1],point_mat[i-1,2],]
  
  movement_r <-  -d
  movement_p <-  round(movement_r/precision,0)
  
  if (sum(movement_p == c(0,0)) == 2) {            # Adjust point slightly so it does not get stuck
    movement_p[,1] <- sample(c(-1,1),2)
      }
  
  point_mat[i,] <- point_mat[i-1,] + movement_p
  
  # Make sure we are in-bounds
  point_mat[i,1] <- min(point_mat[i,1],l)
  point_mat[i,1] <- max(point_mat[i,1],1)
  
  point_mat[i,2] <- min(point_mat[i,2],l)
  point_mat[i,2] <- max(point_mat[i,2],1)
  
  new_point <- c(x_vec[point_mat[i,1]],y_vec[point_mat[i,2]])
  
  s <- new_point - point
  g <- dV_array[point_mat[i,1],point_mat[i,2],] - dV_array[point_mat[i-1,1],point_mat[i-1,2],]
  
  H_inv_new <- (diag(2) - (s %*% t(g))/(t(g) %*% s)[[1]]) %*% H_inv %*% (diag(2) - (g %*% t(s))/(t(g) %*% s)[[1]]) + ((s %*% t(s)))/((t(g) %*% s)[[1]])
  
  H_inv <- H_inv_new
  
  if(sum((point - new_point) == c(0,0)) == 2) {    # If stuck at edge or corner, pick new point
    
    point_mat[i,] <- c(sample(1:length(x_vec),1),
                       sample(1:length(y_vec),1))
    
    H <- HD_array[point_mat[i,1],point_mat[i,2],,]
    H_inv <- solve(H)
    
  }
  
  rm(d)
  rm(H_inv_new)
  rm(s)
  rm(g)
  
  err <- sum(abs(point_mat[i,] - c(2001,2001)))
  
}

end <- Sys.time()

t_bfgs <- end - start

l_bfgs <- i

output_opt <- matrix(c(l_grad,l_conj,l_nr,l_bfgs,t_grad,t_conj,t_nr,t_bfgs),nrow = 4,ncol = 2)

colnames(output_opt) <- c("Iterations","Time")
rownames(output_opt) <- c("Gradient Descent","Conjugate Descent","Newton-Raphson","BFGS")

## Problem 3 - Grid Search Optimization (INCOMPLETE)

# Initialize grid

precision <- .0001
split <- seq(0,1,precision)

create_params <- function(agent_count,
                          benchmark) {

  if (benchmark == TRUE) {
    
    alpha <- rep(1,agent_count)
    lambda <- rep(1,agent_count)
    omega <- matrix(rep(-.5,agent_count*agent_count),nrow = agent_count)
    
  } else {
    
    alpha <- runif(agent_count) + 1
    lambda <- runif(agent_count) + 1
    omega <- matrix(-(runif(agent_count*agent_count) + .5),nrow = agent_count)
    
  }
  
  return(list(alpha = alpha,lambda = lambda,omega = omega))
}

start <- Sys.time()

e_mat <- matrix(1,nrow = 3,ncol = 3)

x_mat.e <- matrix(ncol = 3, nrow = 3)
x_mat.v <- matrix(ncol = 3, nrow = 3)

params_even <- create_params(agent_count = 3, benchmark = TRUE)

params_vary <- create_params(agent_count = 3, benchmark = FALSE)

alpha_even <- params_even$alpha
omega_even <- params_even$omega
lambda_even <- params_even$lambda

alpha_vary <- c(.9,1,1.1)
omega_vary <- -matrix(c(.4,.5,.6,
                        .6,.5,.4,
                        .5,.4,.6),
                      nrow = 3,
                      ncol = 3,
                      byrow = TRUE)
lambda_vary <- c(.9,1,1.1)

for (i in 0:2) {
  
  output <- grid_search(e_mat,alpha = alpha_even,omega = omega_even,lambda = lambda_even,split = split,index = i)
  
  max.row <- output[which(output[,3] == max(output[,3])),]
  
  e <- sum(e_mat[,i+1])
  
  x.col <- c(e * max.row[1],
            (e - (e * max.row[1])) * max.row[2],
            (e - (e * max.row[1])) * (1 - max.row[2]))
  
  x_mat.e[,i+1] <- x.col
  
}

for (i in 0:2) {
  
  output <- grid_search(e_mat,alpha = alpha_vary,omega = omega_vary,lambda = lambda_vary,split = split,index = i)
  
  max.row <- output[which(output[,3] == max(output[,3])),]
  
  e <- sum(e_mat[,i+1])
  
  x.col <- c(e * max.row[1],
             (e - (e * max.row[1])) * max.row[2],
             (e - (e * max.row[1])) * (1 - max.row[2]))
  
  x_mat.v[,i+1] <- x.col
  
}



end <- Sys.time()

end-start

# Demonstrate why brute force grid search will not work for 10 agent-10 good case

V3 <- array(data = runif(10^3),dim = rep(10,3))

object.size(V3)

# V10 <- array(data = runif(10^10),dim = rep(10,10)) # Error: cannot allocate vector of size 74.5 Gb

## One way to do it would be to nest loops and solve "recursively", that is, nest 9 loops (agents - 1, since the consumption of the last agent is pinned down by the other 9) and within each loop keep only the maximum value. This ensures that you only have to keep a vector of size (10001 in the 3-agent case) at all times. But this would require very careful coding to avoid memory issues and still requires 10001^9 calculations. An alternative approach involving solving the analytic derivatives by hand (this function has no mixed second derivatives) and doing gradient or conjugate descent would likely work, and would allow you to evaluate the function only at the necessary points.

## 4 - Solve competitive equilibrium

# a * (x^1+om)/1+om -> 1/a * (p*theta)^(1/om) = x

alpha_mat <- matrix(alpha_even,3,3,byrow = T)

x_mat.e <- round(x_mat.e,3)
x_mat.v <- round(x_mat.v,3)

f.even <- function(point) {
  
  x <- x_mat.e # values from x_mat
  
  alpha_mat <- matrix(alpha_even,3,3,byrow = T) # Values from alpha
  
  e_mat <- matrix(1,3,3) # values set
  
  omega <- omega_even
  
  p <- c(1,point[1:2]) # normalize price 1 to 1
  
  theta <- point[3:5]
  
  ptheta <- p %*% t(theta)
  
  x.guess <- (1/alpha_mat) * ptheta^omega
  
  x.error <- x - x.guess
  
  p.error <- (e_mat %*% p) - (x.guess %*% p)
  
  return(c(x.error,p.error))
  
}

f.vary <- function(point) {
  
  x <- x_mat.v # values from x_mat
  
  alpha_mat <- matrix(alpha_vary,3,3,byrow = T) # Values from alpha
  
  e_mat <- matrix(1,3,3) # values set
  
  omega <- omega_vary
  
  p <- c(1,point[1:2]) # normalize price 1 to 1
  
  theta <- point[3:5]
  
  ptheta <- p %*% t(theta)
  
  x.guess <- (1/alpha_mat) * ptheta^omega
  
  x.error <- colSums(x - x.guess)
  
  p.error <- (e_mat %*% p) - (x.guess %*% p)
  
  return(c(x.error,p.error))
  
}

fsolve(f.even,c(1.1,   # price 2
                1.1,   # price 3
                1.1,   # theta 1
                1.1,   # theta 2
                1.1))  # theta 3

fsolve(f.vary,c(1,   # price 2
                1,   # price 3
                1,   # theta 1
                1,   # theta 2
                1))  # theta 3

end_all <- Sys.time()

ttotal <- end_all - start_all
  













