## ECON8210
## Homework 1
## Written by Quinn Danielson, Doctoral Student Wharton Finance

#Clear environment
rm(list=ls())

start_all <- Sys.time()

#Set working directory
setwd("C:/Users/qdani/OneDrive/Documents/grad_school/coursework/ECON8210/homework/hw1/rcpp_test/")

library(pracma)


## Question 6 -- VFI

# Parameters
params <- c(aalpha = 0.33, bbeta = 0.97, ddelta = 0.1, pphi = 0.05, ggamma = 0.2)
aalpha <- params["aalpha"]
bbeta <- params["bbeta"]
ddelta <- params["ddelta"]
pphi <- params["pphi"]
ggamma <- params["ggamma"]
ttau <- 0.25
z <- 0.0

# Define functions for control variables 

rate <- function(aalpha,z,k,l){
  r <- aalpha * exp(z) * k^(aalpha - 1) * l^(1 - aalpha)
  return(r)
}

wage <- function(aalpha,z,k,l){
  w <- (1 - aalpha) * exp(z) * k^(aalpha) * l^(-aalpha)
  return(w)
}

govt <- function(ttau,w,l){
  g <- ttau * w * l
  return(g)
}

consumption <- function(ttau,w,r,k,l,i){
  c <- (1 - ttau) * w * l + r * k - i
  return(c)
}

inv_ss <- function(k,ddelta){ # (STEADY STATE VERSION)
  i <- k * ddelta
  return(i)
}

adjustment.costs <- function(pphi,il,i){
  ac <- pphi * ((i/il) - 1)^2 * i
  return(ac)
}

k.next <- function(ddelta,k,i,ac){
  k1 <- i + (1 - ddelta) * k - ac
  return(k1)
}

## Define and solve steady-state

sseq <- function(vars){
  
  z <- 0
  ttau <- .25
  aalpha <- 1/3
  ddelta <- .1
  bbeta <- .97
  
  # unpack vars
  k <- vars[1]
  l <- vars[2]
  
  # Calculate control variables
  r <- rate(aalpha,z,k,l)
  w <- wage(aalpha,z,k,l)
  i <- inv_ss(k,ddelta)
  c <- consumption(ttau,w,r,k,l,i)
  
  # FOC errors (objective function wrt. k, l)
  e1 <- 1 - bbeta * (r + 1 - ddelta)
  e2 <- (1/c) * (1 - ttau) * w - l
  
  return(c(e1,e2))
}

ss.soln <- fsolve(sseq,c(1,1))$x

k.ss <- ss.soln[1]
i.ss <- ddelta * k.ss
l.ss <- ss.soln[2]
r.ss <- rate(aalpha,z = 0,k.ss,l.ss)
w.ss <- wage(aalpha,z = 0,k.ss,l.ss)
c.ss <- consumption(ttau = .25,w.ss,r.ss,k.ss,l.ss,i.ss)
g.ss <- govt(ttau = .25,w.ss,l.ss)

V.ss <- log(c.ss) + ggamma * log(g.ss) - (1/2) * l.ss^2

## Bellman function (incomplete)

make.ttau.z.transition.matrix <- function(){
  
  ttau.levels <- 1:3
  z.levels <- 1:5
  
  ttau.grid <- c(0.2,0.25,0.3)
  z.grid <- c(-0.0673,-0.0336,0,0.0336,0.0673)
  
  ttau.probabilities <- matrix(c(0.9,0.1,0,
                                 0.05,0.9,0.05,
                                 0,0.1,0.9),
                               nrow = 3,
                               ncol = 3,
                               byrow = TRUE)
  
  z.probabilities <- matrix(c(0.9727,0.0273,0,0,0,
                              0.0041,0.9806,0.0153,0,0,
                              0,0.0082,0.9836,0.0082,0,
                              0,0,0.0153,0.9806,0.0041,
                              0,0,0,0.0273,0.9727),
                            nrow = 5,
                            ncol = 5,
                            byrow = TRUE)
  
  ttau.z.mapping <- as.matrix(expand.grid(ttau = ttau.levels,z = z.levels))
  
  ll <- nrow(ttau.z.mapping)
  
  ttau.z.transition.matrix <- matrix(NA,
                                     nrow = ll,
                                     ncol = ll)
  for (i in 1:ll) {
    for (j in 1:ll) {
      
      ttau.z.transition.matrix[i,j] <- ttau.probabilities[ttau.z.mapping[i,1],ttau.z.mapping[j,1]] * z.probabilities[ttau.z.mapping[i,2],ttau.z.mapping[j,2]]
      
    }
  }
  
  return(list(ttau.z.mapping,ttau.z.transition.matrix))
}

ttau.z.transition.mapping <- make.ttau.z.transition.matrix()[[1]]
ttau.z.transition.matrix <- make.ttau.z.transition.matrix()[[2]]

set.params <- function(){
  # Parameters for the function
  aalpha <- 0.33
  bbeta <- 0.97
  ddelta <- 0.1
  pphi <- .05
  ggamma <- .2
  
  return(c(aalpha,
           bbeta,
           ddelta,
           pphi,
           ggamma))
}

params <- set.params()

k.grid <- seq(.7 * k.ss, 1.3 * k.ss, length.out = 25)
il.grid <- seq(.5 * i.ss, 1.5 * i.ss, length.out = 5)



# Define bellman operator

bellman <- function(opt.vars, state.vars, params, V, ttau.z.transition.mapping, ttau.z.transition.matrix) {
  
  ttau.list <- c(0.2,0.25,0.3)
  z.list <- c(-0.0673,-0.0336,0,0.0336,0.0673)
  
  # Unpack parameters
  aalpha <- params[1]
  bbeta <- params[2]
  ddelta <- params[3]
  pphi <- params[4]
  ggamma <- params[5]
  
  # Unpack state variables
  k <- state.vars[1]
  il <- state.vars[2]
  ttau <- ttau.list[ ttau.z.transition.mapping[state.vars[3], 1]]
  z <- z.list[ttau.z.transition.mapping[state.vars[3], 2]]
  
  # Unpack "opt.vars"
  l <- opt.vars[1]
  i <- opt.vars[2]
  
  # Calculate control variables
  r <- rate(aalpha, z, k, l)
  w <- wage(aalpha, z, k, l)
  c <- consumption(ttau, w, r, k, l, i)
  g <- govt(ttau, w, l)
  ac <- adjustment.costs(pphi, il, i)
  k1 <- k.next(ddelta, k, i, ac)
  
    # Keep k1 and i within grid range to avoid NaN in interpolation
  k1 <- max(min(k1, max(k.grid)), min(k.grid))
  i <- max(min(i, max(il.grid)), min(il.grid))
  
  # Calculate expected future value with interpolation
  V.next <- vector(length = 15L)
  for (state in 1:15) {
    V.next[state] <- interp2(x = il.grid, y = k.grid, V[,,state], xp = i, yp = k1)
  }
  
  # Calculate EV
  EV <- V.next %*% ttau.z.transition.matrix[state.vars[3], ]
  
  # Compute the value function
  V <- log(c) + ggamma * log(g) - (1 / 2) * l^2 + EV
  
  return(-V)
}


V0 <- array(0, dim = c(25,5,15))

V <- V0

V1 <- array(0, dim = c(25,5,15))

for (es1 in 1:25) {
  for (es2 in 1:5) {
    for (xs1 in 1:15) {
      
      state.vars <- c(k.grid[es1],
                      il.grid[es2],
                      xs1)
      
      opt_result <- optim(
                    par = c(l.ss, i.ss),
                    fn = bellman,
                    state.vars = state.vars,
                    params = params,
                    V = V,
                    ttau.z.transition.mapping = ttau.z.transition.mapping,
                    ttau.z.transition.matrix = ttau.z.transition.matrix,
                    method = "L-BFGS-B",
                    lower = c(0, 0.1887919),  # Lower bounds for l and i
                    upper = c(Inf, Inf)  # Upper bounds for l and i
                    )
      
       V1[es1,es2,xs1] <- -opt_result$value
       
       # After optim call
l_opt <- opt_result$par[1]
i_opt <- opt_result$par[2]
cat("Optimized l:", l_opt, ", i:", i_opt, "\n")
      
      
    }
    
  }
  
}

max(abs(V1 - V))

V <- V1

end_all <- Sys.time()

end_all - start_all

