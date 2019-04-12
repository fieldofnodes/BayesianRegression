# SMSTC - Gibbs sampler for rats example.

# Name of function and input parameter nt = number of iterations

ratsGibbs <- function(nt) {
  
  # Read in the data
  
  data <- matrix(c(
    151, 199, 246, 283, 320,
    145, 199, 249, 293, 354,
    147, 214, 263, 312, 328,
    155, 200, 237, 272, 297,
    135, 188, 230, 280, 323,
    159, 210, 252, 298, 331,
    141, 189, 231, 275, 305,
    159, 201, 248, 297, 338,
    177, 236, 285, 350, 376,
    134, 182, 220, 260, 296,
    160, 208, 261, 313, 352,
    143, 188, 220, 273, 314,
    154, 200, 244, 289, 325,
    171, 221, 270, 326, 358,
    163, 216, 242, 281, 312,
    160, 207, 248, 288, 324,
    142, 187, 234, 280, 316,
    156, 203, 243, 283, 317,
    157, 212, 259, 307, 336,
    152, 203, 246, 286, 321,
    154, 205, 253, 298, 334,
    139, 190, 225, 267, 302,
    146, 191, 229, 272, 302,
    157, 211, 250, 285, 323,
    132, 185, 237, 286, 331,
    160, 207, 257, 303, 345,
    169, 216, 261, 295, 333,
    157, 205, 248, 289, 316,
    137, 180, 219, 258, 291,
    153, 200, 244, 286, 324),nrow=30,byrow=T)
  
  x <- c(8,15,22,29,36)
  
  xbar <- 22
  
  N = 30
  T = 5
  
  # Calculate mean of response data
  
  ybar <- mean(data)
  
  # Calculate other statistics used in the posterior conditionals:
  
  sumxy <- 0
  
  for (i in 1:N) {
    for (j in 1:T) {
      sumxy <- sumxy + (x[j]-xbar)*data[i,j]
    }}
  
  sumxx <- sum((x-xbar)^2)
  
  # Read in the priors:
  
  # For alpha:
  
  alphaprior <- c(0,10^6)
  
  # For beta:
  
  betaprior <- c(0,10^6)
  
  # For sigma^2
  
  varprior <- c(0.001,0.001)
  
  # Parameters for MH updates (Uniform random walk) for alpha, beta and sigma^2:
  
  delta <- c(5,0.5,50)
  
  # Set initial parameter values:
  
  alpha <- 250
  beta <- 6
  sigma2 <- 250
  
  # Put these into vector of parameters
  
  param <- c(alpha,beta,sigma2)
  
  # param[1] = alpha; param[2] = beta; param[3] = sigma^2
  
  # Define vector for storing simulated values (rows = iterations; columns = parameter)
  
  output <- array(0,dim=c(nt,3))
  
  # Output the time to the screen at beginning of simulations
  
  cat("Start time", format(Sys.time(),"%X"), "\n")
  
  # Iterations:
  
  for (t in 1:nt) {
    
    # Update each parameter in turn using posterior conditional distributions:
    
    # Update alpha (from Normal)
    
    meanprop <- (N*T*ybar*alphaprior[2] + alphaprior[1]*param[3])/(N*T*alphaprior[2]+param[3])
    varprop <- param[3]*alphaprior[2]/(N*T*alphaprior[2]+param[3])
    
    param[1] <- rnorm(1,meanprop,sqrt(varprop))
    
    # Update beta (from Normal)
    
    meanprop <- (sumxy*betaprior[2] + betaprior[1]*param[3])/(N*sumxx*betaprior[2]+param[3])
    varprop <- param[3]*betaprior[2]/(N*sumxx*betaprior[2]+param[3])
    
    param[2] <- rnorm(1,meanprop,sqrt(varprop))
    
    # Update sigma^2 from Inverse Gamma (simulate from Gamma and take the inverse)
    
    sumgamma <- 0
    
    for (i in 1:N) {
      for (j in 1:T) {
        sumgamma <- sumgamma + (data[i,j]-param[1]-param[2]*(x[j]-xbar))^2	
        
      }}
    
    gamma1 <- N*T/2 + varprior[1]
    gamma2 <- sumgamma/2 + varprior[2]
    
    temp <- rgamma(1,gamma1,gamma2)
    param[3] <- 1/temp
    
    
    # Store parameter values in "output"
    
    output[t, ] <- param[ ]
    
  }
  
  # Output the time to the screen at end of simulations
  
  cat("Finish time", format(Sys.time(),"%X"), "\n")
  
  # Output the set of parameter values at each iteration stored in "output"
  
  output
  
}

