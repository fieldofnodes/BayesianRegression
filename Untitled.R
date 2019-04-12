
# SMSTC - Gibbs sampler for rats example.
# Name of function and input parameter nt = number of iterations
ratsGibbs <- function(nt) {
  
  # Read in the data
  data <- matrix(c(
    151, 199, 246, 283, 320,
    145, 199, 249, 293, 354,
    147, 214, 263, 312, 328,
    .......................
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
  # For alpha (mean, variance)
  alphaprior <- c(0,10^6)
  # For beta (mean, variance)
  betaprior <- c(0,10^6)
  # For sigma^2 (shape, rate)
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
  # Define vector for storing simulated values
  # (rows = iterations; columns = parameter)
  output <- array(0,dim=c(nt,3))
  
  # Output the time to the screen at beginning of simulations
  cat("Start time", format(Sys.time(),"%X"), "\n")
  # Iterations:
  for (t in 1:nt) {
    # Update each parameter in turn using posterior conditional distributions:
    # Update alpha (from Normal)
    
    meanprop <- (N*T*ybar*alphaprior[2] + alphaprior[1]*param[3])/
      (N*T*alphaprior[2]+param[3])
    varprop <- param[3]*alphaprior[2]/(N*T*alphaprior[2]+param[3])
    param[1] <- rnorm(1,meanprop,sqrt(varprop))
    # Update beta (from Normal)
    meanprop <- (sumxy*betaprior[2] + betaprior[1]*param[3])/
      (N*sumxx*betaprior[2]+param[3])
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
  
  
  
  
  
  
  
  
  
}