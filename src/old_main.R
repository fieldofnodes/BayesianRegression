#main file
library(ggplot2)
library(dplyr)
library(tidyr)
library(runjags)
library(coda)
library(rjags)
require(rjags)
require(coda)
require(invgamma)
library(ggmcmc)



#Include the data file
source(file = "data/data.R")

y <- data[,1]
x1 <- data[,2]
x2 <- data[,3]
n <- nrow(data)

model.inits <- list(invSigma2=1, alpha=1, beta=0, gamma=0)
iterations <- 10000
burnin <- floor(iterations/2)
chains <- 2

model.fit <- jags.model(file="src/Regression2JAGS.txt", data=list(n=n, y=y, x1=x1, x2=x2), inits=model.inits, n.chains = chains)

model.samples <- coda.samples(model.fit, c("alpha", "beta", "gamma", "R2B", "sigma2"), n.iter=iterations)

summary(window(model.samples, start = burnin))

plot(model.samples, trace=FALSE, density = TRUE)   



S <- ggs(model.samples)
str(S)
ggmcmc(S, file="fig/model_simple-diag.pdf", param_page=2)


