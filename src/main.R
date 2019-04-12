#main file
#If you do not have these packages installed uncomment the line below and run.
#install.packages(c("ggplot2","dplyr","tidyr","runjags","coda","rjags","ggmcmc"))
library(ggplot2)
library(dplyr)
library(tidyr)
library(runjags)
library(coda)
library(rjags)
require(rjags)
require(coda)
#require(invgamma)
library(ggmcmc)



#Include the data file
#source(file = "data/data.R") for ease I have included the data here
y = c(13.9,27.3,21.0,3.7,48.4,86.4,7.1,36.6,22.7,11.4,40.7,67.0,83.7,111.9,120.8,11.4,18.2,8.6,34.6,70.4)
x = c(3,6,4,1,9,15,1,7,6,2,9,12,14,17,18,3,5,2,8,12)
x2 = x^2
data = data.frame(y,x,x2)

#We name of the variables.
y <- data[,1]
x1 <- data[,2]
x2 <- data[,3]
n <- nrow(data)

#Initialising the code with the intersect in play and the variance of the noise.
model.inits <- list(invSigma2=1, alpha=1, beta=0, gamma=0)
#Number of iterations 
iterations <- 10000
#Using a conservative burn in period, starting after half of the iterations have been performed.
burnin <- floor(iterations/2)
#Default number of chains in JAGS is 3, but using two for simplicity of view.
chains <- 2

#Fitting model from JAGS code and specifying the variable names according to the data. x2 is x^2 which has been performed in the file = "data/data.R"
model.fit <- jags.model(file="src/Regression2JAGSVaguePrior.txt", data=list(n=n, y=y, x1=x1, x2=x2), inits=model.inits, n.chains = chains)

#Coda samples coerces the data from model.fit to a list of teh MCMC flavour.
model.samples <- coda.samples(model.fit, c("alpha", "beta", "gamma", "sigma2"), n.iter=iterations)

#A basic summary of the results, including the mean of the parameters.
summary(window(model.samples, start = burnin))

#From the packagfe ggmcmc using the grammar of graphics, we have a function called ggs, which coerces the mcmc.list file of "model.samples" to a dataframe of the parameters we are concerned with.
S <- ggs(model.samples)

#To show the structure of S.
str(S)

#If you run this file you will get a pdf named model_simple-diag.pdf, in which case there will be several pages of plots for your perusal.
ggmcmc(S, file="fig/model_simple-diag.pdf", param_page=2)


#This plot will show the convergence of each group over the number of chains.
ggs_traceplot(S)+ggtitle("Convergence of parameters per chain")

#Plotting a comparision of the whole chain to the last 10% of the chain. The distributions are meant to be similar.
ggs_compare_partial(S)+ggtitle("Comparison of the whole chain to the last 10%")


#Check this and write up a report on the results of prior sensitivity.

