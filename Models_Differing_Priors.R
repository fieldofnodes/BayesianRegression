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
library(gridExtra)
library(ggpubr)




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


#Differing priors
model.fit.vague <- jags.model(file="src/Regression2JAGSVaguePrior.txt", data=list(n=n, y=y, x1=x1, x2=x2), inits=model.inits, n.chains = chains)
model.fit.flat <- jags.model(file="src/Regression2JAGSFlatPrior.txt", data=list(n=n, y=y, x1=x1, x2=x2), inits=model.inits, n.chains = chains)
model.fit.weak <- jags.model(file="src/Regression2JAGSWeakVaguePrior.txt", data=list(n=n, y=y, x1=x1, x2=x2), inits=model.inits, n.chains = chains)

#Creating samples
model.samples.vague <- coda.samples(model.fit.vague, c("alpha", "beta", "gamma", "sigma2"), n.iter=iterations)
model.samples.flat <- coda.samples(model.fit.flat, c("alpha", "beta", "gamma", "sigma2"), n.iter=iterations)
model.samples.weak <- coda.samples(model.fit.weak, c("alpha", "beta", "gamma", "sigma2"), n.iter=iterations)

#Summary
summary(window(model.samples.vague, start = burnin))
summary(window(model.samples.flat, start = burnin))
summary(window(model.samples.weak, start = burnin))

#GGS
S.vague <- ggs(model.samples.vague)
S.flat <- ggs(model.samples.flat)
S.weak <- ggs(model.samples.weak)


#PriorVar
S.vague <- mutate(S.vague, PriorVar = "Vague")
S.vague <- mutate(S.vague, PriorVar = "Flat")
S.vague <- mutate(S.vague, PriorVar = "Weak")

#Structure
str(S.vague)
str(S.flat)
str(S.weak)

#Plot to file
ggmcmc(S.vague, file="fig/model_simple-diag.vague.pdf", param_page=2)
ggmcmc(S.flat, file="fig/model_simple-diag.flat.pdf", param_page=2)
ggmcmc(S.weak, file="fig/model_simple-diag.weak.pdf", param_page=2)

S.vague.alpha <- S.vague %>% select(Iteration, Chain, Parameter, value) %>% filter(Parameter == "alpha")
S.flat.alpha <- S.flat %>% select(Iteration, Chain, Parameter, value) %>% filter(Parameter == "alpha")
S.weak.alpha <- S.weak %>% select(Iteration, Chain, Parameter, value) %>% filter(Parameter == "alpha")



a <- ggs_traceplot(S.vague.alpha)
a
b <- ggs_traceplot(S.flat.alpha)
b

ggarrange(a,b)
ggs_traceplot(S)+ggtitle("Convergence of parameters per chain")


ggs_compare_partial(S)+ggtitle("Comparison of the whole chain to the last 10%")




