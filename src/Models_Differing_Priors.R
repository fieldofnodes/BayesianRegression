#This file shows the modelling of differeing priors and also a density plot of each of the
#parameters posterior distributions.
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
library(latex2exp)



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
burnin <- floor((iterations+1000)/2)
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
summary.vague <- summary(window(model.samples.vague, start = burnin))
summary.flat <- summary(window(model.samples.flat, start = burnin))
summary.weak <- summary(window(model.samples.weak, start = burnin))

########################################################################################################################
########################################################################################################################
########################################################################################################################
#Model done and now plotting


#GGS
S.vague <- ggs(model.samples.vague)
S.flat <- ggs(model.samples.flat)
S.weak <- ggs(model.samples.weak)

#Plot to file
ggmcmc(S.vague, file="fig/model_simple-diag.vague.pdf", param_page=2)
ggmcmc(S.flat, file="fig/model_simple-diag.flat.pdf", param_page=2)
ggmcmc(S.weak, file="fig/model_simple-diag.weak.pdf", param_page=2)

#Structure
str(S.vague)
str(S.flat)
str(S.weak)

#PriorVar
S.vague.NEW <- mutate(S.flat, PriorVar = "Vague")
S.flat.NEW <- mutate(S.weak, PriorVar = "Flat")
S.weak.NEW <- mutate(S.vague, PriorVar = "Weak")

#Append Flat->Weak->Vague, by having Flat, Weak and then Vague
S.total <- dplyr::bind_rows(S.flat.NEW,S.weak.NEW,S.vague.NEW)
str(S.total)


alppha <- S.total %>% select(Iteration, Chain, Parameter, value, PriorVar) %>% filter(Parameter == "alpha")
betta <- S.total %>% select(Iteration, Chain, Parameter, value, PriorVar) %>% filter(Parameter == "beta")
gammmmma <- S.total %>% select(Iteration, Chain, Parameter, value, PriorVar) %>% filter(Parameter == "gamma")
sigggma <- S.total %>% select(Iteration, Chain, Parameter, value, PriorVar) %>% filter(Parameter == "sigma2")


#GGPLOT S.total grouped by parameter and chain 
talpha<-ggplot(alppha, aes(value, colour = PriorVar))+geom_density()+ggtitle(TeX('$\\alpha$'))+facet_grid(.~Chain,scales = 'free_x')+scale_x_continuous(limits = c(-10, 20))+scale_colour_discrete(name ="Prior Sensitivity alpha",labels = c('Flat \nU(-100,100)','Vague \nN(0,0.000001)','Weak \nN(0,0.0001)'))
tbeta<-ggplot(betta, aes(value, colour = PriorVar))+geom_density()+ggtitle(TeX('$\\beta$'))+facet_grid(.~Chain,scales = 'free_x')+scale_x_continuous(limits = c(0, 6))+scale_colour_discrete(name ="Prior Sensitivity beta",labels = c('Flat \nU(-100,100)','Vague \nN(0,0.000001)','Weak \nN(0,0.0001)'))
tgamma<-ggplot(gammmmma, aes(value, colour = PriorVar))+geom_density()+ggtitle(TeX('$\\gamma$'))+facet_grid(.~Chain,scales = 'free_x')+scale_x_continuous(limits = c(0.0, 0.4))+scale_colour_discrete(name ="Prior Sensitivity gamma",labels = c('Flat \nU(-100,100)','Vague \nN(0,0.000001)','Weak \nN(0,0.0001)'))
tsigma2<-ggplot(sigggma, aes(value, colour = PriorVar))+geom_density()+ggtitle(TeX('$\\sigma^2$'))+facet_grid(.~Chain,scales = 'free')+scale_x_continuous(limits = c(0.0, 35))+scale_colour_discrete(name ="Prior Sensitivity sigma2",labels = c('Flat \nU(0,20)','Vague \nInvGFamma(0.001,0.001)','Weak \nInvGamma(1,1)'))

#Shows plot in grid and saves to file in "fig/"
plotpriors <- ggarrange(talpha,tbeta,tgamma,tsigma2,ncol = 1, nrow = 4)
ggsave("fig/DifferingPriorsSensitivity.pdf", plot = plotpriors, dpi = "retina")

