#' Reverse catalytic model 
#' With FOI changing by age
#' SA looking at using half the data NL63 and 229e

library(tidyverse)
library(rjags)
library(binom)
library(varhandle)
require(MCMCvis)
library(cowplot)

#Read in the data

datComb <- readRDS("cleanedCode/cleanedData.RDS")

# Select only the data we want
datShao_229e <- subset(datComb, author == "datShao_229e")
datZhou_229e <- subset(datComb, author == "datZhou_229e")
datCav <- subset(datComb, author == "cav")
datZhou_nl63 <- subset(datComb, author == "datZhou_nl63")
datShao_nl63 <- subset(datComb, author == "datShao_nl63")

datComb_Nl63_229e <- datShao_229e %>% 
  bind_rows(datZhou_229e) %>% 
  bind_rows(datCav) %>% 
  bind_rows(datZhou_nl63) %>% 
  bind_rows(datShao_nl63) 


###################
### define model
###################


jcode <- "model{
for (i in 1:6){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data
seropos_est[i] <- ifelse(age[i] < cutoff,

(lambdaShao_22 / (lambdaShao_22 + delta)) * (1 - exp(-(lambdaShao_22 + delta)*age[i])),

((lambdaShao_22 / (lambdaShao_22+delta))* (1-exp(-(lambdaShao_22 +delta)*cutoff)) 
- ((lambdaShao_22*alpha) / ((lambdaShao_22*alpha)+delta) )) * exp(-((lambdaShao_22*alpha) + delta)*(age[i]-cutoff)) 
+ ((lambdaShao_22*alpha) / ((lambdaShao_22*alpha)+delta)))
}
for (i in 7:13){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data
seropos_est[i] <- ifelse(age[i] < cutoff,

(lambdaZhou_22 / (lambdaZhou_22 + delta)) * (1 - exp(-(lambdaZhou_22 + delta)*age[i])),

((lambdaZhou_22 / (lambdaZhou_22+delta))* (1-exp(-(lambdaZhou_22 +delta)*cutoff)) 
- ((lambdaZhou_22*alpha) / ((lambdaZhou_22*alpha)+delta) )) * exp(-((lambdaZhou_22*alpha) + delta)*(age[i]-cutoff)) 
+ ((lambdaZhou_22*alpha) / ((lambdaZhou_22*alpha)+delta)))
}
for (i in 14:19){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data
seropos_est[i] <- ifelse(age[i] < cutoff,

(lambdaCav / (lambdaCav + delta)) * (1 - exp(-(lambdaCav + delta)*age[i])),

((lambdaCav / (lambdaCav+delta))* (1-exp(-(lambdaCav +delta)*cutoff)) 
- ((lambdaCav*alpha) / ((lambdaCav*alpha)+delta) )) * exp(-((lambdaCav*alpha) + delta)*(age[i]-cutoff)) 
+ ((lambdaCav*alpha) / ((lambdaCav*alpha)+delta)))
}
for (i in 20:26){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data

seropos_est[i] <- ifelse(age[i] < cutoff,

(lambdaZhou_nl / (lambdaZhou_nl + delta)) * (1 - exp(-(lambdaZhou_nl + delta)*age[i])),

((lambdaZhou_nl / (lambdaZhou_nl+delta))* (1-exp(-(lambdaZhou_nl +delta)*cutoff)) 
- ((lambdaZhou_nl*alpha) / ((lambdaZhou_nl*alpha)+delta) )) * exp(-((lambdaZhou_nl*alpha) + delta)*(age[i]-cutoff)) 
+ ((lambdaZhou_nl*alpha) / ((lambdaZhou_nl*alpha)+delta)))
}
for (i in 27:32){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data

seropos_est[i] <- ifelse(age[i] < cutoff,

(lambdaShao_nl / (lambdaShao_nl + delta)) * (1 - exp(-(lambdaShao_nl + delta)*age[i])),

((lambdaShao_nl / (lambdaShao_nl+delta))* (1-exp(-(lambdaShao_nl +delta)*cutoff)) 
- ((lambdaShao_nl*alpha) / ((lambdaShao_nl*alpha)+delta) )) * exp(-((lambdaShao_nl*alpha) + delta)*(age[i]-cutoff)) 
+ ((lambdaShao_nl*alpha) / ((lambdaShao_nl*alpha)+delta)))
}


# priors

lambdaShao_22 ~ dgamma(1.2,4)
lambdaZhou_22 ~ dgamma(1.2,4)
lambdaCav ~ dgamma(1.2,4)
lambdaZhou_nl ~ dgamma(1.2,4)
lambdaShao_nl ~ dgamma(1.2,4)

delta ~ dunif(0,5) #uniformative prior
alpha ~ dgamma(5,5)
cutoff ~ dunif(0,20)
}"

paramVector = c("lambdaShao_22", "lambdaZhou_22", "lambdaCav", "lambdaZhou_nl", "lambdaShao_nl", "delta","alpha","cutoff")

mcmc.length=40000
jdat <- list(n.pos= datComb_Nl63_229e$N_positive,
             N=datComb_Nl63_229e$N,
             age=datComb_Nl63_229e$midpoint)

jmod = jags.model(textConnection(jcode), data=jdat, n.chains=4, n.adapt=10000)
update(jmod)

jpos = coda.samples(jmod, paramVector, n.iter=mcmc.length)
# plot(jpos[[1]]) # check convergence
plot(jpos)


mcmcMatrix <- as.matrix(jpos)
mcmcDF = as_tibble(mcmcMatrix)

mcmcDF %>%
  gather() %>% 
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_density()

summary(jpos)
MCMCsummary(jpos, round = 2)

deltaPointEst <- mcmcMatrix[,"delta"] %>% quantile(probs=c(.5,.025,.975))
alphaPointEst <- mcmcMatrix[,"alpha"] %>% quantile(probs=c(.5,.025,.975))
cutoffPointEst <- mcmcMatrix[,"cutoff"] %>% quantile(probs=c(.5,.025,.975))

lambdaPointEstShao22 <- mcmcMatrix[,"lambdaShao_22"] %>% quantile(probs=c(.5,.025,.975))
lambdaPointEstZhou22 <- mcmcMatrix[,"lambdaZhou_22"] %>% quantile(probs=c(.5,.025,.975))
lambdaPointEstCav <- mcmcMatrix[,"lambdaCav"] %>% quantile(probs=c(.5,.025,.975))
lambdaPointEstZhounl <- mcmcMatrix[,"lambdaZhou_nl"] %>% quantile(probs=c(.5,.025,.975))
lambdaPointEstShaonl <- mcmcMatrix[,"lambdaShao_nl"] %>% quantile(probs=c(.5,.025,.975))



paramEstimates = list(lambdaPointEstShao22,
                      lambdaPointEstZhou22,
                      lambdaPointEstCav,
                      lambdaPointEstZhounl,
                      lambdaPointEstShaonl,
                      deltaPointEst,
                      alphaPointEst,
                      cutoffPointEst)

varOutput = c()
for(i in 1:length(paramEstimates)){
  var = paramEstimates[[i]]
  varOut = paste(round(var[[1]],2)," (",round(var[[2]],2)," - ",round(var[[3]],2),")",sep = "")
  varOutput = c(varOutput,varOut)
}

paramDat = data.frame(paramVector,varOutput)

