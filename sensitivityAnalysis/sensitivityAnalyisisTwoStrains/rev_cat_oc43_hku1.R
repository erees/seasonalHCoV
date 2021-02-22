#' Reverse catalytic model 
#' With FOI changing by age
#' SA looking at using half the data HKU1 and OC43

library(tidyverse)
library(rjags)
library(binom)
library(varhandle)
require(MCMCvis)
library(cowplot)

#Read in the data

datComb <- readRDS("cleanedCode/cleanedData.RDS")

# Select only the data we want
datChan <- subset(datComb, author == "chan")
datZhou_hku1 <- subset(datComb, author == "datZhou_hku1")
datZhou_oc43 <- subset(datComb, author == "datZhou_oc43")
datMonto <- subset(datComb, author == "monto66")
datSar <- subset(datComb, author == "sar75")


datComb_oc43_hku1 <- datChan %>% 
  bind_rows(datZhou_hku1) %>% 
  bind_rows(datZhou_oc43) %>% 
  bind_rows(datMonto) %>% 
  bind_rows(datSar) 


###################
### define model
###################


jcode <- "model{
for (i in 1:7){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data
seropos_est[i] <- ifelse(age[i] < cutoff,

(lambdaChan / (lambdaChan + delta)) * (1 - exp(-(lambdaChan + delta)*age[i])),

((lambdaChan / (lambdaChan+delta))* (1-exp(-(lambdaChan +delta)*cutoff)) 
- ((lambdaChan*alpha) / ((lambdaChan*alpha)+delta) )) * exp(-((lambdaChan*alpha) + delta)*(age[i]-cutoff)) 
+ ((lambdaChan*alpha) / ((lambdaChan*alpha)+delta)))
}
for (i in 8:14){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data
seropos_est[i] <- ifelse(age[i] < cutoff,

(lambdaZhou_hk / (lambdaZhou_hk + delta)) * (1 - exp(-(lambdaZhou_hk + delta)*age[i])),

((lambdaZhou_hk / (lambdaZhou_hk+delta))* (1-exp(-(lambdaZhou_hk +delta)*cutoff)) 
- ((lambdaZhou_hk*alpha) / ((lambdaZhou_hk*alpha)+delta) )) * exp(-((lambdaZhou_hk*alpha) + delta)*(age[i]-cutoff)) 
+ ((lambdaZhou_hk*alpha) / ((lambdaZhou_hk*alpha)+delta)))
}
for (i in 15:21){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data

seropos_est[i] <- ifelse(age[i] < cutoff,

(lambdaZhou_oc / (lambdaZhou_oc + delta)) * (1 - exp(-(lambdaZhou_oc + delta)*age[i])),

((lambdaZhou_oc / (lambdaZhou_oc+delta))* (1-exp(-(lambdaZhou_oc +delta)*cutoff)) 
- ((lambdaZhou_oc*alpha) / ((lambdaZhou_oc*alpha)+delta) )) * exp(-((lambdaZhou_oc*alpha) + delta)*(age[i]-cutoff)) 
+ ((lambdaZhou_oc*alpha) / ((lambdaZhou_oc*alpha)+delta)))
}
for (i in 22:27){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data

seropos_est[i] <- ifelse(age[i] < cutoff,

(lambdaMonto / (lambdaMonto + delta)) * (1 - exp(-(lambdaMonto + delta)*age[i])),

((lambdaMonto / (lambdaMonto+delta))* (1-exp(-(lambdaMonto +delta)*cutoff)) 
- ((lambdaMonto*alpha) / ((lambdaMonto*alpha)+delta) )) * exp(-((lambdaMonto*alpha) + delta)*(age[i]-cutoff)) 
+ ((lambdaMonto*alpha) / ((lambdaMonto*alpha)+delta)))
}
for (i in 28:31){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data

seropos_est[i] <- ifelse(age[i] < cutoff,

(lambdaSar / (lambdaSar + delta)) * (1 - exp(-(lambdaSar + delta)*age[i])),

((lambdaSar / (lambdaSar+delta))* (1-exp(-(lambdaSar +delta)*cutoff)) 
- ((lambdaSar*alpha) / ((lambdaSar*alpha)+delta) )) * exp(-((lambdaSar*alpha) + delta)*(age[i]-cutoff)) 
+ ((lambdaSar*alpha) / ((lambdaSar*alpha)+delta)))
}

# priors

lambdaChan ~ dgamma(1.2,4)
lambdaZhou_hk ~ dgamma(1.2,4)
lambdaZhou_oc ~ dgamma(1.2,4)
lambdaMonto ~ dgamma(1.2,4)
lambdaSar ~ dgamma(1.2,4)

delta ~ dunif(0,5) #uniformative prior
alpha ~ dgamma(5,5)
cutoff ~ dunif(0,20)
}"

paramVector = c("lambdaZhou_oc", "lambdaMonto", "lambdaSar", "lambdaChan", "lambdaZhou_hk", "delta","alpha","cutoff")
mcmc.length=40000
jdat <- list(n.pos= datComb_oc43_hku1$N_positive,
             N=datComb_oc43_hku1$N,
             age=datComb_oc43_hku1$midpoint)

jmod = jags.model(textConnection(jcode), data=jdat, n.chains=4, n.adapt=10000)
update(jmod)

jpos = coda.samples(jmod, paramVector, n.iter=mcmc.length)

plot(jpos)


mcmcMatrix <- as.matrix(jpos)

summary(jpos)
MCMCsummary(jpos, round = 2)
mcmcDF = as_tibble(mcmcMatrix)

mcmcDF %>%
  gather() %>% 
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_density()

deltaPointEst <- mcmcMatrix[,"delta"] %>% quantile(probs=c(.5,.025,.975))
cutoffPointEst <- mcmcMatrix[,"cutoff"] %>% quantile(probs=c(.5,.025,.975))
alphaPointEst <- mcmcMatrix[,"alpha"] %>% quantile(probs=c(.5,.025,.975))
lambdaPointEstChan <- mcmcMatrix[,"lambdaChan"] %>% quantile(probs=c(.5,.025,.975))
lambdaPointEstZhouhk <- mcmcMatrix[,"lambdaZhou_hk"] %>% quantile(probs=c(.5,.025,.975))
lambdaPointEstZhouoc <- mcmcMatrix[,"lambdaZhou_oc"] %>% quantile(probs=c(.5,.025,.975))
lambdaPointEstMonto <- mcmcMatrix[,"lambdaMonto"] %>% quantile(probs=c(.5,.025,.975))
lambdaPointEstSar <- mcmcMatrix[,"lambdaSar"] %>% quantile(probs=c(.5,.025,.975))


paramEstimates = list(lambdaPointEstZhouoc,
                      lambdaPointEstMonto,
                      lambdaPointEstSar,
                      lambdaPointEstChan,
                      lambdaPointEstZhouhk,
                      deltaPointEst,
                      alphaPointEst,
                      cutoffPointEst)

## Outputting point estimates for inclusion within tables
varOutput = c()
for(i in 1:length(paramEstimates)){
  var = paramEstimates[[i]]
  varOut = paste(round(var[[1]],2)," (",round(var[[2]],2)," - ",round(var[[3]],2),")",sep = "")
  varOutput = c(varOutput,varOut)
}

paramDat = data.frame(paramVector,varOutput)
print(1/deltaPointEst) %>% round(2)
