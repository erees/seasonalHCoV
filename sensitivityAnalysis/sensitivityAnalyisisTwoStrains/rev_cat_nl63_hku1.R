#' Reverse catalytic model 
#' With FOI changing by age
#' SA looking at using half the data NL63 and HKU1

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
datZhou_nl63 <- subset(datComb, author == "datZhou_nl63")
datShao_nl63 <- subset(datComb, author == "datShao_nl63")

datComb_Nl63_hku1 <- datChan %>% 
  bind_rows(datZhou_hku1) %>% 
  bind_rows(datZhou_nl63) %>% 
  bind_rows(datShao_nl63) 


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

(lambdaZhou_nl / (lambdaZhou_nl + delta)) * (1 - exp(-(lambdaZhou_nl + delta)*age[i])),

((lambdaZhou_nl / (lambdaZhou_nl+delta))* (1-exp(-(lambdaZhou_nl +delta)*cutoff)) 
- ((lambdaZhou_nl*alpha) / ((lambdaZhou_nl*alpha)+delta) )) * exp(-((lambdaZhou_nl*alpha) + delta)*(age[i]-cutoff)) 
+ ((lambdaZhou_nl*alpha) / ((lambdaZhou_nl*alpha)+delta)))
}
for (i in 22:27){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data

seropos_est[i] <- ifelse(age[i] < cutoff,

(lambdaShao_nl / (lambdaShao_nl + delta)) * (1 - exp(-(lambdaShao_nl + delta)*age[i])),

((lambdaShao_nl / (lambdaShao_nl+delta))* (1-exp(-(lambdaShao_nl +delta)*cutoff)) 
- ((lambdaShao_nl*alpha) / ((lambdaShao_nl*alpha)+delta) )) * exp(-((lambdaShao_nl*alpha) + delta)*(age[i]-cutoff)) 
+ ((lambdaShao_nl*alpha) / ((lambdaShao_nl*alpha)+delta)))
}


# priors

lambdaChan ~ dgamma(1.2,4)
lambdaZhou_hk ~ dgamma(1.2,4)
lambdaZhou_nl ~ dgamma(1.2,4)
lambdaShao_nl ~ dgamma(1.2,4)

delta ~ dunif(0,5) #uniformative prior
alpha ~ dgamma(5,5)
cutoff ~ dunif(0,20)
}"

paramVector = c("lambdaChan", "lambdaZhou_hk", "lambdaZhou_nl", "lambdaShao_nl", "delta","alpha","cutoff")

mcmc.length=80000
jdat <- list(n.pos= datComb_Nl63_hku1$N_positive,
             N=datComb_Nl63_hku1$N,
             age=datComb_Nl63_hku1$midpoint)

jmod = jags.model(textConnection(jcode), data=jdat, n.chains=4, n.adapt=30000)
update(jmod)

jpos = coda.samples(jmod, paramVector, n.iter=mcmc.length)
# plot(jpos[[1]]) # check convergence
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

lambdaPointEstZhounl <- mcmcMatrix[,"lambdaZhou_nl"] %>% quantile(probs=c(.5,.025,.975))
lambdaPointEstShaonl <- mcmcMatrix[,"lambdaShao_nl"] %>% quantile(probs=c(.5,.025,.975))


paramEstimates = list(lambdaPointEstChan,
                      lambdaPointEstZhouhk,
                      lambdaPointEstZhounl,
                      lambdaPointEstShaonl,
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
