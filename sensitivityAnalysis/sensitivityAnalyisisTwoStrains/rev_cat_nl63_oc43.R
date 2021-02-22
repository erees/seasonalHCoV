#' Reverse catalytic model 
#' With FOI changing by age
#' SA looking at using half the data NL63 and OC43

library(tidyverse)
library(rjags)
library(binom)
library(varhandle)
require(MCMCvis)
library(cowplot)

#Read in the data

datComb <- readRDS("cleanedCode/cleanedData.RDS")

datZhou_oc43 <- subset(datComb, author == "datZhou_oc43")
datMonto <- subset(datComb, author == "monto66")
datSar <- subset(datComb, author == "sar75")
datZhou_nl63 <- subset(datComb, author == "datZhou_nl63")
datShao_nl63 <- subset(datComb, author == "datShao_nl63")

datComb_Nl63_OC43 <- datZhou_oc43 %>% 
  bind_rows(datMonto) %>% 
  bind_rows(datSar) %>% 
  bind_rows(datZhou_nl63) %>% 
  bind_rows(datShao_nl63)


###################
### define model
###################


jcode <- "model{
for (i in 1:7){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data

seropos_est[i] <- ifelse(age[i] < cutoff,

(lambdaZhou_oc / (lambdaZhou_oc + delta)) * (1 - exp(-(lambdaZhou_oc + delta)*age[i])),

((lambdaZhou_oc / (lambdaZhou_oc+delta))* (1-exp(-(lambdaZhou_oc +delta)*cutoff)) 
- ((lambdaZhou_oc*alpha) / ((lambdaZhou_oc*alpha)+delta) )) * exp(-((lambdaZhou_oc*alpha) + delta)*(age[i]-cutoff)) 
+ ((lambdaZhou_oc*alpha) / ((lambdaZhou_oc*alpha)+delta)))
}
for (i in 8:13){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data

seropos_est[i] <- ifelse(age[i] < cutoff,

(lambdaMonto / (lambdaMonto + delta)) * (1 - exp(-(lambdaMonto + delta)*age[i])),

((lambdaMonto / (lambdaMonto+delta))* (1-exp(-(lambdaMonto +delta)*cutoff)) 
- ((lambdaMonto*alpha) / ((lambdaMonto*alpha)+delta) )) * exp(-((lambdaMonto*alpha) + delta)*(age[i]-cutoff)) 
+ ((lambdaMonto*alpha) / ((lambdaMonto*alpha)+delta)))
}
for (i in 14:17){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data

seropos_est[i] <- ifelse(age[i] < cutoff,

(lambdaSar / (lambdaSar + delta)) * (1 - exp(-(lambdaSar + delta)*age[i])),

((lambdaSar / (lambdaSar+delta))* (1-exp(-(lambdaSar +delta)*cutoff)) 
- ((lambdaSar*alpha) / ((lambdaSar*alpha)+delta) )) * exp(-((lambdaSar*alpha) + delta)*(age[i]-cutoff)) 
+ ((lambdaSar*alpha) / ((lambdaSar*alpha)+delta)))
}
for (i in 18:24){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data

seropos_est[i] <- ifelse(age[i] < cutoff,

(lambdaZhou_nl / (lambdaZhou_nl + delta)) * (1 - exp(-(lambdaZhou_nl + delta)*age[i])),

((lambdaZhou_nl / (lambdaZhou_nl+delta))* (1-exp(-(lambdaZhou_nl +delta)*cutoff)) 
- ((lambdaZhou_nl*alpha) / ((lambdaZhou_nl*alpha)+delta) )) * exp(-((lambdaZhou_nl*alpha) + delta)*(age[i]-cutoff)) 
+ ((lambdaZhou_nl*alpha) / ((lambdaZhou_nl*alpha)+delta)))
}
for (i in 25:30){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data

seropos_est[i] <- ifelse(age[i] < cutoff,

(lambdaShao_nl / (lambdaShao_nl + delta)) * (1 - exp(-(lambdaShao_nl + delta)*age[i])),

((lambdaShao_nl / (lambdaShao_nl+delta))* (1-exp(-(lambdaShao_nl +delta)*cutoff)) 
- ((lambdaShao_nl*alpha) / ((lambdaShao_nl*alpha)+delta) )) * exp(-((lambdaShao_nl*alpha) + delta)*(age[i]-cutoff)) 
+ ((lambdaShao_nl*alpha) / ((lambdaShao_nl*alpha)+delta)))
}


# priors
lambdaZhou_oc ~ dgamma(1.2,4)
lambdaMonto ~ dgamma(1.2,4)
lambdaSar ~ dgamma(1.2,4)
lambdaZhou_nl ~ dgamma(1.2,4)
lambdaShao_nl ~ dgamma(1.2,4)

delta ~ dunif(0,5) #uniformative prior
alpha ~ dgamma(5,5)
cutoff ~ dunif(0,20)
}"

paramVector = c("lambdaZhou_oc", "lambdaMonto", "lambdaSar", "lambdaZhou_nl", "lambdaShao_nl", "delta","alpha","cutoff")

mcmc.length=40000
jdat <- list(n.pos= datComb_Nl63_OC43$N_positive,
             N=datComb_Nl63_OC43$N,
             age=datComb_Nl63_OC43$midpoint)

jmod = jags.model(textConnection(jcode), data=jdat, n.chains=4, n.adapt=10000)
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
lambdaPointEstZhouoc <- mcmcMatrix[,"lambdaZhou_oc"] %>% quantile(probs=c(.5,.025,.975))
lambdaPointEstMonto <- mcmcMatrix[,"lambdaMonto"] %>% quantile(probs=c(.5,.025,.975))
lambdaPointEstSar <- mcmcMatrix[,"lambdaSar"] %>% quantile(probs=c(.5,.025,.975))
lambdaPointEstZhounl <- mcmcMatrix[,"lambdaZhou_nl"] %>% quantile(probs=c(.5,.025,.975))
lambdaPointEstShaonl <- mcmcMatrix[,"lambdaShao_nl"] %>% quantile(probs=c(.5,.025,.975))

paramEstimates = list(lambdaPointEstZhouoc,
                      lambdaPointEstMonto,
                      lambdaPointEstSar,
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
