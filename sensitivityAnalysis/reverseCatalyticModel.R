# Reverse catalytic model 
# FOI estimated by study
# waning jointly estimated across all studies

library(tidyverse)
library(rjags)
library(binom)
library(varhandle)
require(MCMCvis)
library(cowplot)

datComb <- readRDS("dataProcessing/cleanedData.RDS")

###################
### define model
###################

jcode <- "model{
for (i in 1:6){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data
seropos_est[i] = (lambdaShao_22 / (lambdaShao_22 + delta)) * (1-exp(-(lambdaShao_22+delta)*age[i])) 
}
for (i in 7:13){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data
seropos_est[i] = (lambdaZhou_22 / (lambdaZhou_22 + delta)) * (1-exp(-(lambdaZhou_22+delta)*age[i])) 
}
for (i in 14:19){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data
seropos_est[i] = (lambdaCav / (lambdaCav + delta)) * (1-exp(-(lambdaCav+delta)*age[i]))
}
for (i in 20:26){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data
seropos_est[i] = (lambdaChan / (lambdaChan + delta)) * (1-exp(-(lambdaChan+delta)*age[i]))
}
for (i in 27:33){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data
seropos_est[i] = (lambdaZhou_hk / (lambdaZhou_hk + delta)) * (1-exp(-(lambdaZhou_hk+delta)*age[i]))
}
for (i in 34:40){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data
seropos_est[i] = (lambdaZhou_oc / (lambdaZhou_oc + delta)) * (1-exp(-(lambdaZhou_oc+delta)*age[i]))
}
for (i in 41:46){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data
seropos_est[i] = (lambdaMonto / (lambdaMonto + delta)) * (1-exp(-(lambdaMonto+delta)*age[i]))
}
for (i in 47:50){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data
seropos_est[i] = (lambdaSar / (lambdaSar + delta)) * (1-exp(-(lambdaSar+delta)*age[i])) 
}
for (i in 51:57){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data
seropos_est[i] = (lambdaZhou_nl / (lambdaZhou_nl + delta)) * (1-exp(-(lambdaZhou_nl+delta)*age[i])) 
}
for (i in 58:63){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data
seropos_est[i] = (lambdaShao_nl / (lambdaShao_nl + delta)) * (1-exp(-(lambdaShao_nl+delta)*age[i])) 
}


# priors

lambdaShao_22 ~ dnorm(0.3,0.5)
lambdaZhou_22 ~ dnorm(0.3,0.5)
lambdaCav ~ dnorm(0.3,0.5)
lambdaChan ~ dnorm(0.3,0.5)
lambdaZhou_hk ~ dnorm(0.3,0.5)
lambdaZhou_oc ~ dnorm(0.3,0.5)
lambdaMonto ~ dnorm(0.3,0.5)
lambdaSar ~ dnorm(0.3,0.5)
lambdaZhou_nl ~ dnorm(0.3,0.5)
lambdaShao_nl ~ dnorm(0.3,0.5)
delta ~ dunif(0,5)

}"

paramVector = c("lambdaShao_22", "lambdaZhou_22", "lambdaCav", "lambdaChan", "lambdaZhou_hk", "lambdaZhou_oc", "lambdaMonto", "lambdaSar", "lambdaZhou_nl", "lambdaShao_nl","delta")


mcmc.length=50000
jdat <- list(n.pos= datComb$N_positive,
             N=datComb$N,
             age=datComb$midpoint)

jmod = jags.model(textConnection(jcode), data=jdat, n.chains=4,n.adapt = 15000)
jpos = coda.samples(jmod, paramVector, n.iter=mcmc.length)

plot(jpos) # check convergence

summary(jpos)
MCMCsummary(jpos, round = 2)

dic.samples(jmod, n.iter = mcmc.length)

mcmcMatrix <- as.matrix(jpos)
mcmcDF = as_tibble(mcmcMatrix)

# Plotting posterior distributions of all parameters
mcmcDF %>%
  gather() %>% 
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_density()

################################################################################
## Create point estimates for all parameters
################################################################################

deltaPointEst <- mcmcMatrix[,"delta"] %>% quantile(probs=c(.5,.025,.975))
lambdaPointEstShao22 <- mcmcMatrix[,"lambdaShao_22"] %>% quantile(probs=c(.5,.025,.975))
lambdaPointEstZhou22 <- mcmcMatrix[,"lambdaZhou_22"] %>% quantile(probs=c(.5,.025,.975))
lambdaPointEstCav <- mcmcMatrix[,"lambdaCav"] %>% quantile(probs=c(.5,.025,.975))
lambdaPointEstChan <- mcmcMatrix[,"lambdaChan"] %>% quantile(probs=c(.5,.025,.975))
lambdaPointEstZhouhk <- mcmcMatrix[,"lambdaZhou_hk"] %>% quantile(probs=c(.5,.025,.975))
lambdaPointEstZhouoc <- mcmcMatrix[,"lambdaZhou_oc"] %>% quantile(probs=c(.5,.025,.975))
lambdaPointEstMonto <- mcmcMatrix[,"lambdaMonto"] %>% quantile(probs=c(.5,.025,.975))
lambdaPointEstSar <- mcmcMatrix[,"lambdaSar"] %>% quantile(probs=c(.5,.025,.975))
lambdaPointEstZhounl <- mcmcMatrix[,"lambdaZhou_nl"] %>% quantile(probs=c(.5,.025,.975))
lambdaPointEstShaonl <- mcmcMatrix[,"lambdaShao_nl"] %>% quantile(probs=c(.5,.025,.975))

paramEstimates = list(lambdaPointEstShao22,
                      lambdaPointEstZhou22,
                      lambdaPointEstCav,
                      lambdaPointEstChan,
                      lambdaPointEstZhouhk,
                      lambdaPointEstZhouoc,
                      lambdaPointEstMonto,
                      lambdaPointEstSar,
                      lambdaPointEstZhounl,
                      lambdaPointEstShaonl,
                      deltaPointEst)

## Outputting point estimates for inclusion within tables
varOutput = c()
for(i in 1:length(paramEstimates)){
  var = paramEstimates[[i]]
  varOut = paste(round(var[[1]],2)," (",round(var[[2]],2)," - ",round(var[[3]],2),")",sep = "")
  varOutput = c(varOutput,varOut)
}

paramDat = data.frame(paramVector,varOutput)

################################################################################
## Create data for plots
## Samples from mcmc chains for credible intervals
################################################################################

ager=0:80
numSamples = 1000

foiVector = paramVector[1:10]
foiEstimates = paramEstimates[1:10]

for(ii in 1:length(foiVector)){
  outDf <- matrix(NA,nrow=numSamples, ncol = length(ager))
  foiStudy <- foiVector[ii]
  
  for (kk in 1:numSamples ) {
    randomNumber <- floor(runif(1, min = 1, max = nrow(mcmcMatrix)))
    
    lambdaSample <- mcmcMatrix[randomNumber,foiStudy]
    deltaSample <- mcmcMatrix[randomNumber,"delta"]
    
    newRow <- (lambdaSample / (lambdaSample+deltaSample)) *(1 - exp(-ager*(lambdaSample+deltaSample)))
    outDf[kk,] <- newRow
  }
  
  # for each row in the matrix get quantiles
  quantileMatrix <- matrix(NA,nrow=ncol(outDf), ncol = 3)
  for(jj in 1:ncol(outDf)){
    quantiles <- outDf[,jj] %>% quantile(probs=c(.5,.025,.975))
    quantileMatrix[jj,] <- quantiles
  }
  
  foi = foiEstimates[[ii]][1]
  delta = deltaPointEst[1]
  
  meanEs = (foi / (foi+delta)) *(1 - exp(-ager*(foi+delta)))
  
  # Create a dataframe for plotting
  df_upperLower = data.frame(
    midpoint = ager,
    mean = meanEs,
    upper = quantileMatrix[,3],
    lower = quantileMatrix[,2]
  )
  
  assign(paste0("df_", foiVector[ii]), df_upperLower)
}

####################
## Plots
####################

hku1Plot <- ggplot(df_lambdaChan, aes(x=midpoint, y=mean, ymin=lower, ymax=upper)) +
  geom_ribbon(alpha=0.2, fill ="#558C8C")+
  geom_line( color = "#558C8C")+
  geom_point(data=datChan, color = "#558C8C")+
  geom_linerange(data=datChan, color = "#558C8C") +
  geom_line(data=df_lambdaZhou_hk, aes(x = midpoint, y=mean),color = "#C05746") +
  geom_ribbon(data = df_lambdaZhou_hk, alpha=0.2, fill = "#C05746")+
  geom_point(data=datZhou_hku1,color="#C05746")+
  geom_linerange(data=datZhou_hku1,color="#C05746") +
  xlab("Age (years)") + ylab("Proportion Seropositive") +
  scale_x_continuous(breaks=seq(0,100,by=10)) +
  ylim(0,1) +
  theme_bw()
hku1Plot

nl63Plot <- ggplot(df_lambdaShao_nl, aes(x=midpoint, y=mean, ymin=lower, ymax=upper)) +
  geom_ribbon(alpha=0.2, fill ="#075E9D")+
  geom_line( color = "#075E9D")+
  geom_point(data=datShao_nl63, color = "#075E9D")+
  geom_linerange(data=datShao_nl63, color = "#075E9D") +
  geom_line(data=df_lambdaZhou_nl, aes(x = midpoint, y=mean),color = "#C05746") +
  geom_ribbon(data = df_lambdaZhou_nl, alpha=0.2, fill = "#C05746")+
  geom_point(data=datZhou_nl63,color="#C05746")+
  geom_linerange(data=datZhou_nl63,color="#C05746") +
  xlab("Age (years)") + ylab("Proportion Seropositive") +
  scale_x_continuous(breaks=seq(0,100,by=10)) +
  ylim(0,1) +
  theme_bw()
nl63Plot

oc43Plot <- ggplot(df_lambdaZhou_oc, aes(x=midpoint, y=mean, ymin=lower, ymax=upper)) +
  geom_ribbon(alpha=0.2, fill ="#C05746")+
  geom_line( color = "#C05746")+
  geom_point(data=datZhou_oc43,color="#C05746")+
  geom_linerange(data=datZhou_oc43,color="#C05746") +
  geom_line(data=df_lambdaMonto, aes(x = midpoint, y=mean),color = "#2C484E") + #grey
  geom_ribbon(data = df_lambdaMonto, alpha=0.2, fill = "#2C484E")+
  geom_point(data=datMonto,color="#2C484E")+
  geom_linerange(data=datMonto,color="#2C484E") +
  geom_line(data=df_lambdaSar, aes(x = midpoint, y=mean),color = "#066025") + # dark green
  geom_ribbon(data = df_lambdaSar, alpha=0.2, fill = "#066025")+
  geom_point(data=datSar,color="#066025")+
  geom_linerange(data=datSar,color="#066025") +
  xlab("Age (years)") + ylab("Proportion Seropositive") +
  scale_x_continuous(breaks=seq(0,100,by=10)) +
  ylim(0,1) +
  theme_bw()
oc43Plot

plot229e <- ggplot(df_lambdaShao_22, aes(x=midpoint, y=mean, ymin=lower, ymax=upper)) +
  geom_ribbon(alpha=0.2, fill ="#075E9D")+
  geom_line( color = "#075E9D")+
  geom_point(data=datShao_229e,color="#075E9D")+
  geom_linerange(data=datShao_229e,color="#075E9D") +
  geom_line(data=df_lambdaZhou_22, aes(x = midpoint, y=mean),color = "#C05746") +
  geom_ribbon(data = df_lambdaZhou_22, alpha=0.2, fill = "#C05746")+
  geom_point(data=datZhou_229e,color="#C05746")+
  geom_linerange(data=datZhou_229e,color="#C05746") +
  geom_line(data=df_lambdaCav, aes(x = midpoint, y=mean),color = "#C9A92C") +
  geom_ribbon(data = df_lambdaCav, alpha=0.2, fill = "#C9A92C")+
  geom_point(data=datCav,color="#C9A92C")+
  geom_linerange(data=datCav,color="#C9A92C") +
  xlab("Age (years)") + ylab("Proportion Seropositive") +
  scale_x_continuous(breaks=seq(0,100,by=10)) +
  ylim(0,1) +
  theme_bw()

plot229e

plot_grid(plot229e,hku1Plot,nl63Plot,oc43Plot)

