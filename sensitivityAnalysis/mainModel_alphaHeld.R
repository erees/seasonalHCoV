#' Reverse catalytic model 
#' With FOI changing by age and by study
#' Using alpha as the relative chang in FOI in old the older age group
#' Waning, cutoff and alpha jointly estimated across all studies

library(tidyverse)
library(rjags)
library(binom)
library(varhandle)
require(MCMCvis)
library(cowplot)

#Read in clean data
datComb <- readRDS("dataProcessing/cleanedData.RDS")

# Individual datasets required for plotting later
datChan <- subset(datComb, author == "chan")
datZhou_hku1 <- subset(datComb, author == "datZhou_hku1")
datZhou_oc43 <- subset(datComb, author == "datZhou_oc43")
datMonto <- subset(datComb, author == "monto66")
datSar <- subset(datComb, author == "sar75")
datShao_229e <- subset(datComb, author == "datShao_229e")
datZhou_229e <- subset(datComb, author == "datZhou_229e")
datCav <- subset(datComb, author == "cav")
datZhou_nl63 <- subset(datComb, author == "datZhou_nl63")
datShao_nl63 <- subset(datComb, author == "datShao_nl63")

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

(lambdaChan / (lambdaChan + delta)) * (1 - exp(-(lambdaChan + delta)*age[i])),

((lambdaChan / (lambdaChan+delta))* (1-exp(-(lambdaChan +delta)*cutoff)) 
- ((lambdaChan*alpha) / ((lambdaChan*alpha)+delta) )) * exp(-((lambdaChan*alpha) + delta)*(age[i]-cutoff)) 
+ ((lambdaChan*alpha) / ((lambdaChan*alpha)+delta)))
}
for (i in 27:33){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data
seropos_est[i] <- ifelse(age[i] < cutoff,

(lambdaZhou_hk / (lambdaZhou_hk + delta)) * (1 - exp(-(lambdaZhou_hk + delta)*age[i])),

((lambdaZhou_hk / (lambdaZhou_hk+delta))* (1-exp(-(lambdaZhou_hk +delta)*cutoff)) 
- ((lambdaZhou_hk*alpha) / ((lambdaZhou_hk*alpha)+delta) )) * exp(-((lambdaZhou_hk*alpha) + delta)*(age[i]-cutoff)) 
+ ((lambdaZhou_hk*alpha) / ((lambdaZhou_hk*alpha)+delta)))
}
for (i in 34:40){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data

seropos_est[i] <- ifelse(age[i] < cutoff,

(lambdaZhou_oc / (lambdaZhou_oc + delta)) * (1 - exp(-(lambdaZhou_oc + delta)*age[i])),

((lambdaZhou_oc / (lambdaZhou_oc+delta))* (1-exp(-(lambdaZhou_oc +delta)*cutoff)) 
- ((lambdaZhou_oc*alpha) / ((lambdaZhou_oc*alpha)+delta) )) * exp(-((lambdaZhou_oc*alpha) + delta)*(age[i]-cutoff)) 
+ ((lambdaZhou_oc*alpha) / ((lambdaZhou_oc*alpha)+delta)))
}
for (i in 41:46){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data

seropos_est[i] <- ifelse(age[i] < cutoff,

(lambdaMonto / (lambdaMonto + delta)) * (1 - exp(-(lambdaMonto + delta)*age[i])),

((lambdaMonto / (lambdaMonto+delta))* (1-exp(-(lambdaMonto +delta)*cutoff)) 
- ((lambdaMonto*alpha) / ((lambdaMonto*alpha)+delta) )) * exp(-((lambdaMonto*alpha) + delta)*(age[i]-cutoff)) 
+ ((lambdaMonto*alpha) / ((lambdaMonto*alpha)+delta)))
}
for (i in 47:50){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data

seropos_est[i] <- ifelse(age[i] < cutoff,

(lambdaSar / (lambdaSar + delta)) * (1 - exp(-(lambdaSar + delta)*age[i])),

((lambdaSar / (lambdaSar+delta))* (1-exp(-(lambdaSar +delta)*cutoff)) 
- ((lambdaSar*alpha) / ((lambdaSar*alpha)+delta) )) * exp(-((lambdaSar*alpha) + delta)*(age[i]-cutoff)) 
+ ((lambdaSar*alpha) / ((lambdaSar*alpha)+delta)))
}
for (i in 51:57){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data

seropos_est[i] <- ifelse(age[i] < cutoff,

(lambdaZhou_nl / (lambdaZhou_nl + delta)) * (1 - exp(-(lambdaZhou_nl + delta)*age[i])),

((lambdaZhou_nl / (lambdaZhou_nl+delta))* (1-exp(-(lambdaZhou_nl +delta)*cutoff)) 
- ((lambdaZhou_nl*alpha) / ((lambdaZhou_nl*alpha)+delta) )) * exp(-((lambdaZhou_nl*alpha) + delta)*(age[i]-cutoff)) 
+ ((lambdaZhou_nl*alpha) / ((lambdaZhou_nl*alpha)+delta)))
}
for (i in 58:63){ 
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
lambdaChan ~ dgamma(1.2,4)
lambdaZhou_hk ~ dgamma(1.2,4)
lambdaZhou_oc ~ dgamma(1.2,4)
lambdaMonto ~ dgamma(1.2,4)
lambdaSar ~ dgamma(1.2,4)
lambdaZhou_nl ~ dgamma(1.2,4)
lambdaShao_nl ~ dgamma(1.2,4)

delta ~ dunif(0,5) #uniformative prior
alpha ~ dgamma(5,5)
cutoff ~ dunif(0,20)
}"

paramVector = c("lambdaShao_22", "lambdaZhou_22", "lambdaCav", "lambdaChan", "lambdaZhou_hk", "lambdaZhou_oc", "lambdaMonto", "lambdaSar", "lambdaZhou_nl", "lambdaShao_nl","delta","alpha","cutoff")

mcmc.length=100000
jdat <- list(n.pos= datComb$N_positive,
             N=datComb$N,
             age=datComb$midpoint)

jmod = jags.model(textConnection(jcode), data=jdat, n.chains=4, n.adapt=30000)
update(jmod)

jpos = coda.samples(jmod, paramVector, n.iter=mcmc.length)

plot(jpos)
MCMCsummary(jpos, round = 2)
summary(jpos)

mcmcMatrix <- as.matrix(jpos)

mcmcDF = as_tibble(mcmcMatrix)

# Plotting posterior distributions of all parameters
mcmcDF %>%
  gather() %>% 
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_density()

## Export posteriors for alpha (used for Fig. 2)
mcmcAlpha = mcmcDF[,1]
saveRDS(mcmcAlpha, "plots/mcmcAlphaPostAll.RDS")

# Calculate DIC
dic.samples(jmod, n.iter = mcmc.length)

################################################################################
## Create point estimates for all parameters
################################################################################


deltaPointEst <- mcmcMatrix[,"delta"] %>% quantile(probs=c(.5,.025,.975))
alphaPointEst <- mcmcMatrix[,"alpha"] %>% quantile(probs=c(.5,.025,.975))
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
cutoffPointEst <- mcmcMatrix[,"cutoff"] %>% quantile(probs=c(.5,.025,.975)) 

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

################################################################################
## Create data for plots
## Samples from mcmc chains for credible intervals
################################################################################


ager=0:80
cutoff = floor(cutoffPointEst[[1]])
agey = 0:cutoff
ageo = (cutoff+1):80
numSamples = 1000

foiVector = paramVector[1:10]
foiEstimates = paramEstimates[1:10]

for(ii in 1:length(foiVector)){
  outDf <- matrix(NA,nrow=numSamples, ncol = length(ager))
  foiStudy <- foiVector[ii]
  
  
  for (kk in 1:numSamples ) {
    randomNumber <- floor(runif(1, min = 1, max = nrow(mcmcMatrix)))
    
    deltaSample <- mcmcMatrix[randomNumber,"delta"]
    lambdaYoungSample <- mcmcMatrix[randomNumber,foiStudy]
    alphaSample <- mcmcMatrix[randomNumber,"alpha"]
    
    newRowYoung <- (lambdaYoungSample / (lambdaYoungSample+deltaSample)) *(1 - exp(-agey*(lambdaYoungSample+deltaSample)))
    
    value_at_cutoff <- (lambdaYoungSample / (lambdaYoungSample+deltaSample))* (1-exp(-(lambdaYoungSample +deltaSample)*cutoff))
    
    newRowOld <- (value_at_cutoff - (lambdaYoungSample*alphaSample) / ((lambdaYoungSample*alphaSample)+deltaSample) ) * exp(-((lambdaYoungSample*alphaSample) + deltaSample)*(ageo-cutoff)) + (lambdaYoungSample*alphaSample) / ((lambdaYoungSample*alphaSample)+deltaSample)
    
    newRow <- c(newRowYoung, newRowOld)
    
    outDf[kk,] <- newRow 
    
  }
  
  # for each row in the matrix get quantiles
  quantileMatrix <- matrix(NA,nrow=ncol(outDf), ncol = 3)
  for(jj in 1:ncol(outDf)){
    quantiles <- outDf[,jj] %>% quantile(probs=c(.5,.025,.975))
    quantileMatrix[jj,] <- quantiles
  }
  
  
  lambdaYoung = foiEstimates[[ii]][1]
  delta = deltaPointEst[1]
  alpha = alphaPointEst[1]
  
  meanYoung <- (lambdaYoung / (lambdaYoung+delta)) *(1 - exp(-agey*(lambdaYoung+delta)))
  
  meanCutoff <- (lambdaYoung / (lambdaYoung+delta))* (1-exp(-(lambdaYoung +delta)*cutoff))
  
  meanOld <- (meanCutoff - (lambdaYoung*alpha) / ((lambdaYoung*alpha)+delta) ) * exp(-((lambdaYoung*alpha) + delta)*(ageo-cutoff)) + (lambdaYoung*alpha) / ((lambdaYoung*alpha)+delta)
  
  mean <- c(meanYoung, meanOld)
  
  
  # Create a dataframe for plotting
  df_upperLower = data.frame(
    midpoint = ager,
    mean = mean,
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
  ggtitle("HCoV-HKu1") +
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
  ggtitle("HCoV-NL63") +
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
  geom_line(data=df_lambdaSar, aes(x = midpoint, y=mean),color = "#A03E99") + # dark green
  geom_ribbon(data = df_lambdaSar, alpha=0.2, fill = "#A03E99")+
  geom_point(data=datSar,color="#A03E99")+
  geom_linerange(data=datSar,color="#A03E99") +
  xlab("Age (years)") + ylab("Proportion Seropositive") +
  ylim(0,1) +
  scale_x_continuous(breaks=seq(0,100,by=10)) +
  ggtitle("HCoV-OC43") +
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
  ylim(0,1) +
  scale_x_continuous(breaks=seq(0,100,by=10)) +
  ggtitle("HCoV-229E") +
  theme_bw()

plot229e

plot_grid(plot229e,hku1Plot,nl63Plot,oc43Plot)


