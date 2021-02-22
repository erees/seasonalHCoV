#' Reverse catalytic model 
#' With FOI changing by age and by study
#' Using alpha as the relative chang in FOI in old the older age group
#' With cutoff and alpha being allowed to vary by setting 
#' Waning estimated by strain


library(tidyverse)
library(rjags)
library(binom)
library(varhandle)
require(MCMCvis)
library(cowplot)

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
seropos_est[i] <- ifelse(age[i] < cutoffShao,

(lambdaShao_22 / (lambdaShao_22 + delta22)) * (1 - exp(-(lambdaShao_22 + delta22)*age[i])),

((lambdaShao_22 / (lambdaShao_22+delta22))* (1-exp(-(lambdaShao_22 +delta22)*cutoffShao)) 
- ((lambdaShao_22*alphaShao) / ((lambdaShao_22*alphaShao)+delta22) )) * exp(-((lambdaShao_22*alphaShao) + delta22)*(age[i]-cutoffShao)) 
+ ((lambdaShao_22*alphaShao) / ((lambdaShao_22*alphaShao)+delta22)))
}
for (i in 7:13){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data
seropos_est[i] <- ifelse(age[i] < cutoffZhou,

(lambdaZhou_22 / (lambdaZhou_22 + delta22)) * (1 - exp(-(lambdaZhou_22 + delta22)*age[i])),

((lambdaZhou_22 / (lambdaZhou_22+delta22))* (1-exp(-(lambdaZhou_22 +delta22)*cutoffZhou)) 
- ((lambdaZhou_22*alphaZhou) / ((lambdaZhou_22*alphaZhou)+delta22) )) * exp(-((lambdaZhou_22*alphaZhou) + delta22)*(age[i]-cutoffZhou)) 
+ ((lambdaZhou_22*alphaZhou) / ((lambdaZhou_22*alphaZhou)+delta22)))
}
for (i in 14:19){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data
seropos_est[i] <- ifelse(age[i] < cutoffCavMonto,

(lambdaCav / (lambdaCav + delta22)) * (1 - exp(-(lambdaCav + delta22)*age[i])),

((lambdaCav / (lambdaCav+delta22))* (1-exp(-(lambdaCav +delta22)*cutoffCavMonto)) 
- ((lambdaCav*alphaCavMonto) / ((lambdaCav*alphaCavMonto)+delta22) )) * exp(-((lambdaCav*alphaCavMonto) + delta22)*(age[i]-cutoffCavMonto)) 
+ ((lambdaCav*alphaCavMonto) / ((lambdaCav*alphaCavMonto)+delta22)))
}
for (i in 20:26){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data
seropos_est[i] <- ifelse(age[i] < cutoffChan,

(lambdaChan / (lambdaChan + deltaHk)) * (1 - exp(-(lambdaChan + deltaHk)*age[i])),

((lambdaChan / (lambdaChan+deltaHk))* (1-exp(-(lambdaChan +deltaHk)*cutoffChan)) 
- ((lambdaChan*alphaChan) / ((lambdaChan*alphaChan)+deltaHk) )) * exp(-((lambdaChan*alphaChan) + deltaHk)*(age[i]-cutoffChan)) 
+ ((lambdaChan*alphaChan) / ((lambdaChan*alphaChan)+deltaHk)))
}
for (i in 27:33){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data
seropos_est[i] <- ifelse(age[i] < cutoffZhou,

(lambdaZhou_hk / (lambdaZhou_hk + deltaHk)) * (1 - exp(-(lambdaZhou_hk + deltaHk)*age[i])),

((lambdaZhou_hk / (lambdaZhou_hk+deltaHk))* (1-exp(-(lambdaZhou_hk +deltaHk)*cutoffZhou)) 
- ((lambdaZhou_hk*alphaZhou) / ((lambdaZhou_hk*alphaZhou)+deltaHk) )) * exp(-((lambdaZhou_hk*alphaZhou) + deltaHk)*(age[i]-cutoffZhou)) 
+ ((lambdaZhou_hk*alphaZhou) / ((lambdaZhou_hk*alphaZhou)+deltaHk)))
}
for (i in 34:40){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data

seropos_est[i] <- ifelse(age[i] < cutoffZhou,

(lambdaZhou_oc / (lambdaZhou_oc + deltaOc)) * (1 - exp(-(lambdaZhou_oc + deltaOc)*age[i])),

((lambdaZhou_oc / (lambdaZhou_oc+deltaOc))* (1-exp(-(lambdaZhou_oc +deltaOc)*cutoffZhou)) 
- ((lambdaZhou_oc*alphaZhou) / ((lambdaZhou_oc*alphaZhou)+deltaOc) )) * exp(-((lambdaZhou_oc*alphaZhou) + deltaOc)*(age[i]-cutoffZhou)) 
+ ((lambdaZhou_oc*alphaZhou) / ((lambdaZhou_oc*alphaZhou)+deltaOc)))
}
for (i in 41:46){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data

seropos_est[i] <- ifelse(age[i] < cutoffCavMonto,

(lambdaMonto / (lambdaMonto + deltaOc)) * (1 - exp(-(lambdaMonto + deltaOc)*age[i])),

((lambdaMonto / (lambdaMonto+deltaOc))* (1-exp(-(lambdaMonto +deltaOc)*cutoffCavMonto)) 
- ((lambdaMonto*alphaCavMonto) / ((lambdaMonto*alphaCavMonto)+deltaOc) )) * exp(-((lambdaMonto*alphaCavMonto) + deltaOc)*(age[i]-cutoffCavMonto)) 
+ ((lambdaMonto*alphaCavMonto) / ((lambdaMonto*alphaCavMonto)+deltaOc)))
}
for (i in 47:50){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data

seropos_est[i] <- ifelse(age[i] < cutoffSar,

(lambdaSar / (lambdaSar + deltaOc)) * (1 - exp(-(lambdaSar + deltaOc)*age[i])),

((lambdaSar / (lambdaSar+deltaOc))* (1-exp(-(lambdaSar +deltaOc)*cutoffSar)) 
- ((lambdaSar*alphaSar) / ((lambdaSar*alphaSar)+deltaOc) )) * exp(-((lambdaSar*alphaSar) + deltaOc)*(age[i]-cutoffSar)) 
+ ((lambdaSar*alphaSar) / ((lambdaSar*alphaSar)+deltaOc)))
}
for (i in 51:57){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data

seropos_est[i] <- ifelse(age[i] < cutoffZhou,

(lambdaZhou_nl / (lambdaZhou_nl + deltaNl)) * (1 - exp(-(lambdaZhou_nl + deltaNl)*age[i])),

((lambdaZhou_nl / (lambdaZhou_nl+deltaNl))* (1-exp(-(lambdaZhou_nl +deltaNl)*cutoffZhou)) 
- ((lambdaZhou_nl*alphaZhou) / ((lambdaZhou_nl*alphaZhou)+deltaNl) )) * exp(-((lambdaZhou_nl*alphaZhou) + deltaNl)*(age[i]-cutoffZhou)) 
+ ((lambdaZhou_nl*alphaZhou) / ((lambdaZhou_nl*alphaZhou)+deltaNl)))
}
for (i in 58:63){ 
n.pos[i] ~ dbinom(seropos_est[i],N[i]) #fit to binomial data

seropos_est[i] <- ifelse(age[i] < cutoffShao,

(lambdaShao_nl / (lambdaShao_nl + deltaNl)) * (1 - exp(-(lambdaShao_nl + deltaNl)*age[i])),

((lambdaShao_nl / (lambdaShao_nl+deltaNl))* (1-exp(-(lambdaShao_nl +deltaNl)*cutoffShao)) 
- ((lambdaShao_nl*alphaShao) / ((lambdaShao_nl*alphaShao)+deltaNl) )) * exp(-((lambdaShao_nl*alphaShao) + deltaNl)*(age[i]-cutoffShao)) 
+ ((lambdaShao_nl*alphaShao) / ((lambdaShao_nl*alphaShao)+deltaNl)))
}


## priors

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

alphaShao ~ dgamma(5,5)
alphaZhou ~ dgamma(5,5)
alphaCavMonto ~ dgamma(5,5)
alphaChan ~ dgamma(5,5)
alphaSar ~ dgamma(5,5)

cutoffShao ~ dunif(0,20)
cutoffZhou ~ dunif(0,20)
cutoffCavMonto ~ dunif(0,20)
cutoffChan ~ dunif(0,20)
cutoffSar ~ dunif(0,20)

delta22 ~ dunif(0,5) #uniformative prior
deltaHk ~ dunif(0,5) #uniformative prior
deltaNl ~ dunif(0,5) #uniformative prior
deltaOc ~ dunif(0,5) #uniformative prior

}"

################################################################################
# Running model
################################################################################

# number of itterations
mcmc.length=100000
## Specify my data
jdat <- list(n.pos= datComb$N_positive,
             N=datComb$N,
             age=datComb$midpoint)

jmod = jags.model(textConnection(jcode), data=jdat, n.chains=4, n.adapt=30000)
update(jmod)

jpos = coda.samples(jmod, c("lambdaShao_22", 
                            "lambdaZhou_22", 
                            "lambdaCav", 
                            "lambdaChan", 
                            "lambdaZhou_hk", 
                            "lambdaZhou_oc", 
                            "lambdaMonto", 
                            "lambdaSar", 
                            "lambdaZhou_nl",
                            "lambdaShao_nl", 
                            "alphaShao", 
                            "alphaZhou",
                            "alphaCavMonto",
                            "alphaChan",
                            "alphaSar",
                            "delta22",
                            "deltaHk",
                            "deltaOc",
                            "deltaNl",
                            "cutoffShao",
                            "cutoffZhou",
                            "cutoffCavMonto",
                            "cutoffChan",
                            "cutoffSar"), n.iter=mcmc.length)

plot(jpos) ## Check convergence of all chains

summary(jpos)

MCMCsummary(jpos, round = 2) ## Check ESS and Rhat

mcmcMatrix <- as.matrix(jpos)

# Plotting posterior distributions of all parameters
mcmcDF <- as_tibble(mcmcMatrix)
mcmcDF %>%
  gather() %>% 
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_density()

## Export posteriors for alpha (used for Fig. 2)
mcmcAlpha = mcmcDF[,1:5]
saveRDS(mcmcAlpha, "mcmcAlphaPost.RDS")

## Export posteriors for foi
mcmcFoi = mcmcDF[,12:21]
saveRDS(mcmcFoi, "mcmcFoiPost.RDS")

# Calculate DIC
dic.samples(jmod, n.iter = mcmc.length)

################################################################################
## Create point estimates for all parameters
################################################################################

delta22PointEst <- mcmcMatrix[,"delta22"] %>% quantile(probs=c(.5,.025,.975))
deltaOcPointEst <- mcmcMatrix[,"deltaOc"] %>% quantile(probs=c(.5,.025,.975))
deltaHkPointEst <- mcmcMatrix[,"deltaHk"] %>% quantile(probs=c(.5,.025,.975))
deltaNlPointEst <- mcmcMatrix[,"deltaNl"] %>% quantile(probs=c(.5,.025,.975))

alphaPointEstShao22 <- mcmcMatrix[,"alphaShao"] %>% quantile(probs=c(.5,.025,.975)) 
alphaPointEstZhou22 <- mcmcMatrix[,"alphaZhou"] %>% quantile(probs=c(.5,.025,.975))
alphaPointEstCav <- mcmcMatrix[,"alphaCavMonto"] %>% quantile(probs=c(.5,.025,.975))
alphaPointEstChan <- mcmcMatrix[,"alphaChan"] %>% quantile(probs=c(.5,.025,.975))
alphaPointEstZhouHk <- mcmcMatrix[,"alphaZhou"] %>% quantile(probs=c(.5,.025,.975))
alphaPointEstZhouOc <- mcmcMatrix[,"alphaZhou"] %>% quantile(probs=c(.5,.025,.975))
alphaPointEstMonto <- mcmcMatrix[,"alphaCavMonto"] %>% quantile(probs=c(.5,.025,.975))
alphaPointEstSar <- mcmcMatrix[,"alphaSar"] %>% quantile(probs=c(.5,.025,.975))
alphaPointEstZhouNl <- mcmcMatrix[,"alphaZhou"] %>% quantile(probs=c(.5,.025,.975))
alphaPointEstShaoNl <- mcmcMatrix[,"alphaShao"] %>% quantile(probs=c(.5,.025,.975))

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

cutoffPointEstShao22 <- mcmcMatrix[,"cutoffShao"] %>% quantile(probs=c(.5,.025,.975))
cutoffPointEstZhou22 <- mcmcMatrix[,"cutoffZhou"] %>% quantile(probs=c(.5,.025,.975))
cutoffPointEstCav <- mcmcMatrix[,"cutoffCavMonto"] %>% quantile(probs=c(.5,.025,.975))
cutoffPointEstChan <- mcmcMatrix[,"cutoffChan"] %>% quantile(probs=c(.5,.025,.975))
cutoffPointEstZhouhk <- mcmcMatrix[,"cutoffZhou"] %>% quantile(probs=c(.5,.025,.975))
cutoffPointEstZhouoc <- mcmcMatrix[,"cutoffZhou"] %>% quantile(probs=c(.5,.025,.975))
cutoffPointEstMonto <- mcmcMatrix[,"cutoffCavMonto"] %>% quantile(probs=c(.5,.025,.975))
cutoffPointEstSar <- mcmcMatrix[,"cutoffSar"] %>% quantile(probs=c(.5,.025,.975))
cutoffPointEstZhounl <- mcmcMatrix[,"cutoffZhou"] %>% quantile(probs=c(.5,.025,.975))
cutoffPointEstShaonl <- mcmcMatrix[,"cutoffShao"] %>% quantile(probs=c(.5,.025,.975))

print(1/delta22PointEst) %>% round(2)
print(1/deltaOcPointEst) %>% round(2)

paramVector = c("lambdaShao_22", 
                "lambdaZhou_22", 
                "lambdaCav", 
                "lambdaChan", 
                "lambdaZhou_hk", 
                "lambdaZhou_oc", 
                "lambdaMonto", 
                "lambdaSar", 
                "lambdaZhou_nl", 
                "lambdaShao_nl",
                "delta22",
                "deltaOc",
                "deltaHk",
                "deltaNl",
                "alphaShao22", 
                "alphaZhou22",
                "alphaCav",
                "alphaChan",
                "alphaZhouHk",
                "alphaZhouOc",
                "alphaMonto",
                "alphaSar",
                "alphaZhouNl",
                "alphaShaoNl",
                "cutoffShao22",
                "cutoffZhou22",
                "cutoffCav",
                "cutoffChan",
                "cutoffZhouHk",
                "cutoffZhouOc",
                "cutoffMonto",
                "cutoffSar",
                "cutoffZhouNl",
                "cutoffShaoNl")

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
                      delta22PointEst,
                      deltaOcPointEst,
                      deltaHkPointEst,
                      deltaNlPointEst,
                      alphaPointEstShao22,
                      alphaPointEstZhou22,
                      alphaPointEstCav,
                      alphaPointEstChan,
                      alphaPointEstZhouHk,
                      alphaPointEstZhouOc,
                      alphaPointEstMonto,
                      alphaPointEstSar,
                      alphaPointEstZhouNl,
                      alphaPointEstShaoNl,
                      cutoffPointEstShao22,
                      cutoffPointEstZhou22,
                      cutoffPointEstCav,
                      cutoffPointEstChan,
                      cutoffPointEstZhouhk,
                      cutoffPointEstZhouoc,
                      cutoffPointEstMonto,
                      cutoffPointEstSar,
                      cutoffPointEstZhounl,
                      cutoffPointEstShaonl)

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

alphaVector = c("alphaShao","alphaZhou","alphaCavMonto","alphaChan","alphaZhou","alphaZhou","alphaCavMonto","alphaSar","alphaZhou","alphaShao")
alphaEstimates = paramEstimates[12:21]

cutoffEstimates = paramEstimates[22:31]

for(ii in 1:length(foiVector)){
  outDf <- matrix(NA,nrow=numSamples, ncol = length(ager))
  foiStudy <- foiVector[ii]
  alphaStudy <- alphaVector[ii]
  cutoffStudy <- cutoffEstimates[[ii]][1]
  
  cutoff = floor(cutoffStudy)
  agey = 0:cutoff
  ageo = (cutoff+1):80
  
  for (kk in 1:numSamples ) {
    randomNumber <- floor(runif(1, min = 1, max = nrow(mcmcMatrix)))
    
    deltaSample <- mcmcMatrix[randomNumber,"delta"]
    lambdaYoungSample <- mcmcMatrix[randomNumber,foiStudy]
    alphaSample <- mcmcMatrix[randomNumber,alphaStudy]
    
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
  alpha = alphaEstimates[[ii]][1]
  
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

############################################################
## Plots
############################################################

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
  scale_x_continuous(breaks=seq(0,100,by=10)) +
  ylim(0,1) +
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
  scale_x_continuous(breaks=seq(0,100,by=10)) +
  ylim(0,1) +
  ggtitle("HCoV-229E") +
  theme_bw()

plot229e

plot_grid(plot229e,hku1Plot,nl63Plot,oc43Plot)

