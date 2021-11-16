
createFig3a <- function(cutoff, cutoffUpper, cutoffLower,delta,deltaUpper,deltaLower,alpha,alphaUpper,alphaLower){
  # Range of FOI values 
  lambdaValues <- seq(0,1, by = 0.01)
  
  ##################################################
  # Reverse catalytic model with age-varying FOI
  ##################################################
  
  outDf <- data.frame(matrix(ncol = 2, nrow = 0))
  
  ager=0:80
  agey = 0:(cutoff-1)
  ageo = (cutoff):80
  
  for (i in 1:length(lambdaValues)){
    lambdaYoungSample <- lambdaValues[i]
    lambdaOldSample <- lambdaYoungSample * alpha
    
    newRowYoung <- (lambdaYoungSample / (lambdaYoungSample+delta)) * (1 - exp(-agey*(lambdaYoungSample+delta)))
    
    # Adapt equation so matches up at cutoff:
    
    value_at_cutoff <- (lambdaYoungSample / (lambdaYoungSample+delta))* (1-exp(-(lambdaYoungSample +delta)*cutoff))
    
    newRowOld <- (value_at_cutoff - lambdaOldSample / (lambdaOldSample+delta) ) * exp(-(lambdaOldSample + delta)*(ageo-cutoff)) + lambdaOldSample / (lambdaOldSample+delta)
    
    newRow <- c(newRowYoung, newRowOld)
    
    seroprev30 <- newRow[31]
    
    outDf <-rbind(outDf,c(lambdaYoungSample,seroprev30))
  }
  
  outDf <- as_tibble(outDf)
  names <- c("foi", "revCatFoi")
  colnames(outDf) <- names
  
  ## Repeat this for upper and lower CI
  outDf_upper <- data.frame(matrix(ncol = 2, nrow = 0))
  agey = 0:(cutoffUpper-1)
  ageo = (cutoffUpper):80
  
  for (i in 1:length(lambdaValues)){
    lambdaYoungSample <- lambdaValues[i]
    lambdaOldSample <- lambdaYoungSample * alphaUpper
    
    newRowYoung <- (lambdaYoungSample / (lambdaYoungSample+deltaUpper)) * (1 - exp(-agey*(lambdaYoungSample+deltaUpper)))
    
    # Adapt equation so matches up at cutoff:
    
    value_at_cutoff <- (lambdaYoungSample / (lambdaYoungSample+deltaUpper))* (1-exp(-(lambdaYoungSample +deltaUpper)*cutoffUpper))
    
    newRowOld <- (value_at_cutoff - lambdaOldSample / (lambdaOldSample+deltaUpper) ) * exp(-(lambdaOldSample + deltaUpper)*(ageo-cutoffUpper)) + lambdaOldSample / (lambdaOldSample+deltaUpper)
    
    newRow <- c(newRowYoung, newRowOld)
    
    seroprev30 <- newRow[31]
    
    outDf_upper <-rbind(outDf_upper,c(lambdaYoungSample,seroprev30))
  }
  
  
  outDf_lower <- data.frame(matrix(ncol = 2, nrow = 0))
  agey = 0:(cutoffLower-1)
  ageo = (cutoffLower):80
  
  for (i in 1:length(lambdaValues)){
    lambdaYoungSample <- lambdaValues[i]
    lambdaOldSample <- lambdaYoungSample * alphaLower
    
    newRowYoung <- (lambdaYoungSample / (lambdaYoungSample+deltaLower)) * (1 - exp(-agey*(lambdaYoungSample+deltaLower)))
    
    # Adapt equation so matches up at cutoff:
    
    value_at_cutoff <- (lambdaYoungSample / (lambdaYoungSample+deltaLower))* (1-exp(-(lambdaYoungSample +deltaLower)*cutoffLower))
    
    newRowOld <- (value_at_cutoff - lambdaOldSample / (lambdaOldSample+deltaLower) ) * exp(-(lambdaOldSample + deltaLower)*(ageo-cutoffLower)) + lambdaOldSample / (lambdaOldSample+deltaLower)
    
    newRow <- c(newRowYoung, newRowOld)
    
    seroprev30 <- newRow[31]
    
    outDf_lower <-rbind(outDf_lower,c(lambdaYoungSample,seroprev30))
  }
  
  outDf = outDf %>% 
    mutate(upper = outDf_upper$X0.1) %>% 
    mutate(lower = outDf_lower$X0.1)
  
  ##################################################
  # Catalytic model (no waning included)
  ##################################################
  
  ager=0:80
  
  outDf2 <- data.frame(matrix(ncol = 2, nrow = 0))
  
  for (i in 1:length(lambdaValues)){
    lambdaYoungSample <- lambdaValues[i]
    
    newRowYoung <- 1 - exp(-ager*lambdaYoungSample)
    
    seroprev30 <- newRowYoung[31]
    
    # newRowPrev <- c(lambdaYoungSample,seroprev30)
    
    outDf2 <-rbind(outDf2,c(lambdaYoungSample,seroprev30))
  }
  
  outDf2 <- as_tibble(outDf2)
  names <- c("foi", "cat")
  colnames(outDf2) <- names
  
  
  ##################################################
  # Reverse catalytic model 
  ##################################################
  
  outDf3 <- data.frame(matrix(ncol = 2, nrow = 0))
  
  for (i in 1:length(lambdaValues)){
    lambdaYoungSample <- lambdaValues[i]
    
    newRowYoung <- (lambdaYoungSample / (lambdaYoungSample+delta)) * (1 - exp(-ager*(lambdaYoungSample+delta)))
    
    seroprev30 <- newRowYoung[31]
    
    # newRowPrev <- c(lambdaYoungSample,seroprev30)
    
    outDf3 <-rbind(outDf3,c(lambdaYoungSample,seroprev30))
  }
  
  ## Repeat for upper and lower CI
  outDf3_upper <- data.frame(matrix(ncol = 2, nrow = 0))
  
  for (i in 1:length(lambdaValues)){
    lambdaYoungSample <- lambdaValues[i]
    
    newRowYoung <- (lambdaYoungSample / (lambdaYoungSample+deltaUpper)) * (1 - exp(-ager*(lambdaYoungSample+deltaUpper)))
    
    seroprev30 <- newRowYoung[31]
    
    # newRowPrev <- c(lambdaYoungSample,seroprev30)
    
    outDf3_upper <-rbind(outDf3_upper,c(lambdaYoungSample,seroprev30))
  }
  
  outDf3_lower <- data.frame(matrix(ncol = 2, nrow = 0))
  
  for (i in 1:length(lambdaValues)){
    lambdaYoungSample <- lambdaValues[i]
    
    newRowYoung <- (lambdaYoungSample / (lambdaYoungSample+deltaLower)) * (1 - exp(-ager*(lambdaYoungSample+deltaLower)))
    
    seroprev30 <- newRowYoung[31]
    
    # newRowPrev <- c(lambdaYoungSample,seroprev30)
    
    outDf3_lower <-rbind(outDf3_lower,c(lambdaYoungSample,seroprev30))
  }
  
  
  outDf3 <- as_tibble(outDf3)
  x <- c("foi", "revCat")
  colnames(outDf3) <- x
  
  ##################################################
  
  ## Combine models into a dataset
  outDf <- outDf %>% 
    mutate(cat = outDf2$cat) %>% 
    mutate(revCat = outDf3$revCat) %>% 
    mutate(revCatUpper = outDf3_upper$X0.1) %>% 
    mutate(revCatLower = outDf3_lower$X0.1)
  
  
  fig3a = ggplot(outDf, aes(x=foi)) +
    geom_line(aes(y=revCatFoi, colour = "Reverse catalytic & age varying FOI"), size=0.8) +
    geom_ribbon(alpha=0.4, aes(ymin = lower, ymax = upper), fill = "#52796f") +
    geom_line(aes(y=cat, colour = "Catalytic"), size=0.8)  +
    geom_line(aes(y=revCat, colour = "Reverse Catalytic"), size=0.8) +
    geom_ribbon(data = outDf, alpha=0.4, aes(ymin = revCatLower, ymax = revCatUpper), fill = "#003366")+
    scale_colour_manual(name = "", 
                        breaks = c("Catalytic", "Reverse catalytic & age varying FOI","Reverse Catalytic"),
                        values = c("#6a040f", "#52796f","#003366")) + 
    xlab("Force of Infection") + ylab("Proportion seropositive at age 30") +
    theme_bw() +
    theme(legend.position = "bottom")
  
  return(fig3a)
}