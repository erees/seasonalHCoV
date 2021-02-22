createFig3b <- function(cutoff,delta,lambda,alpha){
  
  # Set waning to be 0
  # Interested in the scenario where there is no waning occuring (In order to calculate "never infected")
  delta <- 0
  
  ager=0:80
  agey = 0:(cutoff-1)
  ageo = (floor(cutoff)):80
  
  lambdaYoung <- lambda
  lambdaOld <- lambdaYoung * alpha
  
  meanYoung <- (lambdaYoung / (lambdaYoung+delta)) *(1 - exp(-agey*(lambdaYoung+delta)))
  
  meanCutoff <- (lambdaYoung / (lambdaYoung+delta))* (1-exp(-(lambdaYoung +delta)*cutoff))
  
  meanOld <- (meanCutoff - (lambdaOld) / ((lambdaOld)+delta) ) * exp(-((lambdaOld) + delta)*(ageo-cutoff)) + (lambdaOld) / ((lambdaOld)+delta)
  
  newRow <- c(meanYoung, meanOld)
  
  neverInf <- 1-newRow
  
  
  df <- tibble(age = ager,neverInfected = neverInf, infected = newRow)
  
  dfStacked <- df %>% 
    gather("status","value",2:3,-age)
  
  ### Now start to allow infections 
  ## Would expect infections to occur 1/foi + 1/waning
  # 1/w +1/foi = time after first infection
  # 1/0.93 + 1/1.07 = 2.01
  
  infectedOnce = df$infected
  ## working in integers, so 2.01 -> 2
  infectedTwice = c(0,0,infectedOnce[1:79])
  infectedThree = c(0,0,infectedTwice[1:79])
  infectedFour = c(0,0,infectedThree[1:79])
  
  df <- df %>% 
    mutate(temp = infectedTwice) %>% 
    mutate(temp2 = infectedThree) %>% 
    mutate(infectedFour = infectedFour) %>% 
    mutate(infectedThree = temp2 - infectedFour) %>% 
    mutate(infectedTwice = temp - infectedThree - infectedFour)
  
  df <- df %>% 
    mutate(infectedOnce = infected - infectedTwice - infectedThree - infectedFour)
  
  df1 <- df %>% 
    select(-infected,-temp,-temp2) %>% 
    relocate(infectedOnce, .after = neverInfected) %>% 
    relocate(infectedTwice, .after = infectedOnce) %>% 
    relocate(infectedFour, .after = infectedThree) 
  
  dfStacked2 <- df1 %>% 
    gather("status","value",2:6,-age) %>% 
    mutate(status = factor(status,levels = c("infectedFour", "infectedThree","infectedTwice","infectedOnce","neverInfected")))
  
  plot3b <- ggplot(dfStacked2, aes(x=age, y=value, fill=status)) +
    geom_bar(position="fill", stat="identity") +
    xlim(0,15) +
    labs(fill = "Number of infections") +
    scale_fill_manual(name = "Number of infections",
                      labels = c(" ≥4", "3", "2","1","0"),
                      values = c("#212B31","#354f52","#52796f","#84a98c","#cad2c5")) +
    xlab("Age (years)") + ylab("Proportion infected") +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.direction = "horizontal")
  
  return(plot3b)
}