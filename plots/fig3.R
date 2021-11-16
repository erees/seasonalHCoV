library(tidyverse)
library(cowplot)

source("plots/fig3aFunction.R")
source("plots/fig3bFunction.R")

## Parameter values from reverse cat model with age-varying FOI (alpha, waning and cutoff jointly estimated)
cutoff = 8.49
cutoffUpper = 9.94
cutoffLower = 7.52

delta <- 0.45 #waning
deltaUpper = 0.64
deltaLower = 0.32

alpha = 1.93
alphaUpper = 2.19
alphaLower = 1.69

lambda <-  0.46 # median foi

fig3a <- createFig3a(cutoff = cutoff, 
                     cutoffUpper = cutoffUpper, 
                     cutoffLower = cutoffLower,
                     delta = delta, 
                     deltaUpper = deltaUpper,
                     deltaLower = deltaLower,
                     alpha = alpha,
                     alphaLower = alphaLower,
                     alphaUpper = alphaUpper)

fig3b <- createFig3b(cutoff = cutoff,
                     lambda = lambda,
                     alpha = alpha)


fig3 <- plot_grid(fig3a,fig3b,labels = c('A', 'B'))
