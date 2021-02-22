#' Code used to create Figure 2. 
#' Posterior distributions of alpha 

library(tidyverse)

## Read in alpha chains
mcmcAlpha <- readRDS("plots/mcmcAlphaPost.RDS")
mcmcAlphaAll <- readRDS("plots/mcmcAlphaPostAll.RDS")

# Bind alpha chains into one dataframe
mcmcAlpha <- mcmcAlpha %>% 
  cbind(mcmcAlphaAll)

names <- c("Cavallaro and Monto","Chan","Saranteanu","Shao","Zhou","Combined")
# Assign names
colnames(mcmcAlpha) <- names

mcmcAlpha <- as_tibble(mcmcAlpha)
# Convert to factor
mcmcAlpha = mcmcAlpha %>%
  gather() %>% 
  mutate(key = as.factor(key))

# Change factor order (for purposes of plot)
mcmcAlpha$key <- fct_relevel(mcmcAlpha$key, "Combined", after = 5)

ggplot(mcmcAlpha, aes(x = value, fill = key)) +
  geom_density(alpha=0.50,size=0.3) + 
  stat_function(fun=dgamma, args=list(shape=5, rate=5),size=0.3,linetype="dashed") +
  theme_bw() +
  theme(legend.position = "bottom",legend.title = element_blank()) +
  xlab("Alpha Estimate") +
  ylab("Density") +
  scale_fill_manual(values=c("#52796f","#C9A92C","#A03E99","#075E9D", "#C05746","black")) +
  xlim(0,4.5)

