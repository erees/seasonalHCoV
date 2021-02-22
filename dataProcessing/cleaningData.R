## Cleaning seroprevalence data
## Removing under 1's


library(tidyverse)
library(binom)

dat <- read_csv("41467_2020_18450_MOESM7_ESM-1.csv")

######################################
### 229E
######################################

###################
### shao
###################

datShao_229e <- dat %>%
  filter(PMID == 17889596) %>%
  filter(Outcome == "HCoV-229E") %>%
  select(ageStart,ageEnd,N,N_positive) %>% 
  mutate(author = "datShao_229e")

datShao_229e <- datShao_229e %>%
  mutate(midpoint = c(0.08, 0.2085, 0.375, 0.542, 0.7085, 0.9165, 1,2,3,4,7,14.5)) %>%
  mutate(seropostive = N_positive/N) %>% 
  filter(ageStart >=1)

# Adding confidence intervals to serosamples
datShao_229e[,c("mean","lower","upper")] <- binom.confint(datShao_229e$N_positive, datShao_229e$N,method="exact")[,c("mean","lower","upper")]

###################
### zhou
###################

datZhou_229e <- dat %>%
  filter(PMID == 24040960) %>%
  filter(Outcome == "HCoV-229E") %>%
  filter(Assay == "IFA IgG") %>%
  select(ageStart,ageEnd,N,N_positive)

datZhou_229e <- datZhou_229e %>%
  mutate(midpoint = c(0.75,2,5,10.5,23,35.5,45.5,63)) %>%
  mutate(seropostive = N_positive/N) %>% 
  mutate(author = "datZhou_229e") %>% 
  filter(ageStart >=1)

# Adding confidence intervals to serosamples
datZhou_229e[,c("mean","lower","upper")] <- binom.confint(datZhou_229e$N_positive, datZhou_229e$N,method="exact")[,c("mean","lower","upper")]

###################
### cavallaro
###################

datCav <- dat %>%
  filter(PMID == 5504709) %>%
  filter(Outcome_definition == "seropositivity") %>%
  select(ageStart,ageEnd,N,N_positive)

datCav <- datCav %>%
  mutate(midpoint = c(7.5,12.5,17.5,24,34,44)) %>%
  mutate(seropostive = N_positive/N) %>% 
  mutate(author = "cav")

# Adding confidence intervals to serosamples
datCav[,c("mean","lower","upper")] <- binom.confint(datCav$N_positive, datCav$N,method="exact")[,c("mean","lower","upper")]

######################################
### HKU1
######################################

###################
### chan
###################

datChan <- dat %>%
  filter(PMID == 19342289) %>%
  filter(Outcome == "HCoV-HKU1") %>%
  filter(Assay_cutpoint == "Mean + 3SD (OD>0.495)") %>%
  select(ageStart,ageEnd,N,N_positive)

datChan <- datChan %>% 
  mutate(midpoint = c(5,15.5,25.5,35.5,45.5,55.5,65.5)) %>% 
  mutate(seropostive = N_positive/N) %>% 
  mutate(author = "chan")

# Adding confidence intervals to serosamples
datChan[,c("mean","lower","upper")] <- binom.confint(datChan$N_positive, datChan$N,method="exact")[,c("mean","lower","upper")]

###################
### zhou
###################

datZhou_hku1 <- dat %>%
  filter(PMID == 24040960) %>%
  filter(Outcome == "HCoV-HKU1") %>%
  filter(Assay == "IFA IgG") %>%
  select(ageStart,ageEnd,N,N_positive)

datZhou_hku1 <- datZhou_hku1 %>%
  mutate(midpoint = c(0.75,2,5,10.5,23,35.5,45.5,63)) %>%
  mutate(seropostive = N_positive/N) %>% 
  mutate(author = "datZhou_hku1") %>% 
  filter(ageStart >=1)

# Adding confidence intervals to serosamples
datZhou_hku1[,c("mean","lower","upper")] <- binom.confint(datZhou_hku1$N_positive, datZhou_hku1$N,method="exact")[,c("mean","lower","upper")]

######################################
### OC43
######################################

###################
### zhou
###################

datZhou_oc43 <- dat %>%
  filter(PMID == 24040960) %>%
  filter(Outcome == "HCoV-OC43") %>%
  filter(Assay == "IFA IgG") %>%
  select(ageStart,ageEnd,N,N_positive)

datZhou_oc43 <- datZhou_oc43 %>%
  mutate(midpoint = c(0.75,2,5,10.5,23,35.5,45.5,63)) %>%
  mutate(seropostive = N_positive/N) %>% 
  mutate(author = "datZhou_oc43") %>% 
  filter(ageStart >=1)

# Adding confidence intervals to serosamples
datZhou_oc43[,c("mean","lower","upper")] <- binom.confint(datZhou_oc43$N_positive, datZhou_oc43$N,method="exact")[,c("mean","lower","upper")]

###################
### Monto
###################

datMonto <- dat %>%
  filter(PMID == 4816305) %>%
  filter(TimePeriodStart == 1966) %>%
  select(ageStart,ageEnd,N,N_positive)

datMonto <- datMonto %>%
  mutate(midpoint = c(7,12,17,24.5,34.5,44.5)) %>%
  mutate(seropostive = N_positive/N) %>% 
  mutate(author = "monto66")

# Adding confidence intervals to serosamples
datMonto[,c("mean","lower","upper")] <- binom.confint(datMonto$N_positive, datMonto$N,method="exact")[,c("mean","lower","upper")]

###################
### Sarentenu
###################

datSar1 <- dat %>%
  filter(PMID == 6248465) %>%
  filter(TimePeriodStart == 1975 ) %>%
  select(ageStart,ageEnd,N,N_positive)

datSar2 <- dat %>%
  filter(PMID == 6248465) %>%
  filter(TimePeriodStart == 1976) %>%
  select(ageStart,ageEnd,N,N_positive)

newN = c()
newNp = c()

for(i in 1:4){
  newN = c(newN,datSar1[i,]$N + datSar2[i,]$N)
  newNp = c(newNp,  datSar1[i,]$N_positive + datSar2[i,]$N_positive)
}

datSar1 <- datSar1 %>% 
  mutate(N = newN) %>% 
  mutate(N_positive = newNp)

datSar <- datSar1 %>% 
  mutate(midpoint = c(7.5,19.5,42,69.5)) %>% 
  mutate(seropostive = N_positive/N) %>% 
  mutate(author = "sar75")

# Adding confidence intervals to serosamples
datSar[,c("mean","lower","upper")] <- binom.confint(datSar$N_positive, datSar$N,method="exact")[,c("mean","lower","upper")]

######################################
### NL63
######################################

###################
### shao
###################

datShao_nl63 <- dat %>%
  filter(PMID == 17889596) %>%
  filter(Outcome == "HCoV-NL63") %>%
  select(ageStart,ageEnd,N,N_positive) %>% 
  mutate(author = "datShao_nl63")  

datShao_nl63 <- datShao_nl63 %>%
  mutate(midpoint = c(0.08, 0.2085, 0.375, 0.542, 0.7085, 0.9165, 1,2,3,4,7,14.5)) %>%
  mutate(seropostive = N_positive/N) %>% 
  filter(ageStart >=1)

# Adding confidence intervals to serosamples
datShao_nl63[,c("mean","lower","upper")] <- binom.confint(datShao_nl63$N_positive, datShao_nl63$N,method="exact")[,c("mean","lower","upper")]

###################
### zhou
###################

datZhou_nl63 <- dat %>%
  filter(PMID == 24040960) %>%
  filter(Outcome == "HCoV-NL63") %>%
  filter(Assay == "IFA IgG") %>%
  select(ageStart,ageEnd,N,N_positive)

datZhou_nl63 <- datZhou_nl63 %>%
  mutate(midpoint = c(0.75,2,5,10.5,23,35.5,45.5,63)) %>%
  mutate(seropostive = N_positive/N) %>% 
  mutate(author = "datZhou_nl63") %>% 
  filter(ageStart >=1)

# Adding confidence intervals to serosamples
datZhou_nl63[,c("mean","lower","upper")] <- binom.confint(datZhou_nl63$N_positive, datZhou_nl63$N,method="exact")[,c("mean","lower","upper")]
###################
###################

datComb <- datShao_229e %>% 
  bind_rows(datZhou_229e) %>% 
  bind_rows(datCav) %>% 
  bind_rows(datChan) %>% 
  bind_rows(datZhou_hku1) %>% 
  bind_rows(datZhou_oc43) %>% 
  bind_rows(datMonto) %>% 
  bind_rows(datSar) %>% 
  bind_rows(datZhou_nl63) %>% 
  bind_rows(datShao_nl63) 

saveRDS(datComb, "cleanedCode/cleanedData.RDS")
