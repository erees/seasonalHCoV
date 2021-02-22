# seasonalHCoV

This repository contains the data and code used to produce the results presented in "Estimating the duration of seropositivity of human seasonal coronaviruses using seroprevalence studies".

The raw data used for all the different studies is located in `dataProcessing\41467_2020_18450_MOESM7_ESM-1.csv`. This contains estimates extracted from all the studies.

Two different scripts are used the clean the data:

* `cleaningData.R` (does not include individuals <1 year)
* `cleanrDataAllAges.R` (includes individuals <1 year)

The main model presented in the paper is `mainModel.R`.

Several sensitivity analyses were conducted, and these are located in the folder "sensitvityAnalysis":

* `mainModel_allAges.R` (the same as the main model but includes individuals ages <1 year)
* `mainModel_alphaHeld.R` (Alpha and cutoff are jointly estimated across settings, instead of for all studies)
* `mainModel_strain.R` (Waning is estimated for each strain, instead of jointly estimated across all studies and strains)
* `mainModel_widePrior.R` (the same as the main model, but a less informed prior is used for the FOI)
* `reverseCatalyticModel.R` (Simple reverse catalytic model, with no age-varying FOI).
* SensitivityAnalysisTwoStrains (Folder containing scripts for the two strain analysis where only half the data is used.

There is also folder containing code used to create Figure 2 and Figure 3 in "plots". 
