#PCA AND KS tests of pca score distributions for All, native, and non-native urban plants
#figure 2

####Set up####

#load packages
library(plyr)
library(tidyr)
library(ggplot2)
library(ggfortify)
library(ggpubr)
library(devtools)
library(cowplot)
library(hypervolume)
library(dplyr)
library(visdat)
library(purrr)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(kSamples)


#Load data sets 

#Data from Diaz et al. 2022 and GBIF for urban for speices that have all 6 diaz traits. taxonomically matched to WCVP 
Diaz_Final <- read.csv("Data/Diaz_Gubic_final_tot2.csv",row=1)


#counts of native (n=1564) and non-native(n=1213)
prov_counts<-Diaz_Final %>%
  group_by(provenance_glonaf) %>%
  tally 
prov_counts

#### PCA for all species####
#scale data
z_Diaz_Final<-scale(Diaz_Final[, c("Leaf_area","LMA","Leaf_N","Seed_mass","Stem_density","Height")])

#run PCA
PCA <- prcomp(z_Diaz_Final)

#look at PCA results and extract loadings and scores
summary(PCA) 
PCAvalues   <- data.frame(Species = Diaz_Final$Species_WCVP, provenance = Diaz_Final$provenance_glonaf, growth_form= Diaz_Final$growth_form, PCA$x)# Extract PC axes for plotting
PCAloadings <- data.frame(Variables = rownames(PCA$rotation), PCA$rotation) # Extract loading of the variables

PCAvalues$PC2 <- (PCAvalues$PC2*-1)

pca_scores<-PCAvalues%>%
  mutate(growth_form = case_when(
  growth_form == "herb" ~ "herb",
  growth_form %in% c("shrub", "tree") ~ "woody"
))

write.csv(pca_scores, "Trait_pca_scores.csv")

