#Hypervolume analysis for native and non-native urban species 

####Set up####

#load packages
library(plyr)
library(tidyr)
library(ggplot2)
library(ggfortify)
library(ggpubr)
library(cowplot)
library(hypervolume)
library(dplyr)
library(visdat)
library(purrr)
library(tidyverse)


#Load data sets (these data sets were compiled and cleaned with "Data_harmonization script" script ) 
Diaz_Final <- read.csv("Data/Diaz_Gubic_final_tot2.csv",row=1)

#### PCA for all species####
#scale data
z_Diaz_Final<-scale(Diaz_Final[, c("Leaf_area","LMA","Leaf_N","Seed_mass","Stem_density","Height")])

#run PCA
PCA <- prcomp(z_Diaz_Final)

#look at PCA results and extract loadings and scores
summary(PCA) 
PCAvalues   <- data.frame(Species = Diaz_Final$Species_WCVP, provenance = Diaz_Final$provenance_glonaf, growth_form= Diaz_Final$growth_form, PCA$x)# Extract PC axes for plotting
PCAloadings <- data.frame(Variables = rownames(PCA$rotation), PCA$rotation) # Extract loading of the variables


#get counts for # of native and non-native species
prov_counts<-PCAvalues %>%
            group_by(provenance) %>%
            tally 
prov_counts

####Hypervolume overlap- All Species####

# extract pca values for the first 3 axes (83% of the variation)
nm_native <- filter(PCAvalues, provenance == "native") %>% 
             select(PC1,PC2,PC3)
nm_non_native <- filter(PCAvalues, provenance == "non_native")%>% 
                 select(PC1,PC2,PC3)

# Generate Hypervolumes calculate observed overlap stats
hv_native <- hypervolume(nm_native,name='native')
hv_non_native <- hypervolume(nm_non_native,name='non_native')
observed_All_M1 = hypervolume_overlap_statistics(hypervolume_set(hv_native, hv_non_native, check.memory = FALSE))


# Method 1: Null model from Chen et al. 2024 "boot-strap method"

#combine hypervolumes of native and non-native
combined_sample <- rbind(nm_native, nm_non_native)
combined_HV <- hypervolume(combined_sample)

#this takes over 12 hours to run, you can import the model outputs from an .RDS file and make the figures below
# Create bootstrapped hypervolumes of both sample sizes
#method1_path_size_1564 = hypervolume_resample("urb_1564_boot", combined_HV, "bootstrap", n = 1000)
#method1_path_size_1213 = hypervolume_resample("urb_1213_boot", combined_HV, "bootstrap", n = 1000)

#All_M1 = hypervolume_overlap_test(hv_native, hv_non_native, c(method1_path_size_1564, method1_path_size_1213))
#write.csv(All_M1$distribution, "Outputs/All_M1.csv")


#import hypervolume overlap test results(if you don't run the code above yourself)
All_M1<-readRDS("Data/All_M1.rds")

#extract p-values from test
All_M1$p_values


#plot null distribution of hypervolume overlap for supplement Fig.  S2
All_overlap_plot<-All_M1$plots$sorensen + 
                xlab("Sørensen distance") +
                ylab("Density") +
                ggtitle("All Species") +
                xlim(0.8, 1) +
                theme_bw()


####hypervolume overlap for native and non-native herb species####

#Filter PCA for herb species#
PCAvalues_herb<- PCAvalues%>%
                 filter(growth_form=="herb") 

#get counts for native and non-native herbs
prov_herb_counts<-PCAvalues_herb %>%
                  group_by(provenance) %>%
                  tally 
prov_herb_counts


# extract pca values for the first 3 axes (82% of the variation)

nm_native_herb <- filter(PCAvalues_herb, provenance == "native") %>% 
                  select(PC1,PC2,PC3)
nm_non_native_herb <- filter(PCAvalues_herb, provenance == "non_native")%>% 
                     select(PC1,PC2,PC3)

# Generate Hypervolumes and calculate overlap statistics
hv_native_herb = hypervolume(nm_native_herb,name='native')
hv_non_native_herb = hypervolume(nm_non_native_herb,name='non_native')
observed_herb_M1 = hypervolume_overlap_statistics(hypervolume_set(hv_native_herb, hv_non_native_herb, check.memory = FALSE))


# # Method 1: Null model from Chen et al. 2024 "boot-strap method"
combined_sample_herb = rbind(nm_native_herb, nm_non_native_herb)
combined_HV_herb = hypervolume(combined_sample_herb)

#this takes over 12 hours to run, you can import the model outputs from an .RDS file and make the figures below
#Create bootstrapped hypervolumes of both sample sizes
#method1_path_size_467 = hypervolume_resample("urb_467_boot", combined_HV_herb, "bootstrap", n = 1000, cores=10)
#method1_path_size_707 = hypervolume_resample("urb_707_boot", combined_HV_herb, "bootstrap", n = 1000, cores=10)

#Herb_M1 = hypervolume_overlap_test(hv_native_herb, hv_non_native_herb, c(method1_path_size_467, method1_path_size_707), cores=10)

#write.csv(Herb_M1$distribution, "Outputs/Herb_M1.csv")


#import results of hypervolume overlap test if you don't want to run
Herb_M1<-readRDS("Data/Herb_M1.rds")

#extract pvalues
Herb_M1$p_values

#plot null distribution of hypervolume overlap for supplement fig S3 
herb_overlap_plot<-Herb_M1$plots$sorensen + 
                  xlab("Sørensen distance") +
                  ylab("Density") +
                  ggtitle("Herbaceous Species") +
                  xlim(0.8, 1) +
                  theme_bw()


####Hypervolume overlap for woody/shrub species####

#filter woody/shrub pca scores
PCAvalues_woody<- PCAvalues%>%
                  filter(growth_form=="tree"|growth_form=="shrub") 

#counts of native (n=1002) and non-native(n=534) woody
prov_woody_counts<-PCAvalues_woody %>%
                  group_by(provenance) %>%
                  tally 
prov_woody_counts



# extract pca values for the first 3 axes (83% of the variation)

nm_native_wood <- filter(PCAvalues_woody, provenance == "native") %>% 
                  select(PC1,PC2,PC3)
nm_non_native_wood <- filter(PCAvalues_woody, provenance == "non_native")%>% 
                      select(PC1,PC2,PC3)

# Generate Hypervolumes
hv_native_wood = hypervolume(nm_native_wood,name='native')
hv_non_native_wood = hypervolume(nm_non_native_wood,name='non_native')

observed_woody_M1 = hypervolume_overlap_statistics(hypervolume_set(hv_native_wood, hv_non_native_wood, check.memory = FALSE))


# Method 1: Null model from Chen et al. 2024 "boot-strap method"
combined_sample = rbind(nm_native_wood, nm_non_native_wood)
combined_HV_wood = hypervolume(combined_sample)

#this takes over 12 hours to run, you can import the model outputs from an .RDS file and make the figures below
# Create bootstrapped hypervolumes of both sample sizes
#method1_path_size_1097 = hypervolume_resample("urb_1097_boot", combined_HV_wood, "bootstrap", n = 1000, cores = 10)
#method1_path_size_506 = hypervolume_resample("urb_506_boot", combined_HV_wood, "bootstrap", n = 1000, cores = 10)

#run hypervolume overlap test
#Woody_M1 = hypervolume_overlap_test(hv_native_wood, hv_non_native_wood, c(method1_path_size_1097, method1_path_size_506))

#write.csv(Woody_M1$distribution, "Outputs/Woody_M1.csv")


#import results of hypervolume overlap test if you don't want to run code above due to time
Woody_M1<-readRDS("Data/Woody_M1.rds")

#extract pvalues
Woody_M1$p_values

#plot null distribution of hypervolume overlap for supplement fig S4

wood_overlap_plot<-Woody_M1$plots$sorensen + 
                  xlab("Sørensen distance") +
                  ylab("Density") +
                  ggtitle("Woody Species") +
                  xlim(0, 1) +
                  theme_bw()


