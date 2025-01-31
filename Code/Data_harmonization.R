#Harmonizing data for urban species trait space analysis 

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



#Load and clean Diaz data from Diaz et al. 2022 with WCVP taxonomy
Diaz_WCVP <- read.csv("Data/Diaz_WCVP.csv")

#add a data source column
Diaz_WCVP$source <- 'Diaz'
Diaz_WCVP$dupes = duplicated(Diaz_WCVP$sp_WCVP_binomial)#look for duplicates

#create Diaz data frame log continuous traits)
Diaz_WCVP_clean <- data.frame(
                  Species_WCVP  = Diaz_WCVP$sp_WCVP_binomial,
                  Leaf_area     = log10(Diaz_WCVP$Leaf.Area..mm2.),
                  LMA           = log10(Diaz_WCVP$LMA..g.m2.),
                  Leaf_N        = log10(Diaz_WCVP$Nmass..mg.g.),
                  Seed_mass     = log10(Diaz_WCVP$Seed.Mass..mg.),
                  Stem_density  = log10(Diaz_WCVP$SSD.combined..mg.mm3.),
                  Height        = log10(Diaz_WCVP$Plant.Height..m. ),
                  source        = Diaz_WCVP$source,
                  growth_form   = Diaz_WCVP$growth_form,
                  duplicates    = Diaz_WCVP$dupes)%>%
                  filter(duplicates == "FALSE")%>%
                  select(!duplicates)


# urban species occurences extracted from gbif and gubic data bases(n = 694,661)
gubic_urban_sp<-readRDS("Data/gbif_gubic_list_per_city_final2_20241015.RDS")%>%
                mutate(status_glonaf = replace_na(status_glonaf, "native"))%>%
                mutate(provenance_glonaf= recode(status_glonaf, 'naturalized' = 'non_native', 
                                                 'invasive' = 'non_native', "native"="native",
                                                 "invasive, naturalized"= "non_native"))%>%
                mutate(status_gbif = replace_na(status_gbif, "native"))%>%
                mutate(provenance_gbif= recode(status_gbif, 'naturalized' = 'non_native', 
                                               'invasive' = 'non_native', "native"="native",
                                               "invasive, naturalized"= "non_native"))%>%
                mutate(status_glonaf_gbif = replace_na(status_glonaf_gbif, "native"))%>%
                mutate(provenance_glonaf_gbif= recode(status_glonaf_gbif, 'naturalized' = 'non_native', 
                                                      'invasive' = 'non_native', "native"="native",
                                                      "invasive, naturalized"= "non_native"))%>%
                filter(n_records_urban> 0)%>%
                select(accepted_taxon_name_binomial, provenance_glonaf_gbif,provenance_gbif, provenance_glonaf)




#let's look at one species
test<-gubic_urban_sp%>%
      filter(accepted_taxon_name_binomial== "Erigeron canadensis")             

#these are taxa that either are only native or only non-native in cities (gbiff/gubic data)
provenance_list_distinct<-gubic_urban_sp%>%
                          group_by((accepted_taxon_name_binomial))%>% 
                          mutate(same = +(n_distinct(provenance_glonaf) == 1)) %>% 
                          ungroup%>% 
                          filter(same==1)%>%
                          select(accepted_taxon_name_binomial, provenance_glonaf)%>%
                          unique()%>%
                          filter(provenance_glonaf=="native")#change between native and non-native


#this is the provenance list for all gubic/gbif urban species
provenance_list<-gubic_urban_sp%>%
                 select(accepted_taxon_name_binomial, provenance_glonaf)%>%
                 unique()



#count the number cities each species is native or non-native in and the diff between
prov_counts<-gubic_urban_sp%>%
            select(accepted_taxon_name_binomial, provenance_glonaf)%>%
            gather( key, value, -accepted_taxon_name_binomial) %>%
            group_by(accepted_taxon_name_binomial, key, value) %>%
            tally %>% 
            spread(value, n, fill = 0) %>% 
            mutate(diff = non_native-native)


#make non_native lists 
non_native_total<-as.data.frame(prov_counts)%>% 
                  filter(non_native>0) %>% 
                  select(accepted_taxon_name_binomial)

#global plant traits pulled from databases (try, bien, gift)
global<-read.csv("Data/global_traits.csv", row=1)%>%
        drop_na()%>%
        mutate (source= "global" )


#create global traits data frame (log continuous traits) and format like Diaz et al. 2022 traits
global_traits <- data.frame(
                  Species_WCVP  = global$accepted_wcvp_name_binomial,
                  Leaf_area     = log10(global$leaf_area_mm2 ),
                  LMA           = log10(1/(global$SLA_m2_kg/1000)),#change to lma and g/m2 like in diaz
                  Leaf_N        = log10(global$leaf_N_per_dry_mass_mg_g),
                  Seed_mass     = log10(global$diaspore_mass_mg+0.0000001),
                  Stem_density  = log10(global$wood_density_g_cm3),#g/cm3=mg/mm3 (diaz)
                  Height        = log10(global$plant_height_m),
                  source        = global$source,
                  growth_form   = global$woodiness)%>%
                  mutate(growth_form= recode(growth_form, 'woody' = 'tree', 'non-woody' = 'herb'))


####compare data sets and find shared and unique species

#species with all 6 traits unique to global traits (global traits adds 1025 species with all 6 traits)
unique_global<-as.data.frame(setdiff(global_traits$Species_WCVP,Diaz_WCVP_clean$Species_WCV))
colnames(unique_global)<-"Species_WCVP"


#get traits for species not in Diaz that are not in Diaz 
global_new<-global_traits%>%
            filter(global_traits$Species_WCVP%in% unique_global$Species_WCVP)

#find species from global traits that are in gbif/gubic data but not in Diaz (n = 809 species)
global_urban_unique<-global_new%>%
                     filter((global_new$Species_WCVP)%in% gubic_urban_sp$accepted_taxon_name_binomial)

#get Diaz species that are urban 
Diaz_urban<-Diaz_WCVP_clean%>%
            filter((Diaz_WCVP_clean$Species_WCVP)%in% gubic_urban_sp$accepted_taxon_name_binomial) 

#join urban diaz and global species from trait databases(n=2777 species)
Diaz_Gubic<-bind_rows(Diaz_urban,global_urban_unique)


##ADD PROVENANCE INFORMATION to Diaz_Gubic data set
Diaz_Gubic_native<-provenance_list_distinct%>%
                  filter(provenance_list_distinct$accepted_taxon_name_binomial%in%Diaz_Gubic$Species_WCVP) %>% 
                  mutate(provenance_glonaf= "native")

Diaz_Gubic_non_native<-non_native_total%>%
                      filter(non_native_total$accepted_taxon_name_binomial%in%Diaz_Gubic$Species_WCVP) %>% 
                      mutate(provenance_glonaf= "non_native")


Diaz_Gubic_prov<-bind_rows(Diaz_Gubic_native,Diaz_Gubic_non_native) %>% 
                rename(Species_WCVP= accepted_taxon_name_binomial)

#there are 2777 urban species in the combined diaz and gbif/gubic data (809 additional species from gbiff/gubic)
Diaz_Gubic_final_tot<-Diaz_Gubic %>% 
                    left_join(Diaz_Gubic_prov, by= "Species_WCVP")%>%
                    drop_na()


#save out cleaned data sets for all analyses
write.csv(Diaz_Gubic_final_tot,"Data/Diaz_Gubic_final_tot2.csv")







