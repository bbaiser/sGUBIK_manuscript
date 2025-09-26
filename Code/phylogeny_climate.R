#Get phylogeny and average precip and temp for the 2777 urban species

library(tidyverse)
options(repos = c(rtrees = 'https://daijiang.r-universe.dev', CRAN = 'https://cloud.r-project.org'))
install.packages("rtrees")
library(rtrees)
library(tidyverse)
library(sf)
library(geodata)

#Data
d_pca = read_csv("Data/Trait_pca_scores.csv")#PCA scores for the 2777 urban species



# phylogeny ====

sp_df = sp_list_df(sp_list = d_pca$Species, taxon = "plant")

filter(sp_df, 
       is.na(family))



# went through these species one by one

d_pca = mutate(d_pca, sp 
               = Species,
               
               sp = dplyr::recode(sp, 
                           "Pseudalbizzia niopoides" = 
                             "Albizia niopoides",
                           
                           "Cynanchica aristata" 
                           = "Asperula aristata",
                           
                           "Cynanchica pyrenaica" 
                           = "Asperula pyrenaica",
                           
                           "Thliphthisa purpurea" 
                           = "Asperula purpurea",
                           
                           "Robrichia schomburgkii" 
                           = "Enterolobium schomburgkii",
                           
                           "Heptapleurum actinophyllum" 
                           = "Schefflera actinophylla",
                           
                           "Heptapleurum heptaphyllum" 
                           = "Schefflera heptaphyllum",
                           
                           "Rabelera holostea" 
                           = "Stellaria holostea",
                           
                           "Jarava ichu" 
                           = "Stipa ichu",
                           
                           "Pectinopitys ferruginea" 
                           = "Prumnopitys ferruginea",
                           
                           "Phyllogeiton zeyheri" 
                           = "Berchemia zeyheri"
                           
               ),
               
               tip = str_replace_all(Species, 
                                     " ", "_")
               
)

sp_df = sp_list_df(sp_list 
                   = d_pca$sp, taxon 
                   = "plant")

filter(sp_df, 
       is.na(family)) 
# now they all have the genus in the classification



if(!file.exists("Data/phy.tre")){
  
  phy = get_tree(sp_list 
                 = sp_df, taxon = 
                   "plant")

# update tip labels back to the original names
  
  phy$tip.label = tibble(tip = phy$tip.label)|> 
    
    mutate(sp = str_replace_all(tip, 
                                "_", " ")) |> 
    
    left_join(select(d_pca, sp, Species)) |>
    
    mutate(tip2 = str_replace_all(Species,
                                  " ", 
                                  "_")) |> 
    
    pull(tip2)
  
  
  setdiff(d_pca$tip, phy$tip.label)
  
  
  plot(phy, show.tip.label 
       = F)
  
  
  ape::write.tree(phy, "Data/phy.tre")
  
} else {
  
  phy = ape::read.tree("Data/phy.tre")
  
}





# cities ====

d2 = readRDS("Data/gbif_gubic_list_per_city_final2_20241015.rds")

setdiff(d_pca$Species,
        unique(d2$accepted_taxon_name_binomial))

d2 = filter(d2, accepted_taxon_name_binomial %in% d_pca$Species)

setdiff(unique(d2$accepted_taxon_name_binomial), d_pca$Species)

n_distinct(d2$ID_HDC_G0_2) 
# 554



urban_centers = st_read("Data/urban_centres/")

filter(urban_centers, ID_HDC_ == 1226)

urban_centers2 = 
  filter(urban_centers, ID_HDC_ %in%
           unique(d2$ID_HDC_G0_2))|> 
  
  select(ID_HDC_G0_2 = ID_HDC_) |> distinct()

setdiff(unique(d2$ID_HDC_G0_2), urban_centers2$ID_HDC_G0_2)



d2 = filter(d2, ID_HDC_G0_2 != 1226)



# get worldclim data ====

# monthly data from 1970-2000, averaged

tave = geodata::worldclim_global(var= "tavg", res= 2.5, path= "Data/worldclim")

prec = geodata::worldclim_global(var= "prec", res= 2.5, path= "Data/worldclim")



tave2 = terra::app(tave, fun = mean, na.rm = T)

prec2 = terra::app(prec, fun = sum, na.rm = T)



tave_urban = terra::extract(tave2, urban_centers2, fun= mean, na.rm = T)

summary(tave_urban$mean)

prec_urban = terra::extract(prec2, urban_centers2, fun= mean, na.rm = T)


summary(prec_urban$sum)


urban_env = st_centroid(urban_centers2)|> 
  
  st_coordinates() |> 
  
  as_tibble() |> 
  
  mutate(ID_HDC_G0_2 = urban_centers2$ID_HDC_G0_2,
         
         tmean = tave_urban$mean, precip= prec_urban$sum)

plot(abs(urban_env$Y), urban_env$tmean)

plot(urban_env$tmean, urban_env$precip)



# for each species, find the city it appeared, and take the average

for(i in 1:nrow(d_pca)){
  
  sp_i = d_pca$Species[i]
  
  d2_i = filter(d2, accepted_taxon_name_binomial == sp_i)
  
  urban_env_i = filter(urban_env, ID_HDC_G0_2 %in%
             unique(d2_i$ID_HDC_G0_2))|> 
    
    mutate(lat_abs = abs(Y)) |> 
    
    summarise(ave_lat_abs = mean(lat_abs), 
              
              ave_tmean = mean(tmean),
              
              ave_precip = mean(precip))
  
  d_pca$ave_lat_abs[i]  = urban_env_i$ave_lat_abs
  
  d_pca$ave_tmean[i] = urban_env_i$ave_tmean
  
  d_pca$ave_precip[i] = urban_env_i$ave_precip
  
}

d_pca2<-d_pca

write_csv(d_pca, "Data/d_pca2.csv")
